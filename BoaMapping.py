"""
NAM: BoaMapping.py (module)
DES: contains the BoAMapping and Image classes
"""
__version__=  '$Revision: 2811 $'
__date__=     '$Date: 2015-03-04 15:32:01 +0100 (Wed, 04 Mar 2015) $'

#----------------------------------------------------------------------------------
#----- Import ---------------------------------------------------------------------
#----------------------------------------------------------------------------------
from Numeric import *
from LinearAlgebra import *
import cPickle,string

from boa           import BoaDataAnalyser, BoaMessageHandler, BoaFits, BoaFlagHandler
from boa           import BoaConfig, BoaError
from boa.Bogli     import Plot, MultiPlot
from boa.fortran   import fUtilities, fMap, fStat
from boa.Utilities import Timing, compressNan, tolist_boa, inPolygon, outPolygon
from boa.Utilities import modelBaseEllipticalGaussian,fitBaseEllipticalGaussian

import os.path

#----------------------------------------------------------------------------------
#----- Image Class ----------------------------------------------------------------
#----------------------------------------------------------------------------------

class Image:
    """
    NAM: Image (class)
    DES: An object of this class describes an image and its axis
    """


    # TODO add method getPixel putPixel and/or getRegion setRegion
    def __init__(self):

        self.__MessHand=BoaMessageHandler.MessHand(self.__module__)
        
        self.Header = {
            'Object':   'Unknown',\
            'Telescope':  'None',\
            'FeBe':       'None'\
            }

        self.Data     = []   # The bitmap map (2D arrays)
	self.Weight   = []   # The weights per pixel
	self.Coverage = []   # The coverage per pixel, i.e. weight=1
        
        self.RmsBeam  = 0.   # rms/beam in current unit
	self.BeamSize = 0.   # beam size, in same unit as WCS
        
        self.WCS = {
            # Number of Axis (2 for simple images)
            'NAXIS': 2,\
            \
            # time of observation
            'MJD-OBS': 'None',\
            # frame of reference
            'RADESYS':'FK5',\
            # coordinate epoch
            'EQUINOX': 2000.0,\
            \
            # #of pixel
            'NAXIS1': 0,\
            # Type 
            'CTYPE1': 'X',\
            # pixel scale
            'CDELT1': 1.0,\
            # unit
            'CUNIT1': 'dummy',\
            # Pixel coordinate of reference point
            'CRPIX1': 0.0,\
            # celestial coordinate of reference point
            'CRVAL1': 0.0,\
            \
            'NAXIS2': 0,\
            'CTYPE2': 'Y',\
            'CDELT2': 1.0,\
            'CUNIT2': 'dummy',\
            'CRPIX2': 0.0,\
            'CRVAL2': 0.0,\
            \
            # linear transformation matrix
            'PC1_1': 1,\
            'PC1_2': 0,\
            'PC2_1': 0,\
            'PC2_2': 1}

    def __str__(self):
        """
        DES:Defines a string, shown when the print instruction is used.
        """
        out = "Image of %s with %s on %s (%ix%i pixels)"%\
              (self.Header['Object'],\
               self.Header['FeBe'],\
               self.Header['Telescope'],\
               self.WCS['NAXIS1'],self.WCS['NAXIS2'])

        return out
            
    def wcs2pix(self,X,Y):
        """
        DES: Convert from physical coordinates described by self.WCS
             to pixel coordinates
        INP: float (X,Y) : the physical coordinates to convert from
        OUT: float (i,j) : the pixel coordinates

	     We should switch to libwcs at some point
        """
        if self.WCS['CUNIT1'] == 'arcsec':
            cdeltUnit = 1./3600.
        else:
            cdeltUnit = 1.
        AXIS1 = array([self.WCS['NAXIS1'],self.WCS['CRPIX1'],self.WCS['CDELT1'],
                       self.WCS['CRVAL1'],cdeltUnit])
        AXIS2 = array([self.WCS['NAXIS2'],self.WCS['CRPIX2'],self.WCS['CDELT2'],
                       self.WCS['CRVAL2'],cdeltUnit])
        i,j = fMap.wcs2pix(X,Y,AXIS1,AXIS2)
        if not rank(X) and not rank(Y):
            i = i[0]
            j = j[0]
        return i,j
    
    def wcs2phy(self,i,j):
        """
        DES: Convert from pixel coordinates to physical (world) coordinates
        INP: float (i,j) : the pixel coordinates to convert from
        OUT: float (X,Y) : the physical coordinates

	     We should switch to libwcs at some point
        """
        if self.WCS['CUNIT1'] == 'arcsec':
            cdeltUnit = 1./3600.
        else:
            cdeltUnit = 1.
	AXIS1 = array([self.WCS['NAXIS1'],self.WCS['CRPIX1'],self.WCS['CDELT1'],
                       self.WCS['CRVAL1'],cdeltUnit])
	AXIS2 = array([self.WCS['NAXIS2'],self.WCS['CRPIX2'],self.WCS['CDELT2'],
                       self.WCS['CRVAL2'],cdeltUnit])
        X,Y = fMap.wcs2phy(i,j,AXIS1,AXIS2)
        if not rank(i) and not rank(j):
            X = X[0]
            Y = Y[0]
        return X,Y
    
   
    def computeWCS(self,pixelSize,sizeX=[],sizeY=[],minmax=[]):
        """
        DES: fill main WCS keywords according to pixel size and map limits
        INP: (int) pixelSize = size of pixel in acrsecond
             (float) sizeX = map limits in azimuth, in arcsecond
             (float) sizeY = map limits in elevation, in arcsecond
             (float) minmax = [minAzoff,maxAzoff,minEloff,maxEloff] in this order
        """

        # TODO add shift for 'CRPIX1/2' or 'CRVAL1/2' to be general and allow for
        # offset in  maps

        # determine coordinate limits and map size
        if sizeX == []:
            minAzoff = minmax[0]
            maxAzoff = minmax[1]
        else:
            minAzoff = sizeX[0]
            maxAzoff = sizeX[1]
            
        if sizeY == []:
            minEloff = minmax[2]
            maxEloff = minmax[3]
        else:
            minEloff = sizeY[0]
            maxEloff = sizeY[1]
           
        Xcenter = (minAzoff + maxAzoff)/2.
        Ycenter = (minEloff + maxEloff)/2.
        if self.WCS['CUNIT1'] == 'arcsec':
            cosY = cos(Ycenter*pi/180./3600.)
        else:
            cosY = cos(Ycenter*pi/180.)
            
        dimX = int(ceil(abs(maxAzoff - minAzoff)/pixelSize*cosY +1.))
        dimY = int(ceil(abs(maxEloff - minEloff)/pixelSize +1.))
        
        # Size of the image
        self.WCS['NAXIS1'] = dimX
        self.WCS['NAXIS2'] = dimY
        
        # Size of pixels
        self.WCS['CDELT1'] = (maxAzoff - minAzoff)/abs(maxAzoff - minAzoff)*pixelSize
        self.WCS['CDELT2'] = (maxEloff - minEloff)/abs(maxEloff - minEloff)*pixelSize
        
        # put the reference in the center of the image ...
        self.WCS['CRPIX1'] = (dimX-1)/2.
        self.WCS['CRPIX2'] = (dimY-1)/2.
        
        # ... corresponding to the center of the map
        self.WCS['CRVAL1'] = Xcenter
        self.WCS['CRVAL2'] = Ycenter

    def physicalCoordinates(self):
        """
        DES: return arrays with physical units corresponding to the map
        """
        nbX,nbY = self.WCS['NAXIS1'],self.WCS['NAXIS2']
        resultAz = zeros((nbX,nbY),'f')
        resultEl = zeros((nbX,nbY),'f')
        allY = arange(nbY)
        for i in range(nbX):
            az,el = self.wcs2phy(repeat([i],nbY),allY)
            resultAz[i,::] = az
            resultEl[i,::] = el
        return (resultAz,resultEl)


    def display(self, weight=0, coverage=0, style='idl4', caption='', wedge=1, \
		aspect=0, overplot=0, doContour=0, levels=[], labelContour=0,\
		limitsX = [], limitsY = [], limitsZ = [], \
                showRms=0, rmsKappa=3.5, noerase=0, snmap=0, cell=15,sparse=8):
        """
        DES: show the reconstructed maps in (Az,El)
        INP: (boolean) weigth,coverage : plot the rms or weight map instead of signal map
             (string)  style           : the style used for the color (default idl4)

             (string)  caption       : the caption of the plot (default '')
             (flt array) limitsX/Y/Z : the limits in X/Y/intensity
             (boolean) wedge         : draw a wedge ? (default : yes)
             (boolean) aspect        : keep the aspect ratio (default : yes)
             (boolean) overplot      : should we overplot this image (default : no)
             (boolean) doContour     : draw contour instead of map (default : no)
             (float array) levels    : the levels of the contours
                                       (default : intensity progression)
             (boolean) labelContour  : label the contour (default : no)
             (boolean) showRms       : compute and display rms/beam? (def: no)
             (boolean) noerase       : do not clear the window? (def: false)
             (boolean) snmap         : display a signal-to-noise map in arb. units (def: no)
        """
        
        if snmap:
            toPlot=self.computeSNMap(cell=cell,sparse=sparse)
        else:
            toPlot = self.Data

        WCS    = self.WCS
        
	if weight:
            toPlot = self.Weight
	if coverage:
            toPlot = self.Coverage
        
        labelX=''
        if WCS['CTYPE1'][0:4] == 'OLON':
            labelX = "\gD Az"
        elif WCS['CTYPE1'][0:4] == 'RA--':
            labelX="\ga"
            if WCS['EQUINOX'] == 2000.0:
                labelX += ' (J2000)'
            if WCS['EQUINOX'] == 1950.0:
                labelX += ' (B1950)'
        elif WCS['CTYPE1'][0:4] == 'GLON':
            labelX="Gal. long."
        else:
            labelX = WCS['CTYPE1']
        if WCS['CUNIT1']:
            labelX += " ["+WCS['CUNIT1']+"]"

        labelY=''
        if WCS['CTYPE2'][0:4] == 'OLAT':
            labelY="\gD El"
        elif WCS['CTYPE2'][0:4] == 'DEC-':
            labelY="\gd"
            if WCS['EQUINOX'] == 2000.0:
                labelY+=' (J2000)'
            if WCS['EQUINOX'] == 1950.0:
                labelY+=' (B1950)'
        elif WCS['CTYPE2'][0:4] == 'GLAT':
            labelY="Gal. lat."
        else:
            labelY = WCS['CTYPE2']
        if WCS['CUNIT2']:
            labelY += " ["+WCS['CUNIT2']+"]"

        Plot.draw(toPlot,wedge=wedge, WCS=WCS,\
		  labelX=labelX, labelY=labelY, caption=caption, \
		  limitsX = limitsX, limitsY = limitsY, limitsZ = limitsZ,\
		  nan=1, style=style,aspect=aspect, overplot=overplot,\
		  doContour=doContour,levels = levels, labelContour=labelContour,
                  noerase=noerase)

        if showRms:
            # compute rms/beam
            Xrange, Yrange = [],[]
            if limitsX:
                x1,y1 = self.wcs2pix(limitsX[0],WCS['CRVAL2'])
                x2,y2 = self.wcs2pix(limitsX[1],WCS['CRVAL2'])
                Xrange = [int(x1)+1,int(x2)]
            if limitsY:
                x1,y1 = self.wcs2pix(WCS['CRVAL1'],limitsY[0])
                x2,y2 = self.wcs2pix(WCS['CRVAL1'],limitsY[1])
                Yrange = [int(y1),int(y2)+1]

            self.computeRms(rmsKappa=rmsKappa,
                            limitsX=Xrange,limitsY=Yrange)
        
    def smoothBy(self,Size,norm='peak'):
       """
       DES: Smooth the image with a 2D Gaussian of given FWHM.
            Smoothing is peak-normalised, therefore conserves Jy/beam
            as unit.
       INP: (float) Size : the FWHM of the smoothing gaussian
            (str)   norm : normalize to peak ('peak') or integrated
                           ('int') flux. Default is to nomalize to peak
                           flux.
       """
       
       # Compute a kernel of the same pixelSize as the image
       pixsize = abs(self.WCS['CDELT2'])
       smoothingKernel = Kernel(pixsize,Size)
       self.smoothWith(smoothingKernel)
       # Normalisation of the data plane: pixel unit is flux/beam
       # so, update to the new beam size
       if (norm=='peak'):
           newbeam = sqrt(self.BeamSize**2 + Size**2)
       else:
           newbeam = self.BeamSize
           
       self.Data   *= array((newbeam**2)/ (self.BeamSize**2),'f')
       self.Weight *= array((newbeam**2)/ (self.BeamSize**2),'f')
       self.BeamSize = newbeam

    def smoothWith(self,kernel):
       """
       DES: smooth the image with the given kernel
       INP: (kernel) : the kernel
       """
       
       if abs(self.WCS['CDELT1']) == abs(kernel.WCS['CDELT1']) and \
              abs(self.WCS['CDELT2']) == abs(kernel.WCS['CDELT2']):
           self.Data     = fMap.ksmooth(self.Data,     kernel.Data)
           self.Weight   = fMap.ksmooth(self.Weight,   kernel.Data)
           self.Coverage = fMap.ksmooth(self.Coverage, kernel.Data)


    def blankRegion(self,ccord,radius,outside=0):
        """
        DES: selects a circular region on the map and blanks
             the region, or everything outside the region
        INP: (f list) ccord  : x,y world coordinates of center
             (f)      radius : radius of the region to blank
             (bool)  outside : blank outside region
                   
        """
        dat=array(self.Data)

        # convert radius to number of pixels in x and y directions
        npixx=abs(radius/self.WCS['CDELT1'])
        npixy=abs(radius/self.WCS['CDELT2'])
        
        # convert central coordinates to pixel coordinates
        (cpixx,cpixy)=self.wcs2pix(ccord[0],ccord[1])

        xcoord=copy.deepcopy(dat)
        ycoord=copy.deepcopy(dat)

        for i in range(self.WCS['NAXIS1']):
            for j in range(self.WCS['NAXIS2']):
                xcoord[i,j]=i+1
                ycoord[i,j]=j+1

        pixdist=sqrt(((xcoord-cpixx)**2)/(npixx**2)+\
                     ((ycoord-cpixy)**2)/(npixy**2))
        mask=where(pixdist < 1.,1,0)

        print cpixx,cpixy,npixx,npixy

        if outside:
            mask=where(mask == 0,1,0)

        self.blankOnMask(mask)
        return shape(compress(ravel(mask),ravel(mask)))[0]
        

    def blankOnMask(self,mask):
        """
        DES: cut the map according to an input mask
        INP: (array) mask : input mask 
        """

        dat=copy.deepcopy(array(self.Data))
        weight=copy.deepcopy(array(self.Weight))
        cover=copy.deepcopy(array(self.Coverage))

        putmask(dat,mask,float('NaN'))
        putmask(weight,mask,0.)
        putmask(cover,mask,0)
        
        self.Data=dat#.tolist()
        self.Weight=weight#.tolist()
        self.Coverage=cover#.tolist()

        nset=shape(compress(ravel(mask),ravel(mask)))[0]
        return nset

    def blank(self,below=float('NaN'),above=float('NaN')):
        """
        DES: cut the map below and/or above a threshold
        INP: (f) below : cut below this value
             (f) above : cut above this value
        """

        mapdata=copy.deepcopy(array(self.Data))

        mask=where(bitwise_or((mapdata < below),(mapdata > above)),1,0)

        nset=self.blankOnMask(mask)
        return nset

    def blankSigma(self,below=float('NaN'),above=float('NaN'),snmap=1,cell=15,sparse=8):
        """
        DES: cut the map below and/or above a number of sigmas of the s/n map (default) or the map
        INP: (f) below : cut below this value
             (f) above : cut above this value
             (str) mode : 'sn' to use s/n map (with local rms); 'map' to use actual data (with overall rms)
        """

        if snmap:
            mapdata=self.computeSNMap(cell=cell,sparse=sparse)
        else:
            mapdata=copy.deepcopy(array(self.Data))

        good_data=compressNan(mapdata)[0]

        if snmap:
            mask=where(bitwise_or((mapdata < below),(mapdata > above)),1,0)
        else:
            mean=fStat.f_mean(good_data)
            rms=fStat.f_rms(good_data,mean)
            mask=where(bitwise_or((mapdata < below*rms),(mapdata > above*rms)),1,0)

        nset=self.blankOnMask(mask)
        return nset

    def sigmaClip(self,above=5,below=-5):
        """
        DES: despike (sigma clip) a map
        INP: (f) below : cut below the rms times this value
             (f) above : cut above the rms times this value
        """

        good_data=compressNan(self.Data)[0]
        mean=fStat.f_mean(good_data)
        rms=fStat.f_rms(good_data,mean)

        nset = self.blank(below=below*rms,above=above*rms)
        self._Image__MessHand.info('flagging %s pixels in map' % nset)

        return nset

    def iterativeSigmaClip(self,above=5,below=-5,maxIter=10):
        """
        DES: despike (sigma clip) a map iteratively
        INP: (f) below : cut below the rms times this value
             (f) above : cut above the rms times this value
             (i) maxIter : maximum number of iterations
        """

        despiked=1
        i=0
    
        # Run loop
        while ((despiked > 0) and (i < maxIter)):
            i += 1
            despiked = self.sigmaClip(above=above,below=below)
    #----------------------------------------------------------------------------
    def setValuesOnMask(self,mask,value):
        """
        DES: reassign values to the map according to an input mask
        INP: (array) mask : input mask
             (float) value: the value to be used for reassignment  
        """

        dat=copy.deepcopy(array(self.Data))
        weight=copy.deepcopy(array(self.Weight))
        cover=copy.deepcopy(array(self.Coverage))

        putmask(dat,mask,value)
        putmask(weight,mask,0.)
        putmask(cover,mask,0)
        
        self.Data=dat#.tolist()
        self.Weight=weight#.tolist()
        self.Coverage=cover#.tolist()

        nset=shape(compress(ravel(mask),ravel(mask)))[0]
        return nset
    #----------------------------------------------------------------------------    
    def setValues(self,below=float('NaN'),above=float('NaN'),value=float('NaN')):
        """
        DES: cut the map below and/or above a threshold
        INP: (f) below : cut below this value
             (f) above : cut above this value
        """

        mapdata=copy.deepcopy(array(self.Data))

        mask=where(bitwise_or((mapdata < below),(mapdata > above)),1,0)

        nset=self.setValuesOnMask(mask,value)

        return nset
    #----------------------------------------------------------------------------    
   
    def getPixel(self,nbPix=3):
        """
        DES: allow user to get pixel values using mouse
        INP: (int)  nbPix : size of area to compute average (default 3x3)
        """

        self.__MessHand.info("Click left to get one pixel, mid to get average over " +\
                           str(nbPix*nbPix)+", right to exit (on Data array only)")

        data = self.Data
        WCS  = self.WCS
        
        x,y = 0,0
        char = ''
        while(char != 'X'):
            if (char == 'D'):
                #average over [x1,x2] x [y1,y2]
                offset = (nbPix-1)/2.
                if offset <= i < WCS['NAXIS1']-offset and \
                       offset <= j < WCS['NAXIS2']-offset:
                    x1 = int(i-offset)
                    x2 = int(i+offset+1)
                    y1 = int(j-offset)
                    y2 = int(j+offset+1)
                    
                    listVal = self.Data[x1:x2,y1:y2]
                    listVal = ravel(listVal)
                    
                    listVal, nNan = compressNan([listVal])
                    if len(listVal) > 0:
                        val = sum(listVal)/float(len(listVal))
                        self.__MessHand.info("[%ix%i]-%i: (%6.2f,%6.2f) = %6.2f" \
                                             %(nbPix,nbPix,len(listVal),x,y,val))
                else:
                    self.__MessHand.warning("some pixel(s) are outside the plot")
            elif (char == 'A'):
                # Print a single pixel 
                if 0 <= i < WCS['NAXIS1'] and 0 <= j < WCS['NAXIS2']:
                    val = self.Data[i,j]
                    self.__MessHand.info("("+str("%7.3f" % x)+","+str("%7.3f" % y)+") = "+str(val))
                else:
                    self.__MessHand.warning("pixel outside the plot")

            x,y,char = Plot.getpix(x,y)
            i,j = self.wcs2pix(x,y)
            i = int(i)
            j = int(j)

    #--------------------------------------------------------------------------------
    def zoom(self,mouse=1,style='idl4',wedge=1,\
                limitsZ=[],aspect=0,limitsX=[],limitsY=[],caption=None,\
                doContour=0,levels=[],showRms=1,rmsKappa=3.5):
        """
        DES: allow the user to select a region in the map to zoom in
        INP: (bool) mouse: use the mouse? (default: yes)
               (other parameters: same as display)
        """

        if not self.Data:
            self.__MessHand.error('No map computed yet')
            return

        if not mouse:
            if not limitsX or not limitsY:
                self.MessHand.error('X/Y limits must be given if not using the mouse')
                return
            else:
                limX = limitsX
                limY = limitsY
        else:
            self.__MessHand.info('Use the mouse cursor to select 2 opposite corners')
            WCS  = self.WCS        
            x,y = 0,0
            char = ''
            while char == '':
                x,y,char = Plot.getpix(x,y)
                i,j = self.wcs2pix(x,y)
                if not(0 <= i < WCS['NAXIS1'] and 0 <= j < WCS['NAXIS2']):
                    self.__MessHand.warning("pixel outside the plot")
                    char = ''
            x1,y1 = x,y
            self.__MessHand.longinfo("Corner1: ("+str("%8.4f" % x)+","+str("%8.4f" % y)+")")
            char = ''
            while char == '':
                x,y,char = Plot.getpix(x,y)
                i,j = self.wcs2pix(x,y)
                if not(0 <= i < WCS['NAXIS1'] and 0 <= j < WCS['NAXIS2']):
                    self.__MessHand.warning("pixel outside the plot")
                    char = ''
            x2,y2 = x,y
            self.__MessHand.longinfo("Corner2: ("+str("%8.4f" % x)+","+str("%8.4f" % y)+")")

            # sort limits
            if self.WCS['CTYPE1'][:4] == "OLON":
                # for HO map, X axis goes increasing
                limX = [min([x1,x2]),max([x1,x2])]
            else:
                limX = [max([x1,x2]),min([x1,x2])]                
            limY = [min([y1,y2]),max([y1,y2])]
        if not caption:
            caption = 'Zoom'
        self.__MessHand.info('limitsX=['+str("%8.4f" % max(limX))+','+ \
                             str("%8.4f" % min(limX)+'],')+ \
                             'limitsY=['+str("%8.4f" % min(limY))+','+ \
                             str("%8.4f" % max(limY)+'] '))
        self.display(style=style,
                     limitsX=limX,limitsY=limY,limitsZ=limitsZ,
                     aspect=aspect,caption=caption,wedge=wedge,
                     doContour=doContour,levels=levels,
                     showRms=showRms, rmsKappa=rmsKappa)
            

    #--------------------------------------------------------------------------------
    def extractSource(self,gradient=0, circular=0, radius=-10, \
                      Xpos = 0., Ypos = 0., fixedPos = 0, incl=0., fixIncl=0):
        """
        DES: fit a 2D Gaussian on  a map -
        """
        if not self.Data:
            self.__MessHand.error('No map computed yet')
            return
        
        fwhm = self.BeamSize  # that's in the map coordinates unit
        # If map in EQ system, use central position as starting point
        if self.WCS['CTYPE1'].find('RA') >= 0 or self.WCS['CTYPE2'].find('DEC') >= 0:
            if not Xpos:  # but only if not provided by the user
                Xref = self.WCS['CRVAL1']
            else:
                Xref = Xpos
            if not Ypos:
                Yref = self.WCS['CRVAL2']
            else:
                Yref = Ypos
        else:
            Xref, Yref = Xpos, Ypos
        
        mapArray          = self.Data
        weightArray       = self.Weight
        azArray, elArray  = self.physicalCoordinates()
        
        # put everything into 1D easier for radius compression
        mapArray    = ravel(mapArray)
        weightArray = ravel(weightArray)
        azArray     = ravel(azArray)
        elArray     = ravel(elArray)

        if radius < 0:
            radius = abs(radius)*fwhm

        # Select data to a certain radius from the given X/Y pos
        if radius != 0:
            good_for_fit = less(sqrt((azArray-Xref)**2+(elArray-Yref)**2),radius)
            azArray,n   = fUtilities.compress(azArray,good_for_fit,1)
            azArray     = azArray[:n]
            elArray,n   = fUtilities.compress(elArray,good_for_fit,1)
            elArray     = elArray[:n]
            mapArray,n = fUtilities.compress(mapArray,good_for_fit,1)
            mapArray   = mapArray[:n]
            weightArray,n = fUtilities.compress(weightArray,good_for_fit,1)
            weightArray   = weightArray[:n]

        # now use the data to fit a 2D-Gaussian
        try:
            pointingResult = fitBaseEllipticalGaussian(mapArray,azArray,elArray,\
                                                       err = 1/sqrt(weightArray),\
                                                       fwhm=fwhm,\
                                                       gradient=gradient,\
                                                       circular=circular,\
                                                       Xpos=Xpos,Ypos=Ypos,\
                                                       fixedPos=fixedPos,\
                                                       incl=incl,\
                                                       fixIncl=fixIncl)
        except BoaError, error:
            self.__MessHand.warning('fit did not converge:')
            self.__MessHand.warning(error.msg)
            pointingResult = -1

        # TODO: is it the best thing to return the result?
        return pointingResult
    

    #--------------------------------------------------------------------------
    #
    # methods to compute rms/beam
    #--------------------------------------------------------------------------
    def rmsDistribution(self,cell=3):
        """
        DES: compute and plot the distribution of rms in the map
        INP: (int)  cell : size of cells on which rms are computed (default: 3x3)
        """
        
        data = self.Data
        rmsdistr = copy.deepcopy(self.Data)
        WCS  = self.WCS
        
        allRms = []
        for x in range(WCS['NAXIS1']-cell):
            for y in range(WCS['NAXIS2']-cell):
                listVal = self.Data[x:x+cell,y:y+cell]
                listVal = ravel(listVal)
                    
                listVal, nNan = compressNan([listVal])
                if len(listVal) > 0:
                    mean = fStat.f_mean(listVal)
                    val = fStat.f_rms(listVal,mean)
                    allRms.append(val)
                    
        return allRms

    def rmsMap(self,cell=15,sparse=8):
        """
        DES: compute the distribution of rms in the map 
        INP: (int)   cell : size of cells on which rms are computed (default: 15x15)
             (int) sparse : compute rms only on pixels separated by this number
                            (to save time) (default: 8)
        """
        
        data = self.Data
        rmsdistr = copy.deepcopy(self.Data)*0.0+float('nan')
        WCS  = self.WCS

        off1=int(cell/2.)

        maxx=int(WCS['NAXIS1'])
        maxy=int(WCS['NAXIS2'])

        # find x values to compute rms
        xval=[off1+1]
        done=0
        while (done == 0):
            newx=xval[shape(xval)[0]-1]+sparse
            if (newx < maxx-off1-1):
                xval.extend([newx])
            else:
                done=1
        # same for y
        yval=[off1+1]
        done=0
        while (done == 0):
            newy=yval[shape(yval)[0]-1]+sparse
            if (newy < maxy-off1-1):
                yval.extend([newy])
            else:
                done=1
        
        for x in xval:
            for y in yval:
                
                listVal = self.Data[x-off1:x+off1,y-off1:y+off1]
                listVal = ravel(listVal)
                    
                listVal, nNan = compressNan([listVal])
                if len(listVal) > 0:
                    mean = fStat.f_mean(listVal)
                    val = fStat.f_rms(listVal,mean)
                    rmsdistr[x,y]=val
                
        return rmsdistr

    def computeSNMap(self,cell=15,sparse=8):
        """
        DES: compute a signal-to-noise map from the current map data and weights
        INP: (int)   cell : size of cells on which rms are computed (default: 10x10)
             (int) sparse : compute rms only on pixels separated by this number (to save time) (default: 5)
        """

        # we assume that the local rms in the map is exactly 1./sqrt(weight) times some constant.
        # to determine the constant, we have to first find the local rms in the map
        # if there is a strong source, one should cut it out first.

        dat=copy.deepcopy(array(self.Data))
        weight=copy.deepcopy(array(self.Weight))

        # compute rms map
        rmsmap=ravel(self.rmsMap(cell=cell,sparse=sparse))

        # compute unscaled rms map from weights
        rmsunsc=1./sqrt(ravel(weight))

        # find where both of these maps are ok
        mask1=where(bitwise_and((rmsmap > -10000000.),(rmsunsc > -10000000.)),1,0)

        mask2=where(bitwise_and((rmsmap != 0),(rmsunsc != 0)),1,0)
        mask=mask1*mask2

        rmsmap=compress(mask,rmsmap)
        rmsunsc=compress(mask,rmsunsc)

        # the coefficient should be...
        coeffs=rmsmap/rmsunsc

        #take the median
        coeff=fStat.f_median(coeffs)
        self._Image__MessHand.info('computed rms = %f * (1./sqrt(weight))'%coeff)

        # TBD: remove outliers??

        # return the s/n map
        self.snmap=dat/((1./sqrt((weight)))*coeff)  
        return dat/((1./sqrt((weight)))*coeff)       
        
                    
    def meanDistribution(self,cell=3,limitsX=[],limitsY=[]):
        """
        DES: compute and plot the distribution of means in the map
        INP: (int)  cell : size of cells on which mean values are computed (default: 3x3)
             (i lists) limitsX/Y: optionally define a sub-region (pixel coord)
        """
        
        WCS  = self.WCS
        
        # Number of beams in the map
        nbX = WCS['NAXIS1']-cell
        nbY = WCS['NAXIS2']-cell

        # if sub-region asked, check that coord are within the map
        if limitsX:
            if limitsX[0] < 0:
                limitsX[0] = 0
            if limitsX[1] > nbX:
                limitsX[1] = nbX
        if limitsY:
            if limitsY[0] < 0:
                limitsY[0] = 0
            if limitsY[1] > nbY:
                limitsY[1] = nbY

        theData = copy.copy(self.Data)
        if limitsX:
            theData = theData[limitsX[0]:limitsX[1],::]
            nbX = limitsX[1]-limitsX[0]
        if limitsY:
            theData = theData[:,limitsY[0]:limitsY[1]]
            nbY = limitsY[1]-limitsY[0]

        allMean,nbOk = fStat.meandistribution(theData,nbX,nbY,nbX*nbY,cell)
        allMean = allMean[:nbOk]
        
        return allMean

    def computeRmsBeam(self,cell=3,rmsKappa=3.5,limitsX=[],limitsY=[]):
        """
        DES: compute rms/beam in a map (smoothed at beam resolution)
        INP: (f)     cell: size of one beam in pixel
             (f) rmsKappa: for kappa-sigma clipping before computing rms
             (i lists) limitsX/Y: optionally define a sub-region (pixel coord)
        """
        # Distribution of flux/beam
        m = self.meanDistribution(cell=cell,limitsX=limitsX,limitsY=limitsY)
        # Remove the NaNs at this stage
        mm, nNan = compressNan([m])
        # clip data at more than k-sigma (e.g. a source)
        m2,nbOk = fStat.clipping(mm,rmsKappa)
        m2 = m2[:nbOk]
        self.RmsBeam = fStat.f_stat(m2)[1]
        self.__MessHand.info("r.m.s. / beam = %f"%(self.RmsBeam))

    def computeRms(self,rmsKappa=3.5,limitsX=[],limitsY=[]):
        """
        DES: compute rms/beam in a map (dispersion between pixels)
        INP: (f) rmsKappa: for kappa-sigma clipping before computing rms
             (i lists) limitsX/Y: optionally define a sub-region (pixel coord)
        """
        if limitsX or limitsY:
            if not limitsX:
                localData = self.Data[:,limitsY[0]:limitsY[1]]
            elif not limitsY:
                localData = self.Data[limitsX[0]:limitsX[1],::]
            else:
                localData = self.Data[limitsX[0]:limitsX[1],
                                      limitsY[0]:limitsY[1]]
        else:
            localData = self.Data
        # Distribution of pixel values
        pixels = ravel(localData)
        # Remove the NaNs at this stage
        pixels,nb = compressNan([pixels])
        # clip data at more than k-sigma (e.g. a source)
        m2,nbOk = fStat.clipping(pixels,rmsKappa)
        m2 = m2[:nbOk]
        self.RmsBeam = fStat.f_stat(m2)[3]  # median deviation
        self.__MessHand.info("map r.m.s. = %f"%(self.RmsBeam))

    def submap(self,limitsX=[],limitsY=[]):
        """
        DES: this function returns a map covering a sub-region of the
             initial map
        INP: (f list) limitsX/Y: the limits in world coordinates
        OUT: an object of class Image is returned
        """
        localMap = copy.deepcopy(self)
        if limitsX or limitsY:
            if not limitsX:
                # only limitsY given
                x1,x2 = 0,self.WCS['NAXIS1']
                i,y1 = self.wcs2pix(self.WCS['CRVAL1'],min(limitsY))
                i,y2 = self.wcs2pix(self.WCS['CRVAL1'],max(limitsY))
            elif not limitsY:
                # only limitsX given
                y1,y2 = 0,self.WCS['NAXIS2']
                x1,j = self.wcs2pix(max(limitsX),self.WCS['CRVAL2'])
                x2,j = self.wcs2pix(min(limitsX),self.WCS['CRVAL2'])
            else:
                # both limits given by user
                x1,y1 = self.wcs2pix(max(limitsX),min(limitsY))
                x2,y2 = self.wcs2pix(min(limitsX),max(limitsY))

            # make sure that the limits don't go outside the existing data
            x1 = max([0,int(x1)])
            x2 = min([self.WCS['NAXIS1'],int(x2)+1])
            y1 = max([0,int(y1)])
            y2 = min([self.WCS['NAXIS2'],int(y2)+1])
            localMap.Data     = self.Data[x1:x2,y1:y2]
            localMap.Weight   = self.Weight[x1:x2,y1:y2]
            localMap.Coverage = self.Coverage[x1:x2,y1:y2]
            # update WCS info - need to compute the actual limits
            limX1,limY1 = self.wcs2phy(x1,y1)
            limX2,limY2 = self.wcs2phy(x2,y2)
            minmaxXY = [limX1,limX2,limY1,limY2]
            pixSize = abs(self.WCS['CDELT1'])
            localMap.computeWCS(pixSize,minmax=minmaxXY)

        return localMap
    
    def writeFits(self,outfile='boaMap.fits',
                  overwrite=0,limitsX=[],limitsY=[],intensityUnit="Jy/beam",
                  writeFlux=1,writeWeight=1,writeCoverage=1,
                  writePSF=0,writePST=0,
                  writeRms=0,rmsfile=''):
        """
        DES: store the current map (2D array with WCS info) to a FITS file
        INP: (str)   outfile: output file name (default boaMap.fits)
             (bool) overwrite: overwrite existing file -
                              default = 0: do not overwrite existing file
             (f list) limitsX/Y: optional map limits (in world coordinates)
             (string) intensityUnit: optional unit of the intensity (default: "Jy/beam")
             (bool) writeFlux, writeWeight, writeCoverage: should these planes be
                                            included in the output file? (def. yes)
             (bool) writePSF, writePST: should these planes be
                                            included in the output file? (def. no)
             (bool) writeRms: should another file with rms plane be written? (def. no)
        """

        if os.path.exists(outfile):
            if not overwrite:
                self.__MessHand.error('File %s exists' % outfile)
                return
        
        try:
            dataset = BoaFits.createDataset("!" + outfile)
        except Exception, data:
            self.__MessHand.error('Could not open file %s: %s' % (outfile, data))
            return

        if limitsX or limitsY:
            localMap = self.submap(limitsX=limitsX,limitsY=limitsY)
        else:
            localMap = self
        
        try:
            if writeFlux:
                localMap.__writeImage(dataset, "Intensity", intensityUnit=intensityUnit)
            if writeWeight:
                localMap.__writeImage(dataset, "Weight")
            if writeCoverage:
                localMap.__writeImage(dataset, "Coverage")
            if writePSF:
                localMap.__writeImage(dataset, "PSF")
            if writePST:
                localMap.__writeImage(dataset, "PST")
            dataset.close()
            
        except Exception, data:
            self.__MessHand.error('Could not write data to file %s: %s' % (outfile, data))
            return

        if writeRms:
            localMap = copy.deepcopy(localMap)  # do not overwrite input map
            localMap.Data = array(1.,'f')/sqrt(localMap.Weight)
            if not rmsfile:
                rmsfile = outfile[:-5]+'-rms.fits'
            dataset = BoaFits.createDataset("!" + rmsfile)
            localMap.__writeImage(dataset, "RMS", intensityUnit=intensityUnit)
            dataset.close()
                        
        localMap = 0  # free memory
        dataset = 0
        
    def __writeImage(self, dataset, extname, intensityUnit=""):
        if extname == "Intensity":
            data = self.Data
        elif extname == "Weight":
            data = self.Weight
        elif extname == "Coverage":
            data = self.Coverage
        elif extname == "RMS":
            data = self.Data
        elif extname == 'PST':
            data = self.PST
        elif extname == 'PSF':
            data = self.PSF
        else:
            self.__MessHand.error("Don't know how to write image '%s'" % extname)

        if not data:
            self.__MessHand.info("No data to write image %s" % extname)
            return

        WCS = self.WCS
        
        try:
            dataImage = dataset.createImage(data)
            if not dataImage:
                self.__MessHand.error("Error while creating image '%s'" % extname)
                return

            dataImage.createKeywordDate()
            dataImage.createKeyword("CREATOR","BoA",
                                    comment="Bolometer Data Analysis Project")
            dataImage.createKeyword("EXTNAME",extname,
                                    comment="Type of data contained in this image")

            dataImage.createKeyword("MJD-OBS",WCS["MJD-OBS"],
                                    comment="time of observation")
            dataImage.createKeyword("RADESYS",WCS["RADESYS"],
                                    comment="frame of reference")
            dataImage.createKeyword("EQUINOX",WCS["EQUINOX"],
                                    comment="coordinate epoch")

            dataImage.createKeyword("CTYPE1",WCS["CTYPE1"],
                                    comment="Type of coordinate 1")
            dataImage.createKeyword("CRPIX1",WCS["CRPIX1"]+1,
                                    comment="Reference pixel of coordinate 1")
            dataImage.createKeyword("CDELT1",WCS["CDELT1"],
                                    comment="Increment per pixel of coordinate 1")
            dataImage.createKeyword("CRVAL1",WCS["CRVAL1"],
                                    comment="Value of coordinate 1 at reference point")
            dataImage.createKeyword("CUNIT1",WCS["CUNIT1"],
                                    comment="Unit of coordinate 1")
            
            dataImage.createKeyword("CTYPE2",WCS["CTYPE2"],
                                    comment="Type of coordinate 2")
            dataImage.createKeyword("CRPIX2",WCS["CRPIX2"]+1,
                                    comment="Reference pixel of coordinate 2")
            dataImage.createKeyword("CDELT2",WCS["CDELT2"],
                                    comment="Increment per pixel of coordinate 2")
            dataImage.createKeyword("CRVAL2",WCS["CRVAL2"],
                                    comment="Value of coordinate 2 at reference point")
            dataImage.createKeyword("CUNIT2",WCS["CUNIT2"],
                                    comment="Unit of coordinate 2")

            dataImage.createKeyword("PC1_1",WCS["PC1_1"],
                                    comment="Linear transformation matrix")
            dataImage.createKeyword("PC1_2",WCS["PC1_2"],
                                    comment="Linear transformation matrix")
            dataImage.createKeyword("PC2_1",WCS["PC2_1"],
                                    comment="Linear transformation matrix")
            dataImage.createKeyword("PC2_2",WCS["PC2_2"],
                                    comment="Linear transformation matrix")
            
            dataImage.createKeyword("OBJECT",self.Header["Object"],
                                    comment="Object observed")
            dataImage.createKeyword("TELESCOP",self.Header["Telescope"],
                                    comment="Telescope name")
            dataImage.createKeyword("FEBE",self.Header["FeBe"],
                                    comment="Frontend-backend combination")

            # BUNIT: this assumes that calibration was done!
            if extname == "Intensity":
                dataImage.createKeyword("BUNIT",intensityUnit,
                                        comment="Physical unit of image")
            dataImage.createKeyword("BMAJ",self.BeamSize,
                                    comment="Beam major axis")
            dataImage.createKeyword("BMIN",self.BeamSize,
                                    comment="Beam minor axis")
            dataImage.createKeyword("BPA",0.,
                                    comment="Beam position angle")

            dataImage.writeImage()
            
            self.__MessHand.longinfo("Image '%s' written. Size: %s"
                                   % (extname, str(dataImage.getShape())))

            data = 0  # free memory
            dataImage = 0
            
        except Exception, msg:
            self.__MessHand.error("Error while writing image '%s'" % extname)
            self.__MessHand.error("Exception: '%s'" % msg)
            return
        
    # -------------------------------------------------------------------
    def dumpMap(self,fileName='BoaMap.sav'):
        """
        DES: save an Image instance to a file
        INP: (string) fileName: name of the output file
                      (default = 'BoaMap.sav')
        """
        try:
            f = file(BoaConfig.outDir+fileName,'w')
        except:
            self.__MessHand.error(" permission denied, please change outdir")
            return
        cPickle.dump(self,f,2)
        f.close()
        self.__MessHand.info(" current Image successfully written to %s"%fileName)


#----------------------------------------------------------------------------------
#----- Kernel Class ---------------------------------------------------------------
#----------------------------------------------------------------------------------
class Kernel(Image):
   """
   NAM: Kernel (class)
   DES: define a kernel
   """

   def __init__(self, pixelSize, beamSize):
	"""
	DES: Initialise an instance of a Kernel class
	INP: (float) pixelSize: the physical size of a pixel
             (float) beamSize : the beam FWHM in the same unit
	"""
	
	Image.__init__(self)
	
	self.computeWCS(pixelSize,\
			minmax = array([-1,1,-1,1],Float32)*beamSize*3.)

        # For smoothing, we NEED an odd number of pixels
        # if that's not the case, then generate a slightly larger kernel
        nbX = self.WCS['NAXIS1']
        if not nbX%2:
            self.computeWCS(pixelSize,\
                            minmax = array([-1,1,-1,1],Float32)*beamSize*3. +
                            array([-1,1,-1,1],Float32)*pixelSize/2.)           
        
	kernel_azimuth, kernel_elevation = self.physicalCoordinates()

	#kernel_parameters = array([0.,0.,0.,beamSize**2*16*pi*log(2.),
	kernel_parameters = array([0.,0.,0.,1.,
				   0.,0.,beamSize,beamSize,0.])
	tmpData = modelBaseEllipticalGaussian(kernel_parameters,\
                                              [kernel_azimuth,kernel_elevation])
        self.Data = tmpData.astype(Float32)

	
#----------------------------------------------------------------------------------
#----- Boa Mapping Class ----------------------------------------------------------
#----------------------------------------------------------------------------------
class Map(BoaDataAnalyser.DataAna):
    """
    NAM: Map (class)
    DES: An object of this class is responsible for the restoration of
         mapping data of single or multiple files.
    """
    
    def __init__(self):
        """
        DES: Initialise an instance.
        """
        BoaDataAnalyser.DataAna.__init__(self)
        
        self.__Results = []
	self.Map = Image()
        

#--------------------------------------------------------------------------------
#----- map methods --------------------------------------------------------------
#--------------------------------------------------------------------------------

    def showMap(self,style='idl4',wedge=1,\
                limitsZ=[],aspect=0,limitsX=[],limitsY=[],caption=None,\
                doContour=0,levels=[],showRms=1,rmsKappa=3.5,noerase=0):
        """
        DES: show the reconstructed map in (Az,El) or (Ra,Dec)
        """

        if self.Map.Data == []:
            self.MessHand.error('No map computed yet')
            return
	if not caption:
            caption = self.ScanParam.caption()
            
	self.Map.display(style=style,
			 limitsX=limitsX,limitsY=limitsY,limitsZ=limitsZ,
			 aspect=aspect,caption=caption,wedge=wedge,
			 doContour=doContour,levels=levels,
                         showRms=showRms, rmsKappa=rmsKappa,
                         noerase=noerase)

    #--------------------------------------------------------------------------------

    def displayMap(self, weight=0, coverage=0, style='idl4', caption='', wedge=1, \
		aspect=0, overplot=0, doContour=0, levels=[], labelContour=0,\
		limitsX = [], limitsY = [], limitsZ = [], \
                showRms=0, rmsKappa=3.5):

        if self.Map.Data == []:
            self.MessHand.error('No map computed yet')
            return

        self.Map.display(weight=weight, coverage=coverage, style=style, \
                         caption=caption, wedge=wedge, aspect=aspect, \
                         overplot=overplot, doContour=doContour, levels=levels, \
                         labelContour=labelContour, \
                         limitsX = limitsX, limitsY = limitsY, limitsZ = limitsZ, \
                         showRms=showRms, rmsKappa=rmsKappa )

    #--------------------------------------------------------------------------------

    def smoothMap(self,Size):
       """
       DES: smooth the image with a 2D gaussian of given FWHM
       INP: (float) Size : the FWHM of the smoothing gaussian
       """
       
       if self.Map.Data == []:
           self.MessHand.error('No map computed yet')
           return

       self.Map.smoothBy(Size)

    #--------------------------------------------------------------------------------

    def getPixelFromMap(self,nbPix=3):
        """
        DES: allow user to get pixel values using mouse
        INP: (int)  nbPix : size of area to compute average (default 3x3)
        """
       
        if self.Map.Data == []:
            self.MessHand.error('No map computed yet')
            return

        self.Map.getPixel(nbPix)

    #--------------------------------------------------------------------------------

    def computeRmsFromMap(self,rmsKappa=3.5,limitsX=[],limitsY=[]):
        """
        DES: compute rms/beam in a map (dispersion between pixels)
        INP: (f) rmsKappa: for kappa-sigma clipping before computing rms
             (i lists) limitsX/Y: optionally define a sub-region (pixel coord)
        """
       
        if self.Map.Data == []:
            self.MessHand.error('No map computed yet')
            return

        self.Map.computeRms(rmsKappa,limitsX,limitsY)

#--------------------------------------------------------------------------------
#----- public methods -----------------------------------------------------------
#--------------------------------------------------------------------------------

    def doMap(self,chanList=[], channelFlag=[], plotFlaggedChannels=0,
              dataFlag=[], plotFlaggedData=0, oversamp=2.0, beammap = 0,
              system='HO', sizeX=[], sizeY=[], limitsZ=[], style='idl4',
              wedge=1, smooth=0, noPlot=0, caption=None, aspect=0,
              showRms=1, rmsKappa=3.5, derotate=0, neighbour=0, relative=1):

        """
        DES: reconstruct a map in (Az,El) coordinates combining bolometers
        INP: (int list) chanList = channels to consider
             (integer list) channelFlag : plot data from channels flagged or unflagged accordingly
             (log)  plotFlaggedChannels : channelFlag revers to flagged/unflagged data
             (integer list)    dataFlag : plot data flagged or unflagged accordingly
             (log)      plotFlaggedData : dataFlag revers to flagged/unflagged data
                                          flag   | plotFlagged | Plot..
                                          'None' |  0          | all data
                                          []     |  0          | unflagged data (default)
                                          []     |  1          | data with at least one flag set
                                          1      |  0          | data with flag 1 not set
                                          1      |  1          | data with flag 1 set
                                          [1,2]  |  0          | data with neither flag 1 nor flag 2 set
                                          [1,2]  |  1          | data with either flag 1 or flag 2 set
             (float)    oversamp = oversampling factor (beam fwhm / pixel size). Default=2.
             (log)       beammap = compute a beam map (default: no)
             (str)        system = coord. system, one of 'HO' (Az,El *offsets*) or 'EQ'
                                   (RA, Dec absolute coord.) or 'GAL' (Galactic)
                                   default = 'HO'
                                   optionally 'EQFAST' to do only one rotation
                                   on small maps (faster)
             (list float)  sizeX = limits in Az of the map
             (list float)  sizeY = limits in El of the map
             (logical)     noNan = remove NaN in self.Results?
             (str)         style = color table to use in image
             (logical)    smooth = do we smooth with beam? (default: no)
	     (logical)    noPlot = do not plot the map? (default: no, i.e. yes we plot)
	     (str)       caption = plot caption
             (logical)    aspect = keep aspect ratio? (default: yes)
             (logical)   showRms = compute and print rms/beam? (default: yes)
             (float)    rmsKappa = kappa in kappa-sigma clipping used to compute rms
             (int)      derotate = derotate Nasmyth array by Elevation
             (logical) neighbour = do we divide signal into 4 neighbouring pixels? (def: no)
             (logical)  relative = use bolometer offsets w.r.t. to reference channel
                                   (relative=1, default) or use absolute offsets (relative=0)
        """
        # check channel list; will return all non-flagged channels if chanList = []
        chanList = self.BolometerArray.checkChanList(chanList, \
                                                     flag=channelFlag,getFlagged=plotFlaggedChannels)

        if not(len(chanList)):
            self.MessHand.error("no valid channels")
            return
        chanListIndexes = self.BolometerArray.getChanIndex(chanList)

        myTiming = Timing()
        myTiming.setTime()

        # Select offsets according to system, 'HO' or 'EQ'
        system = string.upper(system)
        if system == 'HO':
            if derotate:
                ElAverage = average(self.ScanParam.El)
                self.BolometerArray.rotateArray(derotate*ElAverage)
            XYOffsets = array([self.ScanParam.get('AzimuthOffset',flag='None'), \
                               self.ScanParam.get('ElevationOffset',flag='None')])
            XOffsets = self.ScanParam.get('AzimuthOffset')
            YOffsets = self.ScanParam.get('ElevationOffset')
        elif (system == 'EQ') or (system=='EQFAST'):
            XYOffsets = array([self.ScanParam.get('RA',flag='None'), \
                               self.ScanParam.get('Dec',flag='None')])
            XOffsets = self.ScanParam.get('RA')
            YOffsets = self.ScanParam.get('Dec')
        elif (system == 'GAL'):
            if not len(self.ScanParam.GLon):
                self.ScanParam.computeGal()
            XYOffsets = array([self.ScanParam.get('Glon',flag='None'), \
                               self.ScanParam.get('Glat',flag='None')])
            XOffsets = self.ScanParam.get('Glon')
            YOffsets = self.ScanParam.get('Glat')
            
 	data        = self.Data
	dataWeights = self.DataWeights

        if relative:
            OffsetsUsed = array(self.BolometerArray.getChanSep(self.BolometerArray.UsedChannels))
        else:
            OffsetsUsed = self.BolometerArray.Offsets
        
        if beammap:
            OffsetsUsed = OffsetsUsed * 0.0
            if system == 'EQ':
                self.MessHand.warning('Doing a beam map in EQ system probably does not make sense')

        # Retrieve the beam and compute pixel size
        fwhm = self.BolometerArray.BeamSize
        pixelSize = fwhm / oversamp

        if (system == 'EQ') or (system == 'EQFAST') or (system == 'GAL'):
            pixelSize = pixelSize/3600.
            fwhm = fwhm/3600.
            
            if (system == 'EQ'):
                if (derotate):
                    rotAngles = array(self.ScanParam.ParAngle)+derotate*array(self.ScanParam.El)
                else:
                    rotAngles = array(self.ScanParam.ParAngle)
            elif (system == 'GAL'):
                #if (derotate):
                #    rotAngles = array(self.ScanParam.ParAngle)+derotate*array(self.ScanParam.El)
                #else:
                if not len(self.ScanParam.GalAngle):
                    self.ScanParam.computeGalAngle()
                rotAngles = array(self.ScanParam.ParAngle) + array(self.ScanParam.GalAngle)

            chanListAzEl = array(self.BolometerArray.UsedChannels)-1
            OffsetsAzEl=array((take(self.BolometerArray.Offsets[0,::],chanListAzEl), \
                                   take(self.BolometerArray.Offsets[1,::],chanListAzEl)))
            refChOffsets=array((self.BolometerArray.RefOffX, self.BolometerArray.RefOffY), 'f')

            # even if system is 'EQ', the following needs to be done in order
            # to allow computation of the extrema of the map
            OffsetsUsed = OffsetsUsed / 3600.
            Yrange = fStat.minmax(YOffsets)
            OffsetsUsed[0,::] /= array(cos((Yrange[0]+Yrange[1])/2. * pi / 180.),'f')
            
        # compute roughly the extrema of the map (taking into account the ScanParam
        # flags but not the BolometerArray Flags nor the Data Flags
	minmaxXY = concatenate(\
		(fStat.minmax(XOffsets)+array([-1,1])*2*pixelSize + \
		 fStat.minmax(take(OffsetsUsed[0,::],chanListIndexes)), \
		 fStat.minmax(YOffsets)+array([-1,1])*2*pixelSize + \
		 fStat.minmax(take(OffsetsUsed[1,::],chanListIndexes)) ) \
		)
        if (system == 'EQ') or (system == 'EQFAST') or (system == 'GAL'):
            # swap X min and max (East to the left)
            minmaxXY[0],minmaxXY[1] = minmaxXY[1],minmaxXY[0]
            
	# Create the resulting image
        Map = self.Map
        Map.__init__()

        Map.Header['Object']    = self.ScanParam.Object
        Map.Header['Telescope'] = self.BolometerArray.Telescope.Name
        Map.Header['FeBe']      = self.BolometerArray.FeBe
        Map.BeamSize = fwhm
        
        if system == 'HO':
            Map.WCS['CTYPE1']   = 'OLON--SFL'
            Map.WCS['CTYPE2']   = 'OLAT--SFL'
            Map.WCS['CUNIT1']   = 'arcsec'
            Map.WCS['CUNIT2']   = 'arcsec'
            cdeltUnit           = 1./3600.
        elif (system == 'EQ') or (system == 'EQFAST'):
            Map.WCS['CTYPE1']   = 'RA---GLS'
            Map.WCS['CTYPE2']   = 'DEC--GLS'
            Map.WCS['CUNIT1']   = 'deg'
            Map.WCS['CUNIT2']   = 'deg'
            cdeltUnit           = 1.
        elif (system == 'GAL'):
            Map.WCS['CTYPE1']   = 'GLON-GLS'
            Map.WCS['CTYPE2']   = 'GLAT-GLS'
            Map.WCS['CUNIT1']   = 'deg'
            Map.WCS['CUNIT2']   = 'deg'
            cdeltUnit           = 1.
        Map.WCS['MJD-OBS']  = self.ScanParam.DateObs
        
        # update main WCS keywords
        Map.computeWCS(pixelSize,sizeX,sizeY,minmaxXY)
	AXIS1 = array([Map.WCS['NAXIS1'],Map.WCS['CRPIX1'],Map.WCS['CDELT1'],
                       Map.WCS['CRVAL1'],cdeltUnit])
	AXIS2 = array([Map.WCS['NAXIS2'],Map.WCS['CRPIX2'],Map.WCS['CDELT2'],
                       Map.WCS['CRVAL2'],cdeltUnit])

        self.MessHand.longinfo("Building a map with dimensions (x,y) = "+\
                           str(int(AXIS1[0]))+','+str(int(AXIS2[0])))

        if Map.WCS['NAXIS1']*Map.WCS['NAXIS2']*32*3/8/1e6 > 600:
            if BoaConfig.online:
                self.MessHand.error('Such maps would require more than 600Mb of memory, aborting...')
                return
            elif not self.MessHand.yesno('do you really want to do such a big map (>600Mb) ?'):
                return
            
        mapData     = zeros((int(AXIS1[0]),int(AXIS2[0])),typecode=Float32)
        mapWeight   = zeros((int(AXIS1[0]),int(AXIS2[0])),typecode=Float32)
        mapCoverage = zeros((int(AXIS1[0]),int(AXIS2[0])),typecode=Float32)
        
        dataFlag = self._removeReservedFlagValues(dataFlag, self.rflags['CHANNEL FLAGGED'])
        if dataFlag==None:
            self.MessHand.error("no valid flags")
            return

        if dataFlag in ['','None']:
            dataFlagMask = ones(shape=self.FlagHandler.getFlags().shape, \
                                typecode=self.FlagHandler.getFlags().typecode())
        else:
            if plotFlaggedData:
                dataFlagMask = self.FlagHandler.isSetMask(dataFlag)
            else:
                dataFlagMask = self.FlagHandler.isUnsetMask(dataFlag)

        if (system=='EQFAST') or ((system=='HO')):
            mapData, mapWeight, mapCoverage = fMap.horizontalprojection(chanListIndexes, data, \
                                                                        dataWeights, \
                                                                        dataFlagMask, 1, \
                                                                        XYOffsets, OffsetsUsed,
                                                                        AXIS1, AXIS2,
                                                                        mapData, mapWeight, mapCoverage)
        if system=='EQ' or system == 'GAL':
            mapData, mapWeight, mapCoverage = fMap.equatorialprojection(chanListIndexes, data, \
                                                                        dataWeights, \
                                                                        dataFlagMask, 1, \
                                                                        XYOffsets, OffsetsAzEl,
                                                                        rotAngles,
                                                                        refChOffsets,
                                                                        AXIS1, AXIS2, neighbour,
                                                                        mapData, mapWeight, mapCoverage)

        Map.Data     = mapData
        Map.Weight   = mapWeight
        Map.Coverage = mapCoverage
        
	self.Map = Map

	# Smoothing after normalisation 
        if smooth:
	    kernel_to_smooth = Kernel(pixelSize,fwhm)
	    self.Map.smoothWith(kernel_to_smooth)

	if not noPlot:
            self.showMap(wedge=wedge,style=style,aspect=aspect,
                         limitsZ=limitsZ,caption=caption,
                         showRms=showRms,rmsKappa=rmsKappa)

        self.MessHand.debug(" map done in "+ str(myTiming))

        if ((system=='HO') and derotate):
            self.BolometerArray.rotateArray(-derotate*ElAverage)   # rotate array back


    #--------------------------------------------------------------------------------
    def flagSource(self,chanList=[],threshold=1.,flag=8,model=None,derotate=0):
        """
        DES: Flag the data according to a model map
        INP: (i list) chanList: the list of channels to work with
             (f)     threshold: the pixel value in input map above which
                                is considered as source
             (i)          flag: the value of flag to set (def: 8)
             (Image object) model: the input model map (with WCS)
                            (default: use current data.Map)
             (i)      derotate: rotate array by El?
        """

        if not model:
            model = self.Map

        if not model.Data:
            self.MessHand.error("no map computed yet, and no model provided")
            return

        chanList = self.BolometerArray.checkChanList(chanList)
        chanListIndexes = self.BolometerArray.getChanIndex(chanList)

        if string.find(model.WCS['CTYPE1'],'GLON') > -1:
            if not len(self.ScanParam.GalAngle):
                if not len(self.ScanParam.GLon):
                    self.ScanParam.computeGal()
                self.ScanParam.computeGalAngle()
            rotAngles = array(self.ScanParam.ParAngle) + array(self.ScanParam.GalAngle)
            XYOffsets = array([self.ScanParam.get('Glon',flag='None'), \
                               self.ScanParam.get('Glat',flag='None')])
        else:
            if (derotate):
                rotAngles = array(self.ScanParam.ParAngle)+derotate*array(self.ScanParam.El)
            else:
                rotAngles = array(self.ScanParam.ParAngle)
            XYOffsets = array([self.ScanParam.get('RA',flag='None'), \
                               self.ScanParam.get('Dec',flag='None')])
        
        chanListAzEl = array(self.BolometerArray.UsedChannels)-1
        OffsetsAzEl=array((take(self.BolometerArray.Offsets[0,::],chanListAzEl), \
                           take(self.BolometerArray.Offsets[1,::],chanListAzEl)))
        refChOffsets=array((self.BolometerArray.RefOffX, self.BolometerArray.RefOffY), 'f')
        AXIS1 = array([model.WCS['NAXIS1'],model.WCS['CRPIX1'],
                       model.WCS['CDELT1'],model.WCS['CRVAL1'],1.])
        AXIS2 = array([model.WCS['NAXIS2'],model.WCS['CRPIX2'],
                       model.WCS['CDELT2'],model.WCS['CRVAL2'],1.])

        # to avoid Segmentation faults with large arrays, one has to allocate the
        # resulting array of flags on the python side
        newflags = zeros(shape(self.Data),'f')
        result = fMap.flagsource(chanListIndexes, self.Data, model.Data,
                                 threshold, flag, XYOffsets, OffsetsAzEl,
                                 rotAngles, refChOffsets, AXIS1, AXIS2, newflags)
        
        mask = fUtilities.as_column_major_storage(result.astype(Int8))
        #self.FlagHandler.setOnMask(mask, flag)
        # here also, to avoid seg. fault, call the fortran routine
        # channel by channel
        for index in chanListIndexes:
            self.FlagHandler.setOnMask(mask[:,index],flag,dim=1,index=index)


    #--------------------------------------------------------------------------------

    def flagSourceOld(self,chanList=[],threshold=1.,flag=8,model=None):
        """
        DES: Flag the data according to a model map
        INP: (i list) chanList: the list of channels to work with
             (f)     threshold: the pixel value in input map above which
                                is considered as source
             (i)          flag: the value of flag to set (def: 8)
             (Image object) model: the input model map (with WCS)
                            (default: use current data.Map)
        """

        if not model:
            model = self.Map

        if not model.Data:
            self.MessHand.error("no map computed yet, and no model provided")
            return

        chanList = self.BolometerArray.checkChanList(chanList)
        chanListIndexes = self.BolometerArray.getChanIndex(chanList)

        XYOffsets = array([self.ScanParam.get('RA',flag='None'), \
                           self.ScanParam.get('Dec',flag='None')])
        
        rotAngles = array(self.ScanParam.ParAngle)
        chanListAzEl = array(self.BolometerArray.UsedChannels)-1
        OffsetsAzEl=array((take(self.BolometerArray.Offsets[0,::],chanListAzEl), \
                           take(self.BolometerArray.Offsets[1,::],chanListAzEl)))
        refChOffsets=array((self.BolometerArray.RefOffX, self.BolometerArray.RefOffY), 'f')
        AXIS1 = array([model.WCS['NAXIS1'],model.WCS['CRPIX1'],
                       model.WCS['CDELT1'],model.WCS['CRVAL1'],1.])
        AXIS2 = array([model.WCS['NAXIS2'],model.WCS['CRPIX2'],
                       model.WCS['CDELT2'],model.WCS['CRVAL2'],1.])
        
        result = fMap.flagsourceold(chanListIndexes, self.Data, model.Data,
                                 threshold, flag,
                                 XYOffsets, OffsetsAzEl, rotAngles, refChOffsets,
                                 AXIS1, AXIS2)
        
        mask = result.astype('1')
        self.FlagHandler.setOnMask(mask, flag)

    #--------------------------------------------------------------------------------

    def addSource(self,model,chanList=[],factor=1.):
        """
        DES: add data to time stream according to a model map
        INP: (i list) chanList: the list of channels to work with
             (f)     factor: multiply by this factor (default 1)
             (Image object) model: the input model map (with WCS)
                            (default: use current data.Map)
        """

        if not model:
            model = self.Map

        if not model.Data:
            self.MessHand.error("no map computed yet, and no model provided")
            return

        chanList = self.BolometerArray.checkChanList(chanList)
        if len(chanList)<1: 
            self.MessHand.error("no valid channel")
            return
        chanListIndexes = self.BolometerArray.getChanIndex(chanList)

        if string.find(model.WCS['CTYPE1'],'GLON') > -1:
            if not len(self.ScanParam.GalAngle):
                if not len(self.ScanParam.GLon):
                    self.ScanParam.computeGal()
                self.ScanParam.computeGalAngle()
            rotAngles = array(self.ScanParam.ParAngle) + array(self.ScanParam.GalAngle)
            XYOffsets = array([self.ScanParam.get('Glon',flag='None'), \
                               self.ScanParam.get('Glat',flag='None')])
        else:
            rotAngles = array(self.ScanParam.ParAngle)
            XYOffsets = array([self.ScanParam.get('RA',flag='None'), \
                               self.ScanParam.get('Dec',flag='None')])

        chanListAzEl = array(self.BolometerArray.UsedChannels)-1
        OffsetsAzEl=array((take(self.BolometerArray.Offsets[0,::],chanListAzEl), \
                           take(self.BolometerArray.Offsets[1,::],chanListAzEl)))
        refChOffsets=array((self.BolometerArray.RefOffX, self.BolometerArray.RefOffY), 'f')
        AXIS1 = array([model.WCS['NAXIS1'],model.WCS['CRPIX1'],
                       model.WCS['CDELT1'],model.WCS['CRVAL1'],1.])
        AXIS2 = array([model.WCS['NAXIS2'],model.WCS['CRPIX2'],
                       model.WCS['CDELT2'],model.WCS['CRVAL2'],1.])

        # get the new data + factor x model array
        tmp = fMap.addsource(chanListIndexes, self.Data, model.Data, \
                             XYOffsets, OffsetsAzEl, rotAngles, refChOffsets, \
                             AXIS1, AXIS2, factor)
        # replace self.Data with updated one
        self.Data = copy.copy(tmp)
        self._DataAna__resetStatistics()
        tmp = 0  # free memory        
    #--------------------------------------------------------------------------------
    
    def zoom(self,mouse=1,style='idl4',wedge=1,\
                limitsZ=[],aspect=0,limitsX=[],limitsY=[],caption=None,\
                doContour=0,levels=[],showRms=1,rmsKappa=3.5):
        """
        DES: allow the user to select a region in the map to zoom in
        INP: (bool) mouse: use the mouse? (default: yes)
               (other parameters: same as showMap)
        """

	if not caption:
            caption = self.ScanParam.caption()
            
        self.Map.zoom(mouse=mouse,style=style,
                      limitsX=limitsX,limitsY=limitsY,limitsZ=limitsZ,
                      aspect=aspect,caption=caption,wedge=wedge,
                      doContour=doContour,levels=levels,
                      showRms=showRms, rmsKappa=rmsKappa)
        
    #--------------------------------------------------------------------------------
    def beamMap(self,chanList=[], \
                channelFlag=[], plotFlaggedChannels=0, \
                dataFlag=[], plotFlaggedData=0, \
                oversamp=2.0,sizeX=[],sizeY=[],\
		style='idl4'):
        """
        DES: build a beam map in (Az,El) coordinates
        INP: (int list) chanList = channels to consider
             (integer list) channelFlag : plot data from channels flagged or unflagged accordingly
             (log)  plotFlaggedChannels : channelFlag revers to flagged/unflagged data
             (integer list)    dataFlag : plot data flagged or unflagged accordingly
             (log)      plotFlaggedData : dataFlag revers to flagged/unflagged data
                                          flag   | plotFlagged | Plot..
                                          'None' |  0          | all data
                                          []     |  0          | unflagged data (default)
                                          []     |  1          | data with at least one flag set
                                          1      |  0          | data with flag 1 not set
                                          1      |  1          | data with flag 1 set
                                          [1,2]  |  0          | data with neither flag 1 nor flag 2 set
                                          [1,2]  |  1          | data with either flag 1 or flag 2 set
             (float) oversamp = oversampling factor (beam fwhm / pixel size). Default=2.
             (list float) sizeX = limits in Az of the map
             (list float) sizeY = limits in El of the map
        """

        # do exactly as a normal map, but don't add channel offsets
        self.doMap(chanList=chanList, \
                   channelFlag=channelFlag, plotFlaggedChannels=plotFlaggedChannels, \
                   dataFlag=dataFlag, plotFlaggedData=plotFlaggedData, \
                   oversamp=oversamp,\
		   beammap=1,sizeX=sizeX,sizeY=sizeY,style=style)


    #----------------------------------------------------------------------------
    def chanMap(self,chanList=[], \
                channelFlag=[], plotFlaggedChannels=0, \
                dataFlag=[], plotFlaggedData=0, \
                oversamp=1.,sizeX=[],sizeY=[],\
                style='idl4',limitsZ=[],center=0,showRms=0,rmsKappa=3.5):
        """
        DES: Compute and plot channel maps in HO offset coordinates
        INP: (i list) chanList = channels to consider
             (integer list) channelFlag : plot data from channels flagged or unflagged accordingly
             (log)  plotFlaggedChannels : channelFlag revers to flagged/unflagged data
             (integer list)    dataFlag : plot data flagged or unflagged accordingly
             (log)      plotFlaggedData : dataFlag revers to flagged/unflagged data
                                          flag   | plotFlagged | Plot..
                                          'None' |  0          | all data
                                          []     |  0          | unflagged data (default)
                                          []     |  1          | data with at least one flag set
                                          1      |  0          | data with flag 1 not set
                                          1      |  1          | data with flag 1 set
                                          [1,2]  |  0          | data with neither flag 1 nor flag 2 set
                                          [1,2]  |  1          | data with either flag 1 or flag 2 set
             (float)  oversamp = oversampling factor (beam fwhm / pixel size). Default=2.
             (2xfloat)   sizeX = limits in Az of the map
             (2xfloat)   sizeY = limits in El of the map
             (str)       style = color table to use in images
             (logical)  center = if set, it will shift each map by the bolometer offsets.
	                Thereby it shifts the source to the center of each channel map.
             (logical) showRms = compute and print rms/beam? (default: no)
             (float)  rmsKappa = kappa in kappa-sigma clipping used to compute rms
        """

        # Check the list of input channels
        chanList = self.BolometerArray.checkChanList(chanList, \
                                                     flag=channelFlag,getFlagged=plotFlaggedChannels)

        if not(len(chanList)):
            self.MessHand.error("no valid channel")
            return

	chanListIndexes = self.BolometerArray.getChanIndex(chanList)

	data        = self.Data
	dataWeights = self.DataWeights
	AzElOffsets = array([self.ScanParam.get('AzimuthOffset',flag='None'), \
			     self.ScanParam.get('ElevationOffset',flag='None')])

     
	OffsetsUsed = array(self.BolometerArray.getChanSep(self.BolometerArray.UsedChannels))
	#OffsetsUsed     = array([zeros(self.BolometerArray.NUsedChannels,Float32),\
        #                         zeros(self.BolometerArray.NUsedChannels,Float32)],Float32)

        if center:
            OffsetsUsed = OffsetsUsed*0

        # Retrieve the beam
        fwhm = self.BolometerArray.BeamSize

        # Compute the pixel size and ...
        pixelSize = fwhm / oversamp
        
        # compute roughly the extrema of the map (taking into account
        # the ScanParam flags but not the BolometerArray Flags or the
        # Data Flags
	
	AzOffsets = self.ScanParam.get('AzimuthOffset')
	ElOffsets = self.ScanParam.get('ElevationOffset')

	minmaxAzEl = concatenate(\
		(fStat.minmax(AzOffsets)+array([-1,1])*2*pixelSize - \
		 fStat.minmax(take(OffsetsUsed[0,::],chanListIndexes))[::-1], \
		 fStat.minmax(ElOffsets)+array([-1,1])*2*pixelSize - \
		 fStat.minmax(take(OffsetsUsed[1,::],chanListIndexes))[::-1] ) \
		)

	# Create a master image
        Map = Image()

        Map.Header['Object']    = self.ScanParam.Object
        Map.Header['Telescope'] = self.BolometerArray.Telescope.Name
        Map.Header['FeBe']      = self.BolometerArray.FeBe

        Map.WCS['CTYPE1']   = 'OLON--SFL'
        Map.WCS['CTYPE2']   = 'OLAT--SFL'
        Map.WCS['CUNIT1']   = 'arcsec'
        Map.WCS['CUNIT2']   = 'arcsec'
        Map.WCS['MJD-OBS']  = self.ScanParam.DateObs
        
        # update main WCS keywords
        Map.computeWCS(pixelSize,sizeX,sizeY,minmaxAzEl)
        cdeltUnit = 1./3600.  # CDELTi are in arcsec
	AXIS1 = array([Map.WCS['NAXIS1'],Map.WCS['CRPIX1'],Map.WCS['CDELT1'],
                       Map.WCS['CRVAL1'],cdeltUnit])
	AXIS2 = array([Map.WCS['NAXIS2'],Map.WCS['CRPIX2'],Map.WCS['CDELT2'],
                       Map.WCS['CRVAL2'],cdeltUnit])

        self.MessHand.info("Building maps with dimensions (x,y) = "+str(int(AXIS1[0]))+','+str(int(AXIS2[0])))

        # Initialise a list of maps
        allMaps = []
        
	ChanRef = self.BolometerArray.RefChannel
      
        # Now start a loop on channels, to build each map
        for chan in range(len(chanList)):
            mapData     = zeros((int(AXIS1[0]),int(AXIS2[0])),typecode=Float32)
            mapWeight   = zeros((int(AXIS1[0]),int(AXIS2[0])),typecode=Float32)
            mapCoverage = zeros((int(AXIS1[0]),int(AXIS2[0])),typecode=Float32)
            
            dataFlag = self._removeReservedFlagValues(dataFlag, self.rflags['CHANNEL FLAGGED'])
            if dataFlag==None:
                self.MessHand.error("no valid flags")
                return

            if dataFlag in ['','None']:
                dataFlagMask = ones(shape=self.FlagHandler.getFlags().shape, \
                                    typecode=self.FlagHandler.getFlags().typecode())
            else:
                if plotFlaggedData:
                    dataFlagMask = self.FlagHandler.isSetMask(dataFlag)
                else:
                    dataFlagMask = self.FlagHandler.isUnsetMask(dataFlag)

            mapData, mapWeight, mapCoverage = fMap.horizontalprojection(chanListIndexes[chan], data, \
                                                                        dataWeights, \
                                                                        dataFlagMask, 1, \
                                                                        AzElOffsets, OffsetsUsed,
                                                                        AXIS1, AXIS2,
                                                                        mapData, mapWeight, mapCoverage)
            # Storing of the map
            allMaps.append(copy.copy(mapData))
            
        # Store all maps in Results attribute
        self.__Results = allMaps
    
        labelX="\gD Az ['']"
        labelY="\gD El ['']"
        
        MultiPlot.draw(chanList,allMaps,WCS=Map.WCS,\
		       wedge=1,labelX=labelX, labelY=labelY,\
		       caption=self.ScanParam.caption(),style=style,\
		       limitsZ=limitsZ,nan=1)
        if showRms:
            cell = oversamp
            nbX = int(AXIS1[0]) - cell
            nbY = int(AXIS2[0]) - cell
            for chan in range(len(chanList)):
                oneMap = allMaps[chan]
                # Distribution of flux/beam
                m,nbOk = fStat.meandistribution(oneMap,nbX,nbY,nbX*nbY,cell)
                m = m[:nbOk]
                # clip data at more than 2-sigma (e.g. a source)
                m2,nbOk = fStat.clipping(m,rmsKappa)
                m2 = m2[:nbOk]
                oneRms = fStat.f_stat(m2)[1]
                self.MessHand.info("Channel %i : r.m.s. / beam = %f"%(chanList[chan],oneRms))

    #----------------------------------------------------------------------------
    def flipOffsets(self,system='eq'):
        """
        DES: change sign of telescope offsets w.r.t. reference position
        INP: (string) system = 'eq' or 'ho', to flip RA/Dec offsets or Az/El
                               offsets (default: 'eq')
        """
        self.ScanParam.flipOffsets(system=system)
        
    #----------------------------------------------------------------------------
    def plotBoloRms(self,smoothFactor=1.5,style='idl4',limitsX=[],limitsY=[],limitsZ=[],
                    caption='',noerase=0):
        """
        DES: plot the array with color scale showing rms
        INP: (float) smoothFactor: the map is smooted by this factor x beam
             style, limits? : plot parameters
        """
        chanList = self.BolometerArray.checkChanList([])
        if len(chanList)<1: 
            self.MessHand.error("no valid channel")
            return

        if not self._DataAna__statisticsDone:
            self._DataAna__statistics()

        bolo = self.BolometerArray
        rms = self.getChanListData('rms',chanList)

        # Create a Map object to produce the plot
        mapRms = Map()
        mapRms.BolometerArray = bolo
        nbInt  = 13
        nbBolo = bolo.NUsedChannels
        mapRms.ScanParam.NInt = nbInt
        mapRms.ScanParam.NObs = 1
        b = bolo.BeamSize * 2. / 3600.
        mapRms.ScanParam.AzOff = array([0.,0.,0.   ,0.  ,0.,-b,-b/2.,b/2.,b ,b/2.,-b/2.,-b/2.,b/2.],'f')
        mapRms.ScanParam.ElOff = array([0.,-b,-b/2.,b/2.,b ,0.,0.   ,0.  ,0.,b/2.,b/2.,-b/2.,-b/2.],'f')
        mapRms.Data        = ones((nbInt,nbBolo),'f')
        mapRms.DataWeights = ones((nbInt,nbBolo),'f')
        mapRms.FlagHandler = BoaFlagHandler.createFlagHandler(zeros((nbInt,nbBolo),Int8))
        mapRms.ScanParam.FlagHandler = BoaFlagHandler.createFlagHandler(zeros(nbInt,Int32))
        for i in xrange(len(rms)):
            mapRms.Data[:,i] *= array(rms[i],'f')

        mapRms.doMap(oversamp=1,noPlot=1,showRms=0)
        mapRms.Map.smoothBy(smoothFactor*bolo.BeamSize)
        # if no limitsZ input, try to compute a good one
        if limitsZ == []:
            med = fStat.f_median(rms)
            sigma = fStat.f_rms(rms,med)
            z1 = max(min(rms),med-sigma)
            z2 = max(max(rms),med+sigma)
            limitsZ = [z1,z2]
        if not caption:
            caption = bolo.FeBe + ' - channel rms'
        mapRms.Map.display(style=style,limitsX=limitsX,limitsY=limitsY,
                           limitsZ=limitsZ,caption=caption,noerase=noerase)
        bolo.plotArray(overplot=1,num=1)
        
    #--------------------------------------------------------------------------------
    def reduce(self,datasetName='',obstoProc=[],update=0,febe='',tau=0.):
        """
        DES: Process a map scan - this method is called by the apexCalibrator
        INP: (string) datasetName: path to the dataset to be reduced
             (i list) obstoProc: list of subscans to consider (default: all)
        """

        if len(obstoProc)==1:
            if type(obstoProc[0]) == type([]): # e.g. obstoProc == [range(4,8)]
                self.read(inFile=datasetName,subscans=obstoProc[0],febe=febe)
            else:
                # cannot work subscan by subscan
                self.read(inFile=datasetName,subscans=range(1,obstoProc[0]+1),febe=febe)
        else:
            self.read(inFile=datasetName,subscans=obstoProc,febe=febe)

        # Automatic flagging of dead channels
        self.flagFractionRms(ratio=5.)   # flag at median(rms)/5 and median*5

        # Conversion to Jy - ToDo: not so clean!
        if string.find(self.BolometerArray.FeBe,'LABOCA') >= 0:
            be = self.BolometerArray.BEGain
            self.Data *= array(be/270. * 6.3E6,'f')
        elif string.find(self.BolometerArray.FeBe,'BOLOSZ') >= 0:
            self.Data *= array(0.135,'f')
        else:
            self.MessHand.warning("Unknown instrument - data not calibrated")

        if tau:
            self.correctOpacity(tau)

        # Remove correlated noise
        self.zeroStart()
        self.medianNoiseRemoval(chanRef=-1,factor=0.99,nbloop=2)
        self.medianNoiseRemoval(chanRef=-2,factor=0.99,nbloop=1)
        # Another order 1 baseline
        self.polynomialBaseline(order=1,subscan=0)
        # despiking
        self.despike(below=-3,above=10)

        # compute weight - hopefully, no too strong source!
        self.computeWeight()
        self.doMap(oversamp=2,system='EQ',noPlot=1)
        self.Map.computeRms()
        
        # compute good limitsZ
        pixels = ravel(self.Map.Data)
        pixels,nb = compressNan([pixels])
        minmax = fStat.minmax(pixels)
        minZ = max([minmax[0],-3.*self.Map.RmsBeam])

        # clip max according to weights
        maxW = fStat.minmax(ravel(self.Map.Weight))[1]
        minW = sqrt(maxW)
        if minW > maxW:
            # happens if maxW < 1.
            minW = maxW/5.
        ok = where(greater(self.Map.Weight,minW),1,0)
        goodmap = take(ravel(self.Map.Data),nonzero(ravel(ok)))
        goodmap,nb = compressNan([goodmap])
        maxZ = max([max(goodmap),5.*self.Map.RmsBeam])
        self.Map.display(limitsZ=[minZ,maxZ],aspect=0,caption=self.ScanParam.caption())

        # compute rms in good weighted part, after sigma-clipping
        m2,nbOk = fStat.clipping(goodmap,3.)
        m2 = m2[:nbOk]
        rms = fStat.f_stat(m2)[3]  # median deviation
        x0,y0 = self.Map.wcs2phy(self.Map.WCS['NAXIS1']/2.,1)
        Plot.xyout(x0,y0,str("map r.m.s. = %6.1f mJy/beam"%(1.E3*rms)),size=2)
        
#--------------------------------------------------------------------------------
def mapSum2(mapList):
    """
    Function (NOT a method) to co-add Image objects.
    Map data, weights and coverage planes are co-added.
    Returns a new Image object, with same WCS and data size.
    
    WARNING: this function assumes that all Image objects correspond
             to the same region of the sky (same size, same center)

    # Example of use:
    scans   = [some list of scan numbers]
    mapList = []  # initialise empty list
    ra1,ra2,de1,de2 = ...  # define limits to be used for all maps
    for s in scans:
        read(str(s))
        <processing of each scan>
        mapping(system='EQ',sizeX=[ra1,ra2],sizeY=[de1,de2])
        mapList.append(data.Map)
    ms = mapSum(mapList)  # co-added Image object
    ms.display()          # can be displayed
    ms.zoom()             # zoom function can be used
    ms.writeFits("output.fits")
    """

    result = copy.deepcopy(mapList[0])

    result.Data[:,::]=0.0
    result.Weight[:,::]=0.0
    result.Coverage[:,::]=0.0
    for m in mapList:
        try:
            mask_bad=where( (m.Data <= 10000000000000.),0,1)
            Numeric.putmask(m.Data,mask_bad,0.0)
            Numeric.putmask(m.Weight,mask_bad,0.0)
            result.Weight=result.Weight+m.Weight
            result.Coverage=result.Coverage+m.Coverage
            result.Data=result.Data+m.Data*m.Weight
        except:
            print 'map could not be added'
    result.Data=result.Data/result.Weight
    mask_bad=where((result.Weight == 0.0),1,0)
    Numeric.putmask(result.Data,mask_bad,float('nan'))

    return result


#--------------------------------------------------------------------------------
def mapSum(mapList):
    """
    Function (NOT a method) to co-add Image objects.
    Map data, weights and coverage planes are co-added.
    Returns a new Image object, with same WCS and data size.
    
    WARNING: this function assumes that all Image objects correspond
             to the same region of the sky (same size, same center)

    # Example of use:
    scans   = [some list of scan numbers]
    mapList = []  # initialise empty list
    ra1,ra2,de1,de2 = ...  # define limits to be used for all maps
    for s in scans:
        read(str(s))
        <processing of each scan>
        mapping(system='EQ',sizeX=[ra1,ra2],sizeY=[de1,de2])
        mapList.append(data.Map)
    ms = mapSum(mapList)  # co-added Image object
    ms.display()          # can be displayed
    ms.zoom()             # zoom function can be used
    ms.writeFits("output.fits")
    """

    result = copy.deepcopy(mapList[0])

    nX = result.WCS['NAXIS1']
    nY = result.WCS['NAXIS2']
    for i in range(nX):
        for j in range(nY):
            Pij = 0.
            Wij = 0.
            Cij = 0.
            for m in mapList:
                if str(m.Data[i,j]) != str(float('nan')):
                    Pij += m.Weight[i,j] * m.Data[i,j]
                    Wij += m.Weight[i,j]
                    Cij += m.Coverage[i,j]
            if Wij:
                result.Data[i,j] = Pij/Wij
                result.Weight[i,j] = Wij
                result.Coverage[i,j] = Cij

    return result

#--------------------------------------------------------------------------------
def mapsumfast(mapList):
    """
    Function (NOT a method) to co-add Image objects.
    Map data, weights and coverage planes are co-added.
    Returns a new Image object, with same WCS and data size.
    
    WARNING: this function assumes that all Image objects correspond
             to the same region of the sky (same size, same center)

    # Example of use:
    scans   = [some list of scan numbers]
    mapList = []  # initialise empty list
    ra1,ra2,de1,de2 = ...  # define limits to be used for all maps
    for s in scans:
        read(str(s))
        <processing of each scan>
        mapping(system='EQ',sizeX=[ra1,ra2],sizeY=[de1,de2])
        mapList.append(data.Map)
    ms = mapSum(mapList)  # co-added Image object
    ms.display()          # can be displayed
    ms.zoom()             # zoom function can be used
    ms.writeFits("output.fits")
    """
    result = copy.deepcopy(mapList[0])
    if len(mapList) == 1:
        return result
    for i in range(len(mapList)-1):
        secondmap=mapList[i+1]
        newmap,newcoverage,newweight=fMap.mapsum(result.Data,result.Weight,result.Coverage,
                                                secondmap.Data,secondmap.Weight,secondmap.Coverage)
        result.Data = newmap
        result.Weight = newweight
        result.Coverage = newcoverage

    return result

#--------------------------------------------------------------------------------

def setValuesPolygon(map, poly=zeros((1,2)), inout='IN', value=0.): 
    """
    DES: function to replace map data inside/outside a polygon with a given value
    INP: (float array) poly : vertices of polygon
         (str)        inout : inside/outside the polygon, one of 'IN' or 'OUT' 
         (float)      value : replace with this value
    OUT: (object)       map : new image object with same wcs and data size
    """
    # check polygon
    if poly.shape[1] != 2: 
        self.MessHand.error("no valid polygon: wrong dimension")
        return
    if poly.shape[0] <= 2:
        self.MessHand.error("no valid polygon: not enough vertices")
        return

    # check inout
    inout = string.upper(inout) 
    if inout not in ['IN','OUT']: 
        self.MessHand.error("Inside or outside the polygon?")
        return   

    result = copy.deepcopy(map)

    if inout=='IN':
        for i in range(result.WCS['NAXIS1']):
            for j in range(result.WCS['NAXIS2']):
                x,y=result.wcs2phy(i,j)   
                if inPolygon(x,y,poly):
                    result.Data[i,j]=value
    elif inout=='OUT':
        for i in range(result.WCS['NAXIS1']):
            for j in range(result.WCS['NAXIS2']):
                x,y=result.wcs2phy(i,j)   
                if outPolygon(x,y,poly):
                    result.Data[i,j]=value    
    return result
