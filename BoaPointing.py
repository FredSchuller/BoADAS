# Copyright (C) 2002-2006
# Max-Planck-Institut fuer Radioastronomie Bonn
# Argelander Institut fuer Astronomie
# Astronomisches Institut der Ruhr-Universitaet Bochum
#
# Produced for the LABOCA project
#
# This library is free software; you can redistribute it and/or modify it under
# the terms of the GNU Library General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option) any
# later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Library General Public License for more
# details.
#
# You should have received a copy of the GNU Library General Public License
# along with this library; if not, write to the Free Software Foundation, Inc.,
# 675 Massachusetts Ave, Cambridge, MA 02139, USA.
#

"""
NAM: BoaPoint.py (module)
DES: contains the BoA Pointing reduction tools
"""
__version__=  '$Revision: 2794 $'
__date__=     '$Date: 2015-03-02 15:30:03 +0100 (Mon, 02 Mar 2015) $'

#----------------------------------------------------------------------------------
#----- Import ---------------------------------------------------------------------
#----------------------------------------------------------------------------------
from boa import BoaMapping
from Numeric import *
import string

from boa           import BoaConfig
from boa.Bogli     import Plot, Forms, BogliConfig
from boa.Utilities import array2list, fitBaseEllipticalGaussian, tolist_boa
from boa.Utilities import compressNan, fitGaussian, modelgaussbase
from boa.BoaError  import BoaError
from boa.fortran   import fUtilities, fMap, fStat

#----------------------------------------------------------------------------------
#----- Boa Pointing Class ----------------------------------------------------------
#----------------------------------------------------------------------------------
class Point(BoaMapping.Map):
    """
    NAM: Point (class)
    DES: An object of this class is responsible for the reduction of
         pointing scan(s)
    """

  
    def __init__(self):
        """
        DES: Initialise an instance.
        """
        BoaMapping.Map.__init__(self)

        self.Maps = []
        self.__Result = []
        self.PointingResult = 0

#--------------------------------------------------------------------------------
#----- public methods -----------------------------------------------------------
#--------------------------------------------------------------------------------

    def iterMap(self,chanList=[],phase=0,flag=0,sizeX=[],sizeY=[]):
        """
        DES: reconstruct a map in (Az,El) coordinates combining bolometers
             and using varying scale to zoom on signal
        INP: (int list) chanList = channels to consider
             (int) phase = phase to plot
             (int) flag  = flag values to consider
             (list float) sizeX = limits in Az of the map
             (list float) sizeY = limits in El of the map
        """

        self.Maps = []  # reinitialise Maps attribute
        # start with one pixel per beam
	self.slowMap(chanList=chanList,phase=phase,flag=flag,sizeX=sizeX,\
		sizeY=sizeY,oversamp=1.)
       
        firstMap = self._Map__Results
	self.Maps.append(firstMap)

	# localize max. of signal
	dimX = shape(firstMap)[0]
	dimY = shape(firstMap)[1]

	# convert to list, to use index
	l = tolist_boa(ravel(firstMap))
        
        # remove NaN values, if any
        allNumbers = []
        for x in l:
              if str(x) != str(float('Nan')):
                    allNumbers.append(x)
                    
        maxi = max(allNumbers)
        mini = min(allNumbers)
        # remove the Nan in this list, otherwise index does not work
        for n in range(len(l)):
            if str(l[n]) == str(float('Nan')):
                   l[n] = mini
	posMax = l.index(maxi)
        
	# convert to x,y
	maxX = int(posMax/dimY)
	maxY = posMax - dimY * maxX + 0.5
        maxX = maxX + 0.5
 	# then to Az, El
        maxAz,maxEl = self.wcs2phy(maxX,maxY)
        # new limits for the second map
        newX1,newY1 = self.wcs2phy(maxX-3.5,maxY-3.5)
        newX2,newY2 = self.wcs2phy(maxX+3.5,maxY+3.5)
        Plot.plot([newX1,newX2,newX2,newX1,newX1],\
                  [newY1,newY1,newY2,newY2,newY1],style='l',overplot=1)
	self.MessHand.info(" after 1st iteration, offsets Az,El ="+str(maxAz)+","+str(maxEl))

	# second iteration:
	# compute a second map +/- 3 beams around max
        self.slowMap(chanList=chanList,phase=phase,flag=flag,\
		sizeX=[newX1,newX2],sizeY=[newY1,newY2],oversamp=3.)
        secondMap = self._Map__Results
	self.Maps.append(secondMap)
	# determine new maximum position
	dimX = shape(secondMap)[0]
	dimY = shape(secondMap)[1]
	# convert to list, to use index
	l = tolist_boa(ravel(secondMap))
        # remove NaN values, if any
        allNumbers = []
        for x in l:
              if str(x) != str(float('Nan')):
                    allNumbers.append(x)
        maxi = max(allNumbers)
        mini = min(allNumbers)
        # remove the Nan in this list, otherwise index does not work
        for n in range(len(l)):
            if str(l[n]) == str(float('Nan')):
                   l[n] = mini
	posMax = l.index(maxi)
	# convert to x,y
	maxX = int(posMax/dimY)
	maxY = posMax - dimY * maxX + 0.5
        maxX = maxX + 0.5
	# then to Az, El
        maxAz,maxEl = self.wcs2phy(maxX,maxY)
        # new limits for the third map
        newX1,newY1 = self.wcs2phy(maxX-3.5,maxY-3.5)
        newX2,newY2 = self.wcs2phy(maxX+3.5,maxY+3.5)
        Plot.plot([newX1,newX2,newX2,newX1,newX1],\
                  [newY1,newY1,newY2,newY2,newY1],style='l',overplot=1)
	self.MessHand.info(" after 2nd iteration, offsets Az,El ="+str(maxAz)+","+str(maxEl))
		
	# third iteration:
	# compute a second map +/- 3 pixels around max
        self.slowMap(chanList=chanList,phase=phase,flag=flag,\
		sizeX=[newX1,newX2],sizeY=[newY1,newY2],oversamp=5.)

	self.Maps.append(self._Map__Results)
		
    def solvePointing(self,chanList=[], gradient=0, circular=0, radius= -5, \
                      Xpos = 0., Ypos = 0., fixedPos = 0, \
                      plot=0, display=1, caption='', aspect=1):
        """
        DES: compute the offset
        INP: (int list) chanList: list of channels to be used (default: all)
             (boolean)  gradient: shall we fit a gradient ? (default: no)
             (boolean)  circular: fit a cricular gaussian instead of an elliptical gaussian
             (float)    radius  : use only bolo inside this radius (negative means multiple of beam) (default 5 beams)
             (float)  Xpos,Ypos : source position if using fixed position
             (boolean) fixedPos : if set, don't fit position, but use Xpos, Ypos
             (boolean)  plot    : do we plot the results? (default: no)
             (boolean)  display    : display the result of the fit (default: yes)    
        OUT: store in self.PoitingResult the results of the fit (i.e. all parameters
             as computed by mpfit routine). If mpfit failed, then self.PoitingResult
             is set to -1
        """

        # Retrieve the data..
        # The az/el off already contains the actual array offsets 
        chanList = self.BolometerArray.checkChanList(chanList)
        theData = self.getChanListData('flux',chanList)
        theWeig = self.getChanListData('weight',chanList)
        dataX = self.getChanListData('azoff',chanList)
        dataY = self.getChanListData('eloff',chanList)

        fwhm = self.BolometerArray.BeamSize

        # Transform that to a 1D array
        dataX = array2list(dataX)
        dataY = array2list(dataY)
        theData = array2list(theData)
        theWeig = array2list(theWeig)

        if radius < 0:
            radius = abs(radius)*fwhm
            
        if radius != 0:
            
            good_for_fit = less(sqrt((dataX-Xpos)**2+(dataY-Ypos)**2),radius)

            dataX   = compress(equal(good_for_fit,1),dataX)
            dataY   = compress(equal(good_for_fit,1),dataY)
            theData = compress(equal(good_for_fit,1),theData)
            theWeig = compress(equal(good_for_fit,1),theWeig)

            # fUtilities.compress crash for array > 2100000
            # dataX,n   = fUtilities.compress(dataX,good_for_fit,1)
            # dataX     = dataX[:n]
            # dataY,n   = fUtilities.compress(dataY,good_for_fit,1)
            # dataY     = dataY[:n]
            # theData,n = fUtilities.compress(theData,good_for_fit,1)
            # theData   = theData[:n]
          
            
        # now use the data to fit a 2D-Gaussian
        
        try:
            self.PointingResult = fitBaseEllipticalGaussian(theData,dataX,dataY,\
                                                            err = 1./sqrt(theWeig),\
                                                            fwhm=fwhm,\
                                                            gradient=gradient,\
                                                            circular=circular,\
                                                            Xpos=Xpos,Ypos=Ypos,\
                                                            fixedPos=fixedPos)
            self.showPointing(plot=plot,display=display,caption=caption,aspect=aspect)
        except BoaError, error:
            self.MessHand.warning('fit did not converge:')
            self.MessHand.warning(error.msg)
            self.PointingResult = -1

    def solvePointingOnMap(self,gradient=0, circular=0, radius=-10, \
                           Xpos = 0., Ypos = 0., fixedPos = 0, \
                           plot=0, display=1, caption='', aspect=1, style='idl4'):
        """
        DES: compute the offset on the data.Map object
        INP: (boolean)  gradient: shall we fit a gradient ? (default: no)
             (boolean)  circular: fit a cricular gaussian instead of an elliptical gaussian
             (float)    radius  : use only bolo inside this radiu
                                  (negative means multiple of beam - default: 10 beams)
             (float)  Xpos,Ypos : source position if using fixed position
             (boolean) fixedPos : if set, don't fit position, but use Xpos, Ypos
             (boolean)  plot    : do we plot the results? (default: no)
             (boolean)  display : display the result of the fit (default: yes)    
        OUT: store in self.PointingResult the results of the fit (i.e. all parameters
             as computed by mpfit routine). If mpfit failed, then self.PoitingResult
             is set to -1

             WARNING : No Smoothing should be applied to the map
             before using this function, or the fitted fwhm will be
             useless, use fine oversamp to make reasonable fit
        """

        if not self.Map.Data:
            self.MessHand.error('No map computed yet')
            return
        
        # fit a 2D-Gaussian on the map
        try:
            self.PointingResult = self.Map.extractSource(gradient=gradient,
                                                         circular=circular,
                                                         radius=radius,
                                                         Xpos=Xpos,Ypos=Ypos,
                                                         fixedPos=fixedPos)
            self.showPointing(plot=plot,display=display,caption=caption,
                              aspect=aspect,style=style)
        except BoaError, error:
            self.MessHand.warning('fit did not converge:')
            self.MessHand.warning(error.msg)
            self.PointingResult = -1

      
    def showPointing(self,plot=1,display=1,noMap=0,caption='',
                     aspect=1,style='idl4',limitsZ=[],noerase=0):
        """
        DES: display results of last solvePointing (in text, and on the map if plot=1)
        INP: (logical) plot  : display the results on a map (default: no)
             (logical) display : display the result on screen (default: yes)
        """

        if not(self.PointingResult):
            self.MessHand.error('No pointing results to be displayed')
            return

        PointingResult = self.PointingResult

        if display:
            # Display out the result
            
            # Dictionnary to make links between keys in the pointing result
            # variable and the displayed name            
            pointingDict = {'gauss_x_offset': 'Delta Az ["]',\
                            'gauss_x_fwhm'  : 'FWHM_1   ["]',\
                            'gauss_y_offset': 'Delta El ["]',\
                            'gauss_y_fwhm'  : 'FWHM_2   ["]',\
                            'gauss_tilt'    : 'Tilt   [deg]',
                            'gauss_peak'    : 'Peak flux   '}
            
            if self.Map.WCS['CTYPE1'].find('RA') >= 0 and self.Map.WCS['CTYPE2'].find('DEC') >= 0:
                pointingDict['gauss_x_offset'] = 'RA     [deg]'
                pointingDict['gauss_y_offset'] = 'Dec    [deg]'
                # if map is equatorial, convert FWHM to arcsec for display
                PointingResult['gauss_x_fwhm']['value'] *= 3600.
                PointingResult['gauss_y_fwhm']['value'] *= 3600.
            if self.Map.WCS['CTYPE1'].find('GLON') >= 0 and self.Map.WCS['CTYPE2'].find('GLAT') >= 0:
                pointingDict['gauss_x_offset'] = 'GLon   [deg]'
                pointingDict['gauss_y_offset'] = 'GLat   [deg]'
                # if map is in Galactic, convert FWHM to arcsec for display
                PointingResult['gauss_x_fwhm']['value'] *= 3600.
                PointingResult['gauss_y_fwhm']['value'] *= 3600.
                
            # Select what will be displayed from the PointingResult
            toBeOutput = ['gauss_peak',\
                          'gauss_x_offset','gauss_y_offset',\
                          'gauss_x_fwhm','gauss_y_fwhm','gauss_tilt']
            
            outStr=""
            for item in toBeOutput:
                outStr+=pointingDict[item] + " = " + \
                         str(PointingResult[item]['value'])+'\n'
            self.MessHand.info(outStr)
            if self.Map.WCS['CTYPE1'].find('RA') >= 0 and self.Map.WCS['CTYPE2'].find('DEC') >= 0:
                # convert FWHM back to deg to compute the ellipsis
                PointingResult['gauss_x_fwhm']['value'] *= 1./3600.
                PointingResult['gauss_y_fwhm']['value'] *= 1./3600.
            if self.Map.WCS['CTYPE1'].find('GLON') >= 0 and self.Map.WCS['CTYPE2'].find('GLAT') >= 0:
                PointingResult['gauss_x_fwhm']['value'] *= 1./3600.
                PointingResult['gauss_y_fwhm']['value'] *= 1./3600.
        
        if plot:
            # extract parameters of the ellipsis to plot
            xoff=PointingResult['gauss_x_offset']['value']
            yoff=PointingResult['gauss_y_offset']['value']
            xfwhm=PointingResult['gauss_x_fwhm']['value']
            yfwhm=PointingResult['gauss_y_fwhm']['value']
            tilt=PointingResult['gauss_tilt']['value']*pi/180.

            if not self.Map:
                self.doMap(style=style)
            else:
                if not noMap:
                    self.showMap(aspect=aspect,style=style,limitsZ=limitsZ,
                                 noerase=noerase,caption=caption)
            # Overlay the pointing result
            Forms.ellipse(xoff,yoff,xfwhm,yfwhm,tilt)


    def arrayParameters(self,chanList=[], gradient=0, circular=0, radius=0, plot=0):
        """
        DES: determine the array parameters from the data
        INP: (i list) chanList : the channel list to be used (default: current list)
             (l)      gradient : remove a background gradient in the data (default: no)
             (l)      circular : fit a cricular gaussian instead of an elliptical gaussian
        """

        # Check the channel list
        chanList = self.BolometerArray.checkChanList(chanList)

        fwhm = self.BolometerArray.BeamSize
        
        arrayParamOffsets = []

        # Loop through all the channels
        for chan in chanList:
            if plot:
                self.doMap([chan],aspect=1)

            # By construction the source should be at 0,0
            self.solvePointing([chan],gradient=gradient,circular=circular,radius=radius,\
                               plot=plot, display=0)
            self._Map__Result = []
            
            arrayParamOffsets.append({'channel' : chan,\
                                      'result'  : self.PointingResult})

        self.arrayParamOffsets = arrayParamOffsets
        
    def updateArrayParameters(self,filename=None):
        """
        DES: Update the Parameters Offsets with the computed values
        INP: (str) filename : optional output file name
        """

        arrayParamOffsets = self.arrayParamOffsets
        Offsets = self.BolometerArray.Offsets
        refChan = self.BolometerArray.RefChannel
        gains   = self.FFCF_Gain
        
        # Find the reference offset
        refOffset = array([0,0])
        for iParam in arrayParamOffsets:
            if iParam['channel'] == refChan:
                refOffset = array([iParam['result']['gauss_x_offset']['value'],\
                                   iParam['result']['gauss_y_offset']['value']])

        if filename:
            fout = file(filename,'w')
            
        for iParam in arrayParamOffsets:
            chan = iParam['channel']
            Offsets[:,chan-1] = (Offsets[:,chan-1] - array([iParam['result']['gauss_x_offset']['value'],\
                                                            iParam['result']['gauss_y_offset']['value']]) \
                                 + refOffset).astype(Float32)
            oldgain = self.FF_Median[chan-1]
            gains[chan-1] = iParam['result']['gain']
            if filename:
                fout.write("%i %8.4f %8.4f %8.4f %8.4f\n"%(chan,gains[chan-1],oldgain,\
                                                           Offsets[0,chan-1],Offsets[1,chan-1]))
        self.BolometerArray.Offsets = Offsets
        self.FFCF_Gain = gains
        if filename:
            fout.close()
    
    #---------------------------------------------------------------------------------
    def reduce(self,datasetName='',obstoProc=[],febe='',baseband=1,
               radius=-2.,update=0, tau=0.):
        """
        DES: Process a Pointing scan - this method is called by the apexCalibrator
        INP: (string) datasetName: path to the dataset to be reduced
             (i list)   obstoProc: list of subscans to consider (default: all)
             (string)        febe: frontend-backend to consider
             (float)       radius: radius to be used for fitting (def: 2xbeam)
             (logical)     update: continue previous scan? (def: no)
             (float)          tau: zenithal opacity to apply
        """
        if len(obstoProc)==1:
            if type(obstoProc[0]) == type([]): # e.g. obstoProc == [range(4,8)]
                self.read(inFile=datasetName,subscans=obstoProc[0],
                          febe=febe,readAzEl0=1,baseband=baseband)
            else:
                # cannot work subscan by subscan
                self.read(inFile=datasetName,subscans=range(1,obstoProc[0]+1),
                          febe=febe,readAzEl0=1,baseband=baseband)
        else:
            self.read(inFile=datasetName,subscans=obstoProc,
                      febe=febe,readAzEl0=1,baseband=baseband)

        # If chopped data, then compute phase diff
        if self.ScanParam.WobUsed:
            self.ScanParam.computeOnOff()
            self._phaseDiff()
            
        if string.upper(self.ScanParam.Object) != "JUPITER" and \
               string.upper(self.ScanParam.Object) != "MARS" and \
               string.upper(self.ScanParam.Object) != "SATURN":
            self.flagFractionRms(ratio=10)
        self.medianBaseline()
        self.flatfield()
        
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
            
        Plot.panels(2,2)
        Plot.nextpage()  # start a new page
        # First plot = signal before skynoise removal
        four = self.BolometerArray.fourpixels()
        self.signal(four,noerase=1,caption=self.ScanParam.caption()+' - Raw signal')

        # Now clean up the signal
        beam = self.BolometerArray.BeamSize
        if string.upper(self.ScanParam.Object) == "JUPITER":
            self.flagPosition(radius=3.*beam,flag=8)
            self.flagFractionRms(ratio=10)
        else:
            self.flagPosition(radius=2.*beam,flag=8)
            self.flagFractionRms(ratio=5)
        self.medianNoiseRemoval(chanRef=-1,factor=0.95,nbloop=5)
        # Repeat median noise splitting channels in 4 groups
        # (corresponds to amplifier boxes for LABOCA)
        nb_all = self.BolometerArray.NChannels
        if nb_all > 12:  # at least 3 bolo per group
            self.medianNoiseRemoval(range(int(nb_all/4)+1),
                                    chanRef=-1,factor=0.98,nbloop=2)
            self.medianNoiseRemoval(range(int(nb_all/4)+1,2*int(nb_all/4)+1),
                                    chanRef=-1,factor=0.98,nbloop=2)
            self.medianNoiseRemoval(range(2*int(nb_all/4)+1,3*int(nb_all/4)+1),
                                    chanRef=-1,factor=0.98,nbloop=2)
            self.medianNoiseRemoval(range(3*int(nb_all/4)+1,nb_all+1),
                                    chanRef=-1,factor=0.98,nbloop=2)
        self.polynomialBaseline(order=1,subscan=0)
        if string.upper(self.ScanParam.Object) == "JUPITER":
            self.flagFractionRms(ratio=10)
        else:
            self.flagFractionRms(ratio=5)
        
        self.despike(below=-3,above=5)
        self.computeWeight()
        self.unflag(flag=8)
        Plot.nextpage()
        four = self.BolometerArray.fourpixels()
        self.signal(four,noerase=1,caption='Skynoise subtracted')

        # Now compute the map and solve for pointing
        self.doMap(oversamp=3.,noPlot=1,
                   sizeX=[-10.*beam,10.*beam],sizeY=[-10.*beam,10.*beam])
        self.Map.computeRms()
        self.solvePointingOnMap(radius=radius,plot=0,display=0)
        Plot.nextpage()
        # compute good Z limits
        pixels = ravel(self.Map.Data)
        pixels,nb = compressNan([pixels])
        minmax = fStat.minmax(pixels)
        minZ = max([minmax[0],-3.*self.Map.RmsBeam])
        if self.PointingResult != -1:
            peak = self.PointingResult['gauss_peak']['value']
            if peak > 3.*self.Map.RmsBeam:
                maxZ = peak
            else:
                maxZ = minmax[1]
        else:
            maxZ = minmax[1]
        if self.PointingResult != -1:
            self.showPointing(noerase=1,limitsZ=[minZ,maxZ],caption=' ')
            offX = self.PointingResult['gauss_x_offset']['value']
            offY = self.PointingResult['gauss_y_offset']['value']
            fwX  = self.PointingResult['gauss_x_fwhm']['value']
            fwY  = self.PointingResult['gauss_y_fwhm']['value']
            Plot.xyout(0,8.5*beam,str("peak = %6.1f Jy/beam"%(peak)),size=2.)
            Plot.xyout(0,-8.*beam,str("FWHM = %4.1f x %4.1f"%(fwX,fwY)),size=2.)
            Plot.xyout(0,-9.5*beam,str("x=%5.1f ; y=%5.1f"%(offX,offY)),size=2.)
        else:
            # when pointing failed, display only the map
            self.Map.display(noerase=1,limitsZ=[minZ,maxZ],caption=' ')
            Plot.xyout(0,-9.5*beam,"no fit",size=2.5)
        Plot.nextpage()

        # Finally show the RMS of all bolos, and compute median NEFD
        # but flag the source first
        if self.PointingResult != -1:
            if peak > 3.*self.Map.RmsBeam:
                self.flagPosition(Az = offX, El = offY,radius = max([fwX,fwY]))
        self.plotBoloRms(noerase=1)
        rms = self.getChanListData('rms')
        medrms = fStat.f_median(rms)
        sampl = fStat.f_median(1./self.ScanParam.get('deltat'))
        nefd = medrms / sqrt(sampl)
        Plot.xyout(0,1.05*max(self.BolometerArray.Offsets[1,::]),
                   str("NEFD = %5.1f mJy/sqrt(Hz)"%(nefd*1.E3)),size=1.5)
        Plot.panels(1,1)
        
    #---------------------------------------------------------------------------------
    def reduceCross(self,datasetName='',obstoProc=[],febe='',baseband=1,update=0):
        """
        DES: Process a Pointing scan observed with cross-OTF 
        INP: (string) datasetName: path to the dataset to be reduced
             (i list)   obstoProc: list of subscans to consider (default: all)
             (string)        febe: frontend-backend to consider
             (logical)     update: continue previous scan? (def: no)
        """
        if len(obstoProc)==1:
            if type(obstoProc[0]) == type([]): # e.g. obstoProc == [range(4,8)]
                self.read(inFile=datasetName,subscans=obstoProc[0],
                          febe=febe,readAzEl0=1,baseband=baseband)
            else:
                # cannot work subscan by subscan
                self.read(inFile=datasetName,subscans=range(1,obstoProc[0]+1),
                          febe=febe,readAzEl0=1,baseband=baseband)
        else:
            self.read(inFile=datasetName,subscans=obstoProc,
                      febe=febe,readAzEl0=1,baseband=baseband)

        # First, subtract zero order baseline and average noise
        self.medianBaseline(subscan=0)

        # Plot signal before subtracting correlated noise
        Plot.panels(2,2)
        Plot.nextpage()  # start a new page
        four = self.BolometerArray.fourpixels()
        self.signal(four,noerase=1,caption=self.ScanParam.caption()+' - Raw signal')

        # If chopped data, then compute phase diff
        if self.ScanParam.WobUsed:
            self.ScanParam.computeOnOff()
            self._phaseDiff()
            
        self.medianBaseline(subscan=0)
        ref = self.BolometerArray.RefChannel
        self.averageNoiseRemoval(chanRef=ref)

        self.despike(below=-5)  # flag negative, if any
        Plot.nextpage()  # start a new page
        self.signal(four,noerase=1,caption='Phase diffed, average noise subtracted')

        # We will store Az and El offsets based on subscans
        offAz,fluxAz = array([],'f'), array([],'f')
        offEl,fluxEl = array([],'f'), array([],'f')
        for i in range(len(self.ScanParam.SubscanNum)):
            subNum = self.ScanParam.SubscanNum[i]
            tmpAz = self.getChanData('azoff',ref,subscans=[subNum])
            tmpEl = self.getChanData('eloff',ref,subscans=[subNum])
            tmpFlux = self.getChanData('flux',ref,subscans=[subNum])
            # Update Az or El data, depending on SCANDIR
            if string.find(self.ScanParam.ScanDir[i],'LON') > -1:
                self.MessHand.info("Subscan %i in Az direction"%(subNum))
                offAz = concatenate((offAz,tmpAz))
                fluxAz = concatenate((fluxAz,tmpFlux))
            elif string.find(self.ScanParam.ScanDir[i],'LAT') > -1:
                self.MessHand.info("Subscan %i in El direction"%(subNum))
                offEl = concatenate((offEl,tmpEl))
                fluxEl = concatenate((fluxEl,tmpFlux))
            else:
                self.MessHand.warning("Subscan %i in undefined direction, skipping..."%(subNum))
            
        fwhm2sigma = 1./(2*sqrt(2*log(2)))

        if len(offAz):
            err = ones(shape(offAz))
            Plot.nextpage()
            Plot.pointSize(5)
            Plot.plot(offAz,fluxAz,noerase=1,style='l',width=1,
                      labelX='Az offset [arcsec]',labelY='Flux [arb. u.]')
            Plot.plot(offAz,fluxAz,overplot=1,style='p')
            m = fitGaussian(offAz,fluxAz,err,const=1)
            if m.status >= 0:
                resAz = m.params[1]
                ampAz = m.params[0]
                widAz = m.params[2] / fwhm2sigma
                xx = min(offAz) + arange(101)*(max(offAz)-min(offAz))/100.
                yy = modelgaussbase(m.params,xx)
                Plot.plot(xx,yy,overplot=1,ci=2,style='l',width=3)
                Plot.plot([resAz,resAz],[0,m.params[0]],overplot=1,style='l',width=5,ci=3)
                Plot.xyout(resAz,0.,"%5.1f"%(resAz),size=2.5)
                y0 = min(fluxAz) + 0.9*(max(fluxAz)-min(fluxAz))
                Plot.xyout(xx[5],y0,"Az",size=4)
            else:
                resAz,ampAz,widAz = -999,0,0
        else:
            resAz,ampAz,widAz = -999,0,0

        if len(offEl):
            err = ones(shape(offEl))
            Plot.nextpage()
            Plot.pointSize(5)
            Plot.plot(offEl,fluxEl,noerase=1,style='l',width=1,
                      labelX='El offset [arcsec]',labelY='Flux [arb. u.]')
            Plot.plot(offEl,fluxEl,overplot=1,style='p')
            m = fitGaussian(offEl,fluxEl,err,const=1)
            if m.status >= 0:
                resEl = m.params[1]
                ampEl = m.params[0]
                widEl = m.params[2] / fwhm2sigma
                xx = min(offEl) + arange(101)*(max(offEl)-min(offEl))/100.
                yy = modelgaussbase(m.params,xx)
                Plot.plot(xx,yy,overplot=1,ci=2,style='l',width=3)
                Plot.plot([resEl,resEl],[0,m.params[0]],overplot=1,style='l',width=5,ci=3)
                Plot.xyout(resEl,0.,"%5.1f"%(resEl),size=2.5)
                y0 = min(fluxEl) + 0.9*(max(fluxEl)-min(fluxEl))
                Plot.xyout(xx[95],y0,"El",size=4)
            else:
                resEl,ampEl,widEl = -999,0,0
        else:
            resEl,ampEl,widEl = -999,0,0
                
        Plot.panels(1,1)

        # Store results in self.PointingResult
        result = {'gauss_x_offset':resAz,'gauss_y_offset':resEl,
                  'width_x':widAz,'width_y':widEl,
                  'ampl_x':ampAz,'ampl_y':ampEl}
        self.PointingResult = result
                  
        
    #---------------------------------------------------------------------------------

    #---------------------------------------------------------------------------------
    def writeModelData(self):
        """
        Generate one line to be written in the .dat file used for
        determining pointing model
        """
        # NB: this will work only if data have been read in with readAzEl0 = 1
        # and after the pointing has been reduced
        param = self.ScanParam
        # we need to find the 1st timestamp that is not flagged
        i0 = 0
        while param.FlagHandler.isSetOnIndex(i0):
            i0 += 1
        sourceAz = param.AntAz0 - param.AzOff[i0]
        sourceEl = param.AntEl0 - param.ElOff[i0]
        encodAz  = param.EncAz0 - param.AzOff[i0]
        encodEl  = param.EncEl0 - param.ElOff[i0]
        PCA      = param.PDeltaCA
        PIE      = param.PDeltaIE
        FCA      = param.FDeltaCA
        FIE      = param.FDeltaIE
        scan     = param.ScanNum
        cCA      = self.PointingResult ['gauss_x_offset']['value'] / 3600.
        cIE      = self.PointingResult ['gauss_y_offset']['value'] / 3600.
        datLine = '%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %i\n' % \
                  (sourceAz, sourceEl, encodAz, encodEl,
                   cCA, cIE, PCA, PIE, FCA, FIE, scan)
        f = file('LABOCA-ABBA_1_2007-05-10.dat','a')
        f.write(datLine)
        f.close

        
