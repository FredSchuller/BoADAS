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
NAM: Utilities.py (module)
DES: contains utilities & methods for fitting functions
"""

__version__ = '$Revision: 2794 $'
__date__  = '$Date: 2015-03-02 15:30:03 +0100 (Mon, 02 Mar 2015) $'

#----------------------------------------------------------------------------------
#----- Import ---------------------------------------------------------------------
#----------------------------------------------------------------------------------
import time, sys, posix, cPickle

from Numeric import *
from LinearAlgebra import *

try:
    import mpfit
except ImportError:
    from boa import mpfit

from boa.fortran  import fUtilities,fStat
from boa.BoaError import BoaError
from boa          import BoaConfig

def ps():
    posix.system('ps -l | grep python | grep -v grep')

#----------------------------------------------------------------------------------
# Timing --------------------------------------------------------------------------
#----------------------------------------------------------------------------------

class Timing:
    """
    NAM: Timing (class)
    DES: easily profile time computation in program
    """

    def __init__(self):
        # the time mark
        self.initTime = time.time()
        self.lastTime = self.initTime
        self.nIter = 0
        
    def __str__(self):
        return "%10.3f seconds" % (time.time()-self.lastTime)

    def getTime(self):
        return (time.time()-self.lastTime)
    
    def setTime(self):
        self.lastTime = time.time()

    def resetTime(self):
        self.initTime = time.time()
        self.lastTime = self.initTime

    def setIter(self,maxiter=1):
        self.nIter = maxiter
        
    def timeLeft(self,iter=0):
        remainTime  = (time.time()-self.initTime)/(iter+1)*(self.nIter-iter)
        remainMin = int(remainTime/60.)
        print str("Time left : %3im %4.1fs"%(remainMin,remainTime-60*remainMin)),"\r",
        sys.stdout.flush()


class ProgressBar:
    """ 
    NAM : progressBar (class)
    DES : Creates a text-based progress bar.
    """

    def __init__(self, minValue = 0, maxValue = 100, totalWidth=80):
        self.progBar = "[]"   # This holds the progress bar string
        self.min = minValue
        self.max = maxValue
        self.span = maxValue - minValue
        self.width = totalWidth
        self.amount = 0       # When amount == max, we are 100% done 
        self.updateAmount(0)  # Build progress bar string

    def updateAmount(self, newAmount = 0):
        if newAmount < self.min: newAmount = self.min
        if newAmount > self.max: newAmount = self.max
        self.amount = newAmount

        # Figure out the new percent done, round to an integer
        diffFromMin = float(self.amount - self.min)
        percentDone = (diffFromMin / float(self.span)) * 100.0
        percentDone = int(round(percentDone))

        # Figure out how many hash bars the percentage should be
        allFull = self.width - 2
        numHashes = (percentDone / 100.0) * allFull
        numHashes = int(round(numHashes))

        # Build a progress bar with an arrow of equal signs; special cases for
        # empty and full
        if numHashes == 0:
            self.progBar = "[>%s]" % (' '*(allFull-1))
        elif numHashes == allFull:
            self.progBar = "[%s]" % ('='*allFull)
        else:
            self.progBar = "[%s>%s]" % ('='*(numHashes-1),
                                        ' '*(allFull-numHashes))

        # figure out where to put the percentage, roughly centered
        percentPlace = (len(self.progBar) / 2) - len(str(percentDone)) 
        percentString = str(percentDone) + "%"

        # slice the percentage into the b1ar
        self.progBar = ''.join([self.progBar[0:percentPlace], percentString,
                                self.progBar[percentPlace+len(percentString)::]
                                ])

    def __str__(self):
        return str(self.progBar)

    def __call__(self, value):
        """ Updates the amount, and writes to stdout. Prints a carriage return
            first, so it will overwrite the current line in stdout."""
        print '\r',
        self.updateAmount(value)
        sys.stdout.write(str(self))
        sys.stdout.flush()



#----------------------------------------------------------------------------------
# Misc Utilities..  ---------------------------------------------------------------
#----------------------------------------------------------------------------------

def safeExp(x):
    """
    NAM: safeExp (function)
    DES: correct a bug in Numeric that raise an expection when
         computing exponential of small numbers, this take a lot of time !
         but faster thant converting to nummarray compute the exp and back to numeric !!
    """
#    return array(numarray.exp(x))

    condition = less(x,-745.)
    tmp = choose(condition,(x,1.))
    return choose(condition,(exp(tmp),0.))

def min2D(a):
    """
    NAM: min2D (function)
    DES: return the minimum value from a list of 1D Numeric arrays
    """
    nb = len(a)          # number of arrays
    mini = min(a[0])     # first minimum
    for i in range(1,nb):
        if min(a[i]) < mini:
            mini = min(a[i])
    return mini
    
def max2D(a):
    """
    NAM: max2D (function)
    DES: return the maximum value from a list of 1D Numeric arrays
    """
    nb = len(a)          # number of arrays
    maxi = max(a[0])     # first minimum
    for i in range(1,nb):
        if max(a[i]) > maxi:
            maxi = max(a[i])
    return maxi

def array2list(a):
    """
    NAM: array2list (function)
    DES: convert a list of 1D arrays to a single 1D array
    """
    
    result = []
    nbRow = len(a)
    for i in range(nbRow):
        result.extend(a[i])
    return array(result)


def tolist_boa(a):
    """
    NAM: tolist_boa (function)
    DES: replacement for tolist (method) from Numeric (which is leaking memory)
    """
    ## note that this function will become obsolete when and if
    ## there is a fix for the memory leak problem in Numeric 

    if len(a.shape)==1:
        l = list(a)
    else:
        l = []
        for sa in a:
            l += [ tolist_boa(sa) ]
    return l
    

#--------------------------------------------------------------------------------
def compressNan(array):
  """
  DES: return an array without nan
  INP: (array)  array : input array
  OUT: (1D array) values of the previous array without Nan
  """
  
  inputArray = concatenate(copy.copy(array))

  nanIndices = []
  nNan = 0
  
  for index in range(len(inputArray)):
    if str(inputArray[index]) == str(float('Nan')):
      nanIndices.append(0)
      nNan +=1
    else:
      nanIndices.append(1)

  return compress(nanIndices,inputArray), nNan

#--------------------------------------------------------------------------------
def compress2d(array,indexes):
    """
    DES: return a 2D sub array based on indexes
    INP: (f)   array : square input array
         (i) indexes : the indexes to take from the array
    """
    return take(take(array,(indexes)),(indexes),1)

#--------------------------------------------------------------------------------
def lCompressNan(array,listArray):
    """
    DES: remove the Nan of an array in a list of array
    INP:    array : test array for the Nan
         (l array): the list of array to compress
    """

    # search for nan
    nanIndices = []
    nNan = 0
  
    for index in range(len(array)):
        if str(array[index]) == str(float('Nan')):
            nanIndices.append(0)
            nNan +=1
        else:
            nanIndices.append(1)

    # construct the output array

    outputArray = []
    for index in range(len(listArray)):
        outputArray.append(compress(nanIndices,listArray[index]))

    return outputArray
#----------------------------------------------------------------------------------
def prettyPrintList(inputList):
    """
    DES:  Pretty print a list avoiding useless entries
    INP: (l) inputList    : the input list, does not need to be sorted
    OUT: (s) outputString : the resulting string 
    """
    theList = list(inputList)
    theList.sort()
    
    startItem   = theList[0]
    currentItem = theList[0]
    
    outputString = ""
        
    for i in range(1,len(theList)):
        if theList[i]-currentItem > 1:
            if currentItem - startItem !=0:
                outputString = outputString + str(startItem) + "-" + str(currentItem) + "; "
            else:
                outputString = outputString + str(startItem) + "; "

            currentItem = theList[i]
            startItem   = theList[i]
        else:
            currentItem = theList[i]
            

    if currentItem - startItem !=0:
        outputString = outputString + str(startItem) + "-" + str(currentItem)
    else:
        outputString = outputString + str(startItem)
            
    
    return outputString
#----------------------------------------------------------------------------------
def stripFitsExtension(filename):
    """
    DES:  Strip any fits extension from a filename
    INP: (s) filename : the input filename
    OUT: (s) output   : the resulting string 
    """
    output = filename

    if output[-5::] == '.fits':
        output = output[:-5]
    if output[-8::] == '.fits.gz':
        output = output[:-8]
    
    return output

#----------------------------------------------------------------------------------
# Fitting Functions ---------------------------------------------------------------
#----------------------------------------------------------------------------------



# ---- Parabola -------------------------------------------------------------------
def fitParabola(x,y,err):
    """
    NAM: fitParabola (method)
    DES: fits parabola to the data using mpfit
    INP: (float) x = x data
         (float) y = y data
    """
    #p=detStartParaParabola(x,y)
    p=[0.0,0.0,0.0]
    parinfo = [{'value':0.,'mpprint':0}]*3 
    for i in range(3): parinfo[i]['value']=p[i]
    fa={'x':x, 'y':y, 'err':err}
    m=mpfit.mpfit(parabola,p,parinfo=parinfo,functkw=fa,quiet=1)
    #if (m.status <= 0):
    #    raise BoaError, str("mpfit failed: %s"%(m.errmsg))
    return m


# ---------------------------------------------------------------------------------
def detStartParaParabola(x,y):
    """
    NAM: defStartParaParabola (method)
    DES: define the proper start parameter to fit a parabola
    INP: (float) x = x data
         (float) y = y data
    """
    p=[0.0,0.0,0.0]; s=[]; t=[]
    for k in range(0, 5): s.append(sum(x**k))
    for k in range(0, 3): t.append(sum(x**k*y))
    d_0=s[0]*t[1]-s[1]*t[0]
    d_1=s[1]*s[2]-s[0]*s[3]
    d_2=s[1]*s[3]-s[2]**2
    d_3=s[1]*t[2]-s[2]*t[1]
    d_4=s[2]*s[3]-s[1]*s[4]
    d_5=s[0]*s[2]-s[1]**2
    p[2]=(d_3*d_5-d_0*d_2)/(d_1*d_2-d_4*d_5)
    p[1]=(d_3+p[2]*d_4)/d_2
    p[0]=(t[0]-p[1]*s[1]-p[2]*s[2])/s[0]
    return p

# ---------------------------------------------------------------------------------
def parabola(p,fjac=None,x=None,y=None,err=None):
    """
    NAM: parabola
    DES: function used by mpfit to fit a parabola
    """
    model=modelparabola(p,x)
    status = 0
    return([status, (y-model)/err]) 


def modelparabola(p,x):
    """
    NAM: modelparabola
    DES: compute a model parabola at position x for a given set of parameters p
    """
    return(p[0]+p[1]*x+p[2]*x**2)


# ---------------------------------------------------------------------------------
def skydip(p,fjac=None,x=None,y=None,err=None):
    """
    NAM: skydip
    DES: signal as a function of elevation, aka skydip
    """
    model=modelSkydip(p,x)
    status = 0
    return([status, (y-model)/err]) 

def modelSkydip(p,x):
    """
    DES: model function for fitting skydip
         full model, with 5 parameters:
         t(x)=(1-opt)*tcabin+opt*tatm*((1-exp(-tauz/sin(pi*x/180.0)))*feff+(1-feff))
         p[0]=opt, p[1]=tcabin, p[2]=tatm, p[3]=tauz, p[4]=feff
    """
    tmp1 = -1.*p[3]/sin(x*pi/180.)
    tmp2 = p[0]*p[2]*((1.-safeExp(tmp1))*p[4] + (1.-p[4]))
    return ((1.-p[0])*p[1] + tmp2)

def fitSkydip(x,y,err,val0,fixT=1):
    """
    DES: fits a skydip signal-elevation function
         only 2 parameters fitted: opt, tauz
    INP: (f array) x = x data
         (f array) y = y data
         (f array) err=errors on y values
         (5xf list) val0 = first guess values, in this order:
                           [coupling, Tcabin, Tatm, tau_z, F_eff]
    """
    parname = ['Coupling','Tcabin','Tatm','tau_z','F_eff']
    p = val0
    parinfo = []
    for i in range(5):
        parinfo.extend([{'parname': parname[i], \
                         'value': p[i], \
                         'fixed': 0, \
                         'limits' : [0.,0.],\
                         'limited': [1,1]}])
    parinfo[0]['limits'] = [0.,1.]
    parinfo[3]['limits'] = [0.,3.]  # don't observe at tau>3 !!
    if fixT:
        parinfo[1]['fixed'] = parinfo[4]['fixed'] = parinfo[2]['fixed'] = 1
    else:
        parinfo[1]['fixed'] = parinfo[4]['fixed'] = 1
    parinfo[1]['limits'] = [200,320]
    parinfo[2]['limits'] = [100,320]
    parinfo[4]['limits'] = [0.,1.]
    
    fa={'x':x, 'y':y, 'err':err}
    m=mpfit.mpfit(skydip,p,parinfo=parinfo,functkw=fa,quiet=1,debug=0)
    if (m.status <= 0):
        raise BoaError, str("mpfit failed: %s"%(m.errmsg))
    return m

# ---- Gaussian -------------------------------------------------------------------
def fitGaussian(x,y,err,const=0):
    """
    DES: fits a Gaussian to the data using mpfit
    INP: (float) x = x data
         (float) y = y data
         (float) err = array with errors on y
         (log) const = should we include a constant term?
    """
    
    p=[1.0,1.0,1.0]
    if const:
        p.append(0.)
        
    # try to guess the parameters quick-and-dirty
    p[0]=max(y)
    weights=y/sum(y)
    p[1]=sum(x*weights)
    cutoff_mask = Numeric.where((array(y) > max(y)/2.),1,0)
    cutoff = Numeric.compress(cutoff_mask, array(y))
    p[2]=max(x)*float(len(cutoff))/float(len(y))
    if const:
        p[3] = fStat.f_median(y)

    parinfo = [{'value':0.,'mpprint':0}]*3 
    for i in range(3):
        parinfo[i]['value']=p[i]
    if const:
        parinfo.append({'value':0.,'mpprint':0})
        parinfo[3]['value']=p[3]
        
    fa={'x':x, 'y':y, 'err':err}
    try:
        if const:
            m=mpfit.mpfit(gaussbase,p,parinfo=parinfo,functkw=fa,quiet=1)
        else:
            m=mpfit.mpfit(gauss,p,parinfo=parinfo,functkw=fa,quiet=1)
        if (m.status <= 0):
            raise BoaError, str("mpfit failed: %s"%(m.errmsg))
        
    except:
        m.params=[0.,1.,1.]
        if const:
            m.append(0.)
        #self.MessHand.warning("could not fit gaussian")
    return m

# ---------------------------------------------------------------------------------
def gauss(p,fjac=None,x=None,y=None,err=None):
    """
    NAM: gauss
    DES: function used by mpfit to fit a gaussian
    """
    model=modelgauss(p,x)
    status = 0
    return([status, (y-model)/err]) 


def modelgauss(p,x):
    """
    NAM: modelgauss
    DES: compute a model gaussian at position x for a given set of parameters p
    """

    norm=p[0]
    mu=p[1]
    sigma=p[2]

    exponent=((x-mu)/sigma)
    
    #return( (norm/(2.*pi*sigma))*safeExp( -(exponent**2)/2.  ) )
    return( norm*safeExp( -(exponent**2)/2.  ) )

def gaussbase(p,fjac=None,x=None,y=None,err=None):
    """
    NAM: gaussbase
    DES: function used by mpfit to fit a gaussian + constant term
    """
    model=modelgaussbase(p,x)
    status = 0
    return([status, (y-model)/err]) 


def modelgaussbase(p,x):
    """
    NAM: modelgaussbase
    DES: compute a model gaussian + constant term at position x
         for a given set of parameters p: [amplitude, mean, sigma, constant]
    """

    norm=p[0]
    mu=p[1]
    sigma=p[2]
    const=p[3]

    exponent=((x-mu)/sigma)
    
    #return( (norm/(2.*pi*sigma))*safeExp( -(exponent**2)/2.  ) )
    return norm*safeExp( -(exponent**2)/2.) + const

#----- 2D Gauss + first order base surface --------------------------------------

def modelBaseEllipticalGaussian(p,position):
    """
    NAM: model2Dgauss
    DES: compute a 2D gaussian defined by the parameter p wihtin the position
    position should be a list of 2 arrays of same dimensions defining the map
    """
    #p.name = ["gauss_cont","gauss_cont_x","gauss_cont_y", \
    #          "gauss_peak","gauss_x_offset","gauss_y_offset", \
    #          "gauss_x_FWHM","gauss_y_FWHM","gauss_tilt"]
    
    x,y = position
    cont,cont_x,cont_y,g_int,x_offset,y_offset,x_fwhm,y_fwhm,tilt = p
    
    fwhm2sigma = 1./(2*sqrt(2*log(2)))
    
    sigma_x = x_fwhm*fwhm2sigma
    sigma_y = y_fwhm*fwhm2sigma
    
    x_x_offset = x-x_offset
    y_y_offset = y-y_offset
    
    cos_tilt = cos(tilt)
    sin_tilt = sin(tilt)
    gauss_int = 2*pi*sigma_x*sigma_y
    
    xp = x_x_offset*cos_tilt-y_y_offset*sin_tilt
    yp = x_x_offset*sin_tilt+y_y_offset*cos_tilt
    U = (xp/sigma_x)**2+(yp/sigma_y)**2
    return cont + cont_x*x + cont_y*y + g_int/gauss_int*safeExp(-U/2)

#    model = fUtilities.modelbaseellipticalgaussian(p,position)
#    return reshape(model,shape(position[0]))

# ---------------------------------------------------------------------------------
def baseEllipticalGaussian(p,fjac=None,x=None,y=None,err=None):
    """
    NAM: baseEllipticalGaussian
    DES: function used by mpfit to fit a 2D gaussian+base
         (5 elmts array) p : parameters of the gaussian (see modelBase2Dgauss)
         (2d array)      x : position of  the pixels on the map
                             "x" = x[0] and "y" = x[1]
         (2d array)      y : the map to fit
                             y.shape should be (len(x[0]),len(x[1]))
    """
    
    model = fUtilities.modelbaseellipticalgaussian(p,x)
    status=0
   
    return([status, ravel((y-model)/err)])

# ---------------------------------------------------------------------------------
def baseCircularGaussian(p,fjac=None,x=None,y=None,err=None):
    """
    NAM: baseCircularGaussian
    DES: function used by mpfit to fit a Circular gaussian+base
         (5 elmts array) p : parameters of the gaussian
         (2d array)      x : position of  the pixels on the map
                             "x" = x[0] and "y" = x[1]
         (2d array)      y : the map to fit
                             y.shape should be (len(x[0]),len(x[1]))
    """

    lp = concatenate([p,[p[-1],0]])
    model = fUtilities.modelbaseellipticalgaussian(lp,x)    
    status=0
   
    return([status, ravel((y-model)/err)])



# ---------------------------------------------------------------------------------
def fitBaseEllipticalGaussian(mapArray, x, y, err=1.0, fwhm=11.0, gradient=1, circular=0,\
                              Xpos=0., Ypos=0., fixedPos=0, incl=0., fixIncl=0):
    """
    NAM: fitBaseEllipticalGaussian (method)
    DES: fits a 2D Gaussian + 1st order base surface
    INP: (arrays) (x,y,mapArray,err)  : the data to fit (arrays of same dimension(s))
         (2els arrays) sizeX/Y : alternative to the x/y array, this limit the size of the
                                 mapArray given as regular gridding between the center
                                 of the two extreme pixels
         (float)         fwhm  : the first guess for the fwhm
         (logical)   gradient  : should we also fit a gradient in the map (default no) ?
         (logical)   circular  : fit a circular gaussian instead of a elliptical gaussian
         (float)     Xpos,Ypos : source position if using fixed position
         (logical)    fixedPos : if set, don't fit position, but use Xpos, Ypos
         (float)         incl  : inclination angle
         (logical)    fixIncl  : fix the inclination angle? default no (0)
         
    OUT: a dictionnary containning the results of the fit
         check 'status' and 'errmsg' to see if the fit was done correctly
         then for each parameters (see the parname variable below) you have the
         'value' 'error' and 'limits' for the fit
    """

    # Test for dimension
    if not ((len(x.shape) == 2 or len(x.shape) == 1) and \
             len(x.shape) == len(y.shape) == len(mapArray.shape)):
        print "Error : Arrays must have the same dimension (1D or 2D)"
        return
    
    # In case we have 2D arrays, put everything into 1D
    x = ravel(x)
    y = ravel(y)
    mapArray = ravel(mapArray)
    err = ravel(err)

    if err[0] == 1.0 and len(err) ==1:
        # Unweighted
        err = ones(mapArray.shape)
        # err = ones(mapArray.shape)*sqrt(abs(mapArray))

    # Search and remove NaN
    lX, lY, lMapArray, lErr = lCompressNan(mapArray,[x,y,mapArray,err])
    # At this point lX/Y/Z/Err are 1D array with only the mesured data points

    # the value to fit
    fa={'x':array([lX,lY]), 'y': lMapArray, 'err':lErr}

    # fitting function by default elliptical gaussian
    fitFunction = baseEllipticalGaussian
        
    # define the parameters to fit
    parname = ["continuum","continuum_x","continuum_y", \
               "gauss_peak","gauss_x_offset","gauss_y_offset", \
               "gauss_x_fwhm","gauss_y_fwhm","gauss_tilt"]
    
    parinfo=[]
    # set the values
    for i in range(len(parname)):
        parinfo.extend([{'parname': parname[i], \
                           'value': 0., \
                           'fixed': 0, \
                         'limits' : [0.,0.],\
                         'limited': [0,0]}])

    # set limits on position
    parinfo[4]['limited'] = [1,1]
    parinfo[4]['limits'] = [min(lX),max(lX)]
    parinfo[4]['value'] = Xpos
    parinfo[5]['limited'] = [1,1]
    parinfo[5]['limits'] = [min(lY),max(lY)]
    parinfo[5]['value'] = Ypos

    # check that X/Y pos are in the limits
    for i in [4,6]:
        if not parinfo[i]['limits'][0] < parinfo[i]['value'] < parinfo[i]['limits'][1]:
            parinfo[i]['value'] = (parinfo[i]['limits'][1]-parinfo[i]['limits'][0])/2

    # fwhm has to be positive so (let say 1/10 of the given fwhm to
    # avoid division by 0) let say also that the fhwm can not be
    # greater than the total map

    # TODO: the definition below is not good: since the gaussian can rotate,
    # gauss_x_FWHM is not bound to the x axis
    parinfo[6]['value']   = fwhm
    parinfo[6]['limited'] = [1,1]
    #parinfo[6]['limits']  = [fwhm/5,max(lX)-min(lX)]
    #parinfo[6]['limits']  = [fwhm/2.,2.*fwhm]
    parinfo[6]['limits']  = [fwhm*0.75,5.*fwhm]  # allow for extended source, but no spike
    parinfo[7]['value']   = fwhm
    parinfo[7]['limited'] = [1,1]
    #parinfo[7]['limits']  = [fwhm/5,max(lY)-min(lY)]
    #parinfo[7]['limits']  = [fwhm/2.,2.*fwhm]
    parinfo[7]['limits']  = [fwhm*0.75,5.*fwhm]  # allow for extended source, but no spike

    # a tilt is always bound in a circle so
    #parinfo[8]['limited'] = [1,1]
    #parinfo[8]['limits'] = [-2.*pi,2.*pi]
    # except that mpfit doesn't converge right...
    parinfo[8]['limited'] = [0,0]
   
    # the peak flux of the source is always positive, this is not SZ !
    parinfo[3]['limited'] = [1,0]
    parinfo[3]['limits']  = [0,0]

    # in case a circular gaussian was asked
    if circular: 
        fitFunction = baseCircularGaussian
        parinfo = parinfo[0:7]

    # if we need to fit a gradient, then first fit it and then retrieve the gaussian
    if gradient:
        # Fix everything for the gaussian
        for i in range(3,len(parinfo)):
            parinfo[i]['fixed'] = 1

        # 0 flux for the gaussian -> pure gradient
        parinfo[3]['value'] = 0

        # set up the first guess array (take the value from the
        # parinfo array)
        p = []
        for i in range(len(parinfo)):
            p.append(parinfo[i]['value'])
        p = array(p)

        result = mpfit.mpfit(fitFunction,p,parinfo=parinfo,functkw=fa,quiet=1,\
                      ftol=1.e-2, xtol=1.e-2, gtol=1.e-2)
        if (result.status <= 0):
            raise BoaError, str("mpfit failed: %s"%(result.errmsg))

        # set the values of the gradient found for the global fit NO !
        # the fitted gradient should be refitted, it is only used to
        # retrive the parameters of the gaussian, otherwise we go into
        # troubles
        for i in range(3):
            parinfo[i]['value'] = result.params[i]

        # compute the residual map after the gradient fit (to find
        # the first guess for the gaussian)

        params = result.params.tolist()
        if circular:
            params.append(result.params[6]) # "gauss_y_fwhm"
            params.append(0)               # "gauss_tilt"
            
        params = array(params)

        model = modelBaseEllipticalGaussian(params,[lX,lY])
        
        # use a variable with a different name otherwhise values in
        # 'fa' will also be changed
        lMapArray = lMapArray-model

        # release everything for the gaussian
        for i in range(3,len(parinfo)):
            parinfo[i]['fixed'] = 0

        # continuum level is biased in case of a weak gradient+strong
        # gaussian so leave it to 0
        
        parinfo[0]['value'] = 0

    # if asked no to fit a gradient, then do not fit it !
    if not gradient:
        for i in range(1,3):
            parinfo[i]['fixed'] = 1
            parinfo[i]['value'] = 0
        
    # We can simply try to retrieve the position of the gaussian:
    maxPos = nonzero(equal(lMapArray,max(lMapArray)))[0]
    # If fixed position, then use it:
    if fixedPos:
         parinfo[4]['value'] = Xpos
         parinfo[4]['fixed'] = 1
         parinfo[5]['value'] = Ypos
         parinfo[5]['fixed'] = 1
    else:
        # Even if not fixed, if initial guesses are given, use them
        if Xpos:
            parinfo[4]['value'] = Xpos
        else:
            parinfo[4]['value'] = lX[maxPos]
        if Ypos:
            parinfo[5]['value'] = Ypos
        else:
            parinfo[5]['value'] = lY[maxPos]

    # fix inclination angle if needed
    if fixIncl:
        parinfo[8]['value'] = incl
        parinfo[8]['fixed'] = 1


    # Suppose we have circular gaussian and the right maximum there,
    # estimate the fwhm
    
#    inside_fwhm = nonzero(greater(lMapArray,(max(lMapArray)-min(lMapArray))/2))
#    distance = sqrt((take(lX,inside_fwhm)-lX[maxPos])**2+\
#                    (take(lY,inside_fwhm)-lY[maxPos])**2)
#    fwhm = max(distance)
    
    # Use it as a first guess
    parinfo[6]['value'] = fwhm

    if not circular:
        parinfo[7]['value'] = fwhm

    # use max in map as first guess for peak flux
    #parinfo[3]['value'] = lMapArray[maxPos]*2*pi*fwhm**2
    parinfo[3]['value'] = lMapArray[maxPos]

    # set up the first guess array (take the value from the parinfo array)
    p = []
    for i in range(len(parinfo)):
        p.append(parinfo[i]['value'])
    p = array(p,Float)
        
    m=mpfit.mpfit(fitFunction,p,parinfo=parinfo,functkw=fa,quiet=1,fastnorm=1,\
                  maxiter=50)# ,\
    #                  ftol=1.e-10, xtol=1.e-9, gtol=1.e-10)
    if (m.status <= 0 or m.status >=4):
        raise BoaError, str("mpfit failed: %s"%(m.errmsg))

    if circular:
        m.params = concatenate([m.params,[m.params[-1],0]])

    result = {'status': m.status, \
              'errmsg': m.errmsg, \
              'params': m.params}

    for i in range(len(parinfo)):
        result[parname[i]] = {'value'  : m.params[i], \
                             'error'   : m.perror[i], \
                             'fixed'   : parinfo[i]['fixed'], \
                             'limits'  : parinfo[i]['limits'],\
                             'limited' : parinfo[i]['limited']}
    if circular:
        result[parname[7]] = { 'value'  : m.params[6], \
                               'error'   : m.perror[6], \
                               'fixed'   : parinfo[6]['fixed'], \
                               'limits'  : parinfo[6]['limits'],\
                               'limited' : parinfo[6]['limited']}
        
        result[parname[8]] = { 'value'  : 0, \
                               'error'   : 0, \
                               'fixed'   : 0, \
                               'limits'  : [0,0],\
                               'limited' : [0,0]}

    # convert tilt angle to degree
    result['gauss_tilt']['value'] *= 180./pi  

    # Make sure that FWHM1 is major axis, and FWHM2 is minor axis...
    if not fixIncl:
        if result['gauss_x_fwhm']['value'] < result['gauss_y_fwhm']['value']:
            result['gauss_x_fwhm'],result['gauss_y_fwhm'] = result['gauss_y_fwhm'],result['gauss_x_fwhm']
            result['gauss_tilt']['value'] += 90.

    # ... and force tilt angle to be within -90 and +90 deg
    ang = result['gauss_tilt']['value']
    result['gauss_tilt']['value'] = (ang+90.)%180. - 90.
    
    # Compute integrated Gaussian
    fwhm2sigma = 1./(2*sqrt(2*log(2)))
    sigma_x = result['gauss_x_fwhm']['value']*fwhm2sigma
    sigma_y = result['gauss_y_fwhm']['value']*fwhm2sigma
    result['gauss_int'] = {'limited' : [0, 0], \
                           'fixed' : 1, \
                           'limits' : [0.0, 0.0], \
                           'value' : result['gauss_peak']['value']*2*pi*sigma_x*sigma_y, \
                           'error' : 0.0}
    
    return result

# ---------------------------------------------------------------------------------

def cropped_circular_gaussian(p,position,threshold=3):
    """
    NAM: cropped_circular_gaussian
    DES: compute a cropped circular gaussian with intensity=1
         defined by the parameter p wihtin the position an a given threshold given in n*'sigma'
         position should be a list of 2 arrays of the same dimension defining the map
    """
    # p.name = ["gauss_x_offset","gauss_y_offset", \
    #           "gauss_fwhm"]

    x,y = position
    x_offset,y_offset,fwhm = p

    sigma_squared = fwhm**2/(8*log(2))

    dist = ((x-x_offset)**2+(y-y_offset)**2)/sigma_squared

    returned_array = zeros(dist.shape)
    
    for i in range(dist.shape[0]):
        if max(dist[i,::]) < threshold:
            returned_array[i,::] = safeExp(-dist[i,::]/2)
            
    return returned_array

def gaussian(r2,sig2):
    """
    DES: Compute value of a Gaussian function
    INP: r2 = _array_ of distances^2, sig2 = sigma^2, related to Gaussian width
    """
    tmp = -1.*r2/(2.*sig2)
    # Put values lesser than -500. to 1., to avoid Overflow (actually
    # underflow) errors with Numeric.exp - the actual limit is -745.1
    #array1 = choose(less(tmp,-500.),(tmp,1.))
    #array1 = exp(array1)
    #array1 = choose(less(tmp,-500.),(array1,0.))
    #return exp(array1) / (2.*pi*sig2)
    return exp(tmp) / (2.*pi*sig2)

def distsq(x1,y1,x2,y2):
    """
    NAM: distsq (function)
    DES: returns distance squared between two points
    INP: (float) x1,y1,x2,y2: coordinates of the two points
    OUT: (float) distance^2
    """
    return (x1-x2)**2 + (y1-y2)**2

def solvePoly(order,dataX,dataY):
    """
    NAM: solvePoly (function)
    DES: perform polyomial interpolation: solve linear system
         dataY = P_n(dataX)
    INP: (int) order : polynomial degree
         (flt arrays) dataX/Y : system to solve
    OUT: (flt array) coeff : polynomial coefficients
    """
    #### TODO: use existing NumPy or Fortran package!!!
    
    result = []
    if order == 1:
      try:
        result.append((dataY[-2]-dataY[-1])/(dataX[-2]-dataX[-1]))
        result.append(dataY[-2]-result[0]*dataX[-2])
      except ZeroDivisionError:
        result = [0.,dataY[-1]]
    #elif order == 2:
      # interpolate parabola...
    # nothing to do if order = 0
    return result
  
#----------------------------------------------------------------------------
# store array attribute of a given class to column major
#----------------------------------------------------------------------------
def as_column_major_storage(classIn):
    """
    DES: save all the attribute as column major to avoid copy in fortran
    """
    attrDict = vars(classIn)
    attrName = attrDict.keys()
    nbAttr = len(attrName)
    
    for attribute in attrName:
        if isinstance(attrDict[attribute],arraytype) and \
               attrDict[attribute] and \
               not fUtilities.has_column_major_storage(attrDict[attribute]) :
            print attribute
            attrDict[attribute] = fUtilities.as_column_major_storage(attrDict[attribute])

#----------------------------------------------------------------------------
# print attribute list of an object
#----------------------------------------------------------------------------
def attrStr(object,badAttributes=[]):
    """
    DES: return a string representing the attributes of the object
    OPT: (str list) badAttributes : list of attributes to remove from the output
    """
    
    attrDic = vars(object)
    attrName = attrDic.keys()
    nbAttr = len(attrName)
    attrName.sort()
    
    # remove badAttributes that can cause trouble
    for badattribute in badAttributes:
        if badattribute in attrName:
            attrName.remove(badattribute)
                    
    out = str(object.__class__)+" object, with "+str(nbAttr)+" attributes\n\n"
    for a in attrName:
        d = attrDic[a]
        typAttr = type(d)
        if (typAttr == type(array(()))):
            dim = shape(d)
            typElement = d.typecode()
        else:
            try:
                dim = len(d)
            except TypeError:
                dim = 0
                
            if (dim and typAttr != dict):
                typElement = type(d[0])
            elif (typAttr == dict):
                typElement = "Keys"
            else:
                typElement = "None"
                        
        out = out+"\t %20s %14s"%(a,typAttr)
        out = out+"\t with %15s elements of type %1s\n"%(dim,typElement)


    return out

        
# ---------------------------------------------------------------------------------
# Utilities related with Fourier transforms
# ---------------------------------------------------------------------------------

def Cr2p(c):
    """
    NAM: Cr2p (function)
    DES: convert complex numbers in rectangular form to polar (mod,arg) form
    INP: (complex)     c : complex number or array
    OUT: [float,float]   : module and phase
    """
    
    amp   = sqrt(c.real**2+c.imag**2)
    phase = arctan2(c.imag,c.real)
    return amp,phase

def Cp2r(amp,phase):
    """
    NAM: Cp2r (function)
    DES: convert complex numbers in polar form to rectangular form (real,imag)
    IN:  [float,float]   : module and phase
    OUT: (complex)     c : complex number or array
    """
    real = amp*cos(phase)
    imag = amp*sin(phase)

    c = zeros(len(amp),Complex)
    c.real = real
    c.imag = imag
    
    return c


# ---------------------------------------------------------------------------------
# Utilities related to PCA
# ---------------------------------------------------------------------------------

def principalComponentAnalysis (rawdata, order):
    """
    NAM: principalComponentAnalysis (function)
    DES: find the principal components of an array
    IN:  (array) data as an NxM array
                 where M - number of channels
                       N - number of time samples
         (order) number of principal components to return
    OUT: data with principal components removed
    """

    # function should eventually be moved to fortran

    m,n=shape(rawdata)
    #m - number of channels
    #n - number of time samples

    rawDataAdjust, originalMean = adjustDataPCA(rawdata)

    # free some memory
    rawdata={}

    # compute corrlation matrix
    corr=fUtilities.matrixmultiply(rawDataAdjust, transpose(rawDataAdjust))

    #eigenvals,eigenvect=LinearAlgebra.eigenvectors(corr)
    eigenvals,eigenvect=eigenvectors(corr)
    corr={}

    #if (order < 0):
    #    # this step should be adaptive, not sure how yet...
    #    order=10

    eigenvals=abs(eigenvals)
    eig_index=argsort(eigenvals)
    
    featureVector=transpose(take(eigenvect,eig_index[0:m-order]))

    # reconstruct the desired feature...
    desiredFeature=fUtilities.matrixmultiply(transpose(featureVector),rawDataAdjust)

    # free some memory
    rawDataAdjust={}

    # ... and reconstruct the data
    rowOrigData=fUtilities.matrixmultiply(featureVector,\
                desiredFeature)+originalMean

    return transpose(rowOrigData), take(eigenvals,eig_index), take(eigenvect,eig_index)

#----------------------------------------------------------------------------
    
def adjustDataPCA (rawdata):
    """
    NAM: adjustDataPCA (function)
    DES: adjust data array by mean for PCA analysis
    IN:  (array) data as an NxM array
                 where M - number of channels
                       N - number of time samples
    OUT: shifted data and mean
    """
    
    m,n=shape(rawdata)
    #m - number of channels
    #n - number of time samples

    # adjust data by mean of each channel
    mean=sum(rawdata,1)/n
    originalMean=transpose(transpose(rawdata*0.0)+mean)
    rawDataAdjust=rawdata-originalMean

    return rawDataAdjust, originalMean

#----------------------------------------------------------------------------
# functions to get tau and calib correction at a given time
def getCalCorr(refmjd,method,calFile):
    """
    NAM: getCalCorr
    DES: get calibration correction factor at a given time from a file
    INP: (f) refmjd: the time (MJD) requested
         (str) method: method used to compute the calibration, can be:
               'linear': linear interpolation between two closest points
               anything else: returns the closest point
         (str) calFile: file name where MJDs and calCorr values are stored
    OUT: (f) returns the calibration correction factor at the required time
    """

    try:
        f = file(calFile)
    except IOError:
        self.__MessHand.error("could not open file %s"%(calFile))
        return

    # read and process CAL file
    param = f.readlines()
    f.close()
    scannumber, date, calmjd, corr, opacitycorr = [], [], [], [], []   # local lists to store MJD and TAU
        
    for i in range(len(param)):	        # -1: skip last line
        if param[i][0] != '!':              # skip comments
            tmp = string.split(param[i])
            scannumber.append(string.atof(tmp[0]))
            calmjd.append(string.atof(tmp[2]))
            corr.append(string.atof(tmp[3]))
            opacitycorr.append(string.atof(tmp[4]))
                
                
    mjd=array(calmjd)
    calcorr=array(corr)
      

    entries=len(mjd)

    mindiff = 1000000.
    mindiffafter = 1000000.
    mindiffbefore = 1000000.
    cbok = 0
    caok = 0
    if method == 'linear':
        for i in range(entries):
            timediff = mjd[i]-refmjd
            if timediff < 0:
                if -1*timediff < mindiffbefore:
                    calbefore = calcorr[i]
                    timebefore = mjd[i]
                    mindiffbefore = -1*timediff
                    cbok = 1
            if timediff > 0:
                if timediff < mindiffafter:
                    calafter = calcorr[i]
                    timeafter = mjd[i]
                    mindiffafter = timediff
                    caok = 1
            if timediff == 0.000:
                resultcalcorr = calcorr[i]
                return resultcalcorr
                
        if cbok == 1 and caok == 1:   
            resultcalcorr = ((calafter-calbefore)/(timeafter-timebefore))*(refmjd-timebefore)+calbefore
        else:
            if cbok == 1:
                resultcalcorr = calbefore
            else:
                resultcalcorr = calafter
                    
    else:
        for i in range(entries):
            timediff = sqrt((mjd[i]-refmjd)**2)
            if timediff < mindiff:
                resultcalcorr = calcorr[i]
                mindiff = timediff

                   
    return resultcalcorr
  
def getTau(refmjd,method,tauFile):
    """
    NAM: getTau
    DES: get zenith opacity (tau) at a given time from a file
    INP: (f) refmjd: the time (MJD) requested
         (str) method: method used to compute the tau, can be:
               'linear': linear interpolation between two closest points
               anything else: returns the closest point
         (str) tauFile: file name where MJDs and tau values are stored
    OUT: (f) returns the tau at the required time
    """

    try:
        f = file(tauFile)
    except IOError:
        self.__MessHand.error("could not open file %s"%(tauFile))
        return

    # read and process TAU file
    param = f.readlines()
    f.close()
    scannumber, date, taumjd, opacity = [], [], [], []   # local lists to store MJD and TAU
    
    for i in range(len(param)):	        # -1: skip last line
        if param[i][0] != '!':              # skip comments
            tmp = string.split(param[i])
            scannumber.append(string.atof(tmp[0]))
            taumjd.append(string.atof(tmp[2]))
            opacity.append(string.atof(tmp[3]))
                
                
    mjd=array(taumjd)
    tau=array(opacity)
      
    entries=len(mjd)

    mindiff = 1000000.
    mindiffafter = 1000000.
    mindiffbefore = 1000000.
    tbok = 0
    taok = 0
    if method == 'linear':
        for i in range(entries):
            timediff = mjd[i]-refmjd
            if timediff < 0:
                if -1*timediff < mindiffbefore:
                    taubefore = tau[i]
                    timebefore = mjd[i]
                    mindiffbefore = -1*timediff
                    tbok = 1
            if timediff > 0:
                if timediff < mindiffafter:
                    tauafter = tau[i]
                    timeafter = mjd[i]
                    mindiffafter = timediff
                    taok = 1
            if timediff == 0.000:
                resulttau = tau[i]
                return resulttau      
        if tbok == 1 and taok == 1:   
            resulttau = ((tauafter-taubefore)/(timeafter-timebefore))*(refmjd-timebefore)+taubefore
        else:
            if tbok == 1:
                resulttau = taubefore
            else:
                resulttau = tauafter
                    
    else:
        for i in range(entries):
            timediff = sqrt((mjd[i]-refmjd)**2)
            if timediff < mindiff:
                resulttau = tau[i]
                mindiff = timediff

                   
    return resulttau
  
def newRestoreData(fileName='BoaData.sav'):
    """ 
    DES: restore a DataEntity object previously saved in a file, and
    set it as the currData attribute of BoaB
    INP: (string) fileName: name of the input file
    optional - default value = 'BoaData.sav'
    """
    #fileName = self.outDir+fileName
    try:
        f = file(fileName)
    except IOError:
        messages.error(" could not open file %s"%(fileName))
        return
    tmp = cPickle.load(f)
    f.close()

    if hasattr(tmp,'DataFlags'):
        tmp.FlagHandler = BoaFlagHandler.createFlagHandler(tmp.DataFlags.astype(Int8))
        tmp.DataFlags = None
    
        tmp.ScanParam.FlagHandler = BoaFlagHandler.createFlagHandler(tmp.ScanParam.Flags.astype(Int32))
        tmp.ScanParam.Flags = None
    
        tmp.BolometerArray.FlagHandler = BoaFlagHandler.createFlagHandler(tmp.BolometerArray.Flags.astype(Int32))
        tmp.BolometerArray.Flags = None

        

    
    #tmp.FillF90()
    return tmp

# ---------------------------------------------------------------------------------
# Utilities related to polygons
# ---------------------------------------------------------------------------------

def inPolygon(x,y,poly):
    """
    DES: check whether point (x,y) is inside a polygon 
    INP: (float)        x/y : coordinates of point  
          (float array) poly : vertices of polygon
    OUT: (int)              : 0=outside, 1=inside polygon
    """
    counter=0
    p1=poly[0]
    for i in range(1,len(poly)+1): 
      p2=poly[i % len(poly)]
      if y > min(p1[1],p2[1]):
        if y <= max(p1[1],p2[1]): 
          if x <= max(p1[0],p2[0]): 
            if p1[1] != p2[1]: 
              xinters=(y-p1[1])*(p1[0]-p2[0])/(p1[1]-p2[1])+p1[0]
              if p1[0] == p2[0] or x <= xinters: 
                counter=counter+1
      p1=p2
    return counter % 2

  #--------------------------------------------------------------------------------

def outPolygon(x,y,poly):
    """
    DES: check whether point (x,y) is outside a polygon 
    INP: (float)        x/y : coordinates of point  
         (float array) poly : vertices of polygon
    OUT: (int)              : 1=outside, 0=inside polygon
    """
    counter=0
    p1=poly[0]
    for i in range(1,len(poly)+1): 
      p2=poly[i % len(poly)]
      if y > min(p1[1],p2[1]):
        if y <= max(p1[1],p2[1]): 
          if x <= max(p1[0],p2[0]): 
            if p1[1] != p2[1]: 
              xinters=(y-p1[1])*(p1[0]-p2[0])/(p1[1]-p2[1])+p1[0]
              if p1[0] == p2[0] or x <= xinters: 
                counter=counter+1
      p1=p2
    return (counter+1) % 2 

