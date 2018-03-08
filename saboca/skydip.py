#import time, sys, posix, cPickle

#from Numeric import *
from LinearAlgebra import *


try:
    import mpfit
except ImportError:
    from boa import mpfit

from boa.fortran  import fUtilities,fStat
from boa.BoaError import BoaError
from boa          import BoaConfig


def safeExp(x):
    """
    NAM: safeExp (function)
    DES: correct a bug in Numeric that raise an expection when
         computing exponential of small numbers, this take a lot of time !
         but faster thant converting to nummarray compute the exp and back to numeric !!
    """
    condition = less(x,-745.)
    tmp = choose(condition,(x,1.))
    return choose(condition,(exp(tmp),0.))


def skydipn(p,fjac=None,x=None,y=None,err=None):
    """
    NAM: skydip
    DES: signal as a function of elevation, aka skydip
    """
    model=modelSkydipn(p,x)
    status = 0
    return([status, (y-model)/err])


def modelSkydipn(p,x):
    """
    DES: model function for fitting skydip
         full model, with 3 parameters:
         t(x)=tatm*((1-exp(-tauz/sin(pi*x/180.0)))*feff+(1-feff))
         p[0]=tatm, p[1]=tauz, p[2]=feff
    """
    tmp1 = -1.*p[1]/sin(x*pi/180.)
    return p[0]*(1.-safeExp(tmp1))+p[2]


def fitSkydipn(x,y,err,val0):
    """
    DES: fits a skydip signal-elevation function
         3 parameters fitted: opt, tauz, feff
    INP: (f array) x = x data
         (f array) y = y data
         (f array) err=errors on y values
         (3xf list) val0 = first guess values, in this order:
                           [Tatm, tau_z, F_eff]
    """
    parname = ['Tatm','tau_z','F_eff']
    p = val0
    parinfo = []
    for i in range(3):
        parinfo.extend([{'parname': parname[i], \
                         'value': p[i], \
                         'fixed': 0, \
                         'limits' : [0.,0.],\
                         'limited': [0,0]}])
    #parinfo[0]['limits'] = [100.,520.]
    parinfo[0]['limited'] = [0,0]
    #parinfo[1]['limits'] = [0.,10.]  # don't observe at tau>3 !!
    #parinfo[2]['limits'] = [0.2,1.0]
    
    fa={'x':x, 'y':y, 'err':err}
    m=mpfit.mpfit(skydipn,p,parinfo=parinfo,functkw=fa,quiet=1,debug=0)
    if (m.status <= 0):
        raise BoaError, str("mpfit failed: %s"%(m.errmsg))
    return m


def correctBlind(data):
   # get full timestream, also at flagged timestamps,
   # to be able to modify the data.Data array
   unflagC([38,39])
   
   blind1 = copy.copy(data.getChanData('flux',38,flag='None'))
   blind2 = copy.copy(data.getChanData('flux',39,flag='None'))
   ok = data.BolometerArray.checkChanList([])
   nb = len(ok)
   for i in range(nb):
      num = data.BolometerArray.getChanIndex(ok[i])[0]
      correc1 = 1.0 * (blind1-blind1[0])
      correc2 = 1.0 * (blind2-blind2[0])
      correction = (correc1 + correc2)/2.  # average of both blind bolos
      correction = correction.astype('f')
      data.Data[:,num] -= correction
