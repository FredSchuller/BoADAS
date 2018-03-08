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

def skydipPlus(p,fjac=None,x=None,y=None,err=None):
    """
    NAM: skydipPlus
    DES: function called by mpfit to fit skydip
    """
    model = modelSkydipPlus(p,x)
    status = 0
    return([status, (y-model)/err])

def skydipFull(p,fjac=None,x=None,y=None,err=None):
    """
    NAM: skydipFull
    DES: function called by mpfit to fit skydip
    """
    model = modelSkydipFull(p,x)
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
    # return p[0]*((1.-safeExp(tmp1))*p[2] + (1.-p[2]))
    # changed 2017-07-03: the second term should be (1-eta)*T_amb and not (1-eta)*T_atm
    T_amb = data.ScanParam.T_amb
    return p[0]*(1.-safeExp(tmp1))*p[2] + (1.-p[2])*T_amb


def modelSkydipPlus(p,x):
    """
    DES: Full model function for fitting skydip:
         T(x) = Feff * Tatm * (1 - exp(-tau/sin(x))) + (1 - Feff) * Tamb
         and:
         S(x) = T(x) / Jy2K + offset

         Parameters to be fitted:
         p[0] = Tatm
         p[1] = tau
         p[2] = Feff
         p[3] = Jy2K
         p[4] = offset
    """
    tmp1 = -1.*p[1]/sin(x*pi/180.)
    T_amb = data.ScanParam.T_amb
    tmp2 = p[0] * p[2] * (1.-safeExp(tmp1)) + (1.-p[2]) * T_amb
    
    return tmp2 / p[3] + p[4]

def modelSkydipFull(p,x):
    """
    DES: Even fuller model function for fitting skydip:
         T(x) = C * (Feff * Tatm * (1 - exp(-tau/sin(x))) + (1 - Feff) * Tamb) + (1-C) * T_cabin
         and:
         S(x) = T(x) / Jy2K + offset

         Parameters to be fitted:
         p[0] = Tatm
         p[1] = tau
         p[2] = Feff
         p[3] = Jy2K
         p[4] = offset
         p[5] = C
    """
    tmp1 = -1.*p[1]/sin(x*pi/180.)
    T_amb = data.ScanParam.T_amb
    tmp2 = p[0] * p[2] * (1.-safeExp(tmp1)) + (1.-p[2]) * T_amb
    # add coupling term
    T_cabin = 283.  # 10 deg C
    tmp3 = p[5] * tmp2 + (1.-p[5]) * T_cabin
    return tmp3 / p[3] + p[4]


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
    #parinfo[0]['limited'] = [0,0]
    #parinfo[1]['limits'] = [0.,10.]  # don't observe at tau>3 !!
    #parinfo[2]['limits'] = [0.2,1.0]
    
    fa={'x':x, 'y':y, 'err':err}
    m=mpfit.mpfit(skydipn,p,parinfo=parinfo,functkw=fa,quiet=1,debug=0)
    if (m.status <= 0):
        raise BoaError, str("mpfit failed: %s"%(m.errmsg))
    return m

def fitSkydipPlus(x,y,err,val0):
    """
    DES: fits a skydip signal-elevation function
         5 parameters fitted:
         p[0] = Tatm
         p[1] = tau
         p[2] = Feff
         p[3] = Jy2K
         p[4] = offset
    """
    parname = ['Tatm','tau_z','F_eff','Jy2K','offset']
    p = val0
    parinfo = []
    for i in range(5):
        parinfo.extend([{'parname': parname[i], \
                         'value': p[i], \
                         'fixed': 0, \
                         'limits' : [0.,0.],\
                         'limited': [0,0]}])
    
    fa={'x':x, 'y':y, 'err':err}
    m=mpfit.mpfit(skydipPlus,p,parinfo=parinfo,functkw=fa,quiet=1,debug=0)
    if (m.status <= 0):
        raise BoaError, str("mpfit failed: %s"%(m.errmsg))
    return m


def fitSkydipFull(x,y,err,val0):
    """
    DES: fits a skydip signal-elevation function
         6 parameters fitted:
         p[0] = Tatm
         p[1] = tau
         p[2] = Feff
         p[3] = Jy2K
         p[4] = offset
         p[5] = coupl
    """
    parname = ['Tatm','tau_z','F_eff','Jy2K','offset','coupl']
    p = val0
    parinfo = []
    for i in range(6):
        parinfo.extend([{'parname': parname[i], \
                         'value': p[i], \
                         'fixed': 0, \
                         'limits' : [0.,0.],\
                         'limited': [0,0]}])
    
    fa={'x':x, 'y':y, 'err':err}
    m=mpfit.mpfit(skydipFull,p,parinfo=parinfo,functkw=fa,quiet=1,debug=0)
    if (m.status <= 0):
        raise BoaError, str("mpfit failed: %s"%(m.errmsg))
    return m




def fnSkydip(x,Tatm,Feff,tau):
    Tspill = 283.
    tmp = -1.*tau / sin(x*pi/180.)
    result = Tatm * Feff * (1.-exp(tmp)) + (1.-Feff) * Tspill
    # result = Tatm * Feff * (1.-exp(tmp)) + offset
    return result
