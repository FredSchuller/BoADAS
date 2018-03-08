##########################################################
#
# boaSaboca.py
#
# Boa script for pipeline reduction of SABOCA data
#
# Created: 2009/xx/xx - copy of LABOCA functions
#
#
# who       when        what
# -----------------------------------------------------
# mdumke    2010-05-01  Added invertSomeChannels 
#                           & getSabocaInvert functions
# aweiss    2010-05-15  Added calcsensitivity
# fschuller 2010-06-07  Use BOA_HOME_SABOCA
##########################################################
from boa.fortran import fStat

# Read SABOCA specific definitions
if not os.getenv('BOA_HOME_SABOCA'):
   raise 'Environment variable BOA_HOME_SABOCA undefined'
sabocadir = os.getenv('BOA_HOME_SABOCA') + '/'
execfile(sabocadir+'saboca-cabling.py')

# Function to retrieve the appropriate RCP file as a function of MJD
def getSabocaRCP(mjdref):
    try:
        f = file(os.path.join(os.getenv('BOA_HOME_LABOCA'), '../saboca/Saboca-rcps.dat'))
    except IOError:
        self.__MessHand.error("could not open Saboca-rcps.dat")	
        return

    param = f.readlines()	
    f.close()

    mjdstart = []
    mjdstop = []
    rcp=[]
    for i in range(len(param)):
        if param[i][0] not in ['!','#','n']:    # skip comments
	    tmp = string.split(param[i])
	    mjdstart.append(string.atof(tmp[0]))
	    mjdstop.append(string.atof(tmp[1]))
	    rcp.append(tmp[2])

    mjdstart = array(mjdstart)
    mjdstop =  array(mjdstop)
    rcpname= -1
	
    for i in range(len(mjdstart)):
        low = mjdstart[i]
	high = mjdstop[i]
	if low <= mjdref and mjdref < high:
	    rcpname = rcp[i]

    return rcpname


# Function to get list of cross-talkers as a function of MJD
def getSabocaCross(mjdref):	
    try:
        f = file(os.path.join(os.getenv('BOA_HOME_LABOCA'), '../saboca/Saboca-cross.dat'))
    except IOError:
        self.__MessHand.error("could not open Saboca-cross.dat")	
        return

    param = f.readlines()	
    f.close()

    mjdstart = []
    mjdstop = []
    cr=[]
    for i in range(len(param)):
        if param[i][0] not in ['!','#','n']:    # skip comments
	    tmp = string.split(param[i])
	    mjdstart.append(string.atof(tmp[0]))
	    mjdstop.append(string.atof(tmp[1]))
	    cr.append(tmp[2])

    mjdstart = array(mjdstart)
    mjdstop =  array(mjdstop)
    	
    cross= []
    crossreturn = []
    for i in range(len(mjdstart)):
        low = mjdstart[i]
	high = mjdstop[i]
	if low <= mjdref and mjdref < high:
	    cross = cr[i]

    tmp=cross.split(',')
    for a in tmp:
        crossreturn.append(int(a))

    return crossreturn

# Function to invert some channels
def invertSomeChannels(iclist=[]):
    for i in range(len(iclist)):
        num=iclist[i]-1
        data.Data[:,num] *= array(-1,'f')

# Function to get list of channels to be inverted as a function of MJD
def getSabocaInvert(mjdref):
    try:
        f = file(os.path.join(os.getenv('BOA_HOME_LABOCA'), '../saboca/Saboca-invert.dat'))
    except IOError:
        self.__MessHand.error("could not open Saboca-invert.dat")	
        return

    param = f.readlines()	
    f.close()

    mjdstart = []
    mjdstop = []
    inv=[]
    for i in range(len(param)):
        if param[i][0] not in ['!','#','n']:    # skip comments
	    tmp = string.split(param[i])
	    mjdstart.append(string.atof(tmp[0]))
	    mjdstop.append(string.atof(tmp[1]))
	    inv.append(tmp[2])

    mjdstart = array(mjdstart)
    mjdstop =  array(mjdstop)
    	
    invert = []
    invertreturn = []
    for i in range(len(mjdstart)):
        low = mjdstart[i]
	high = mjdstop[i]
	if low <= mjdref and mjdref < high:
	    invert = inv[i]

    tmp=invert.split(',')
    for a in tmp:
        invertreturn.append(int(a))
    if invertreturn == [0]:
        invertreturn = []

    return invertreturn

# Function to compute mean array sensitivity
def calcsensitivity():
    usechan = data.BolometerArray.checkChanList([])
    nrbolos = len(usechan)
    weight()
    rms = array((data.getChanListData('rms')),'f')
    adt=(data.ScanParam.get('deltat'))
    dt=fStat.f_median(adt)
    sensitivity=rms*sqrt(dt)*1000.
    array_sensitivity = sum(sensitivity/(rms**2))/sum(1./(rms**2))
    return  array_sensitivity
