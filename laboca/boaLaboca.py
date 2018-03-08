##########################################################
#
# boaLaboca.py
#
# Boa script for pipeline reduction of LABOCA data
#
# Created: 2007/9/08 (copy of older scripts)
# Last modif: 2010/6/02
#
# Author: F. Schuller, A. Weiss
#
##########################################################

# Read LABOCA specific definitions
execfile(os.getenv('BOA_HOME_LABOCA')+'/cabling.py')

from boa.fortran import fStat
from boa.Utilities import Timing
myTiming = Timing()

def reduceLaboca(data,
                 model=None, flag0=5., flag1=5,
                 corBlind=1, corHe3=0,
                 fac1=0.95, nbl1=2, fac2=0.94, nbl2=3,
                 flag2=5, spike1=10,
                 fac3=0.94, nbl3=3, fac4=0.90, nbl4=2,
                 flag3=5, src=4,
                 below=0.1, hiref=0.25, order=1,
                 flag1p5=0):
    """
    DES: generic function to process Laboca data
    INP: Mandatory: 'data' = the data object to work on
         Optional: many arguments with default values.

         IMPORTANT NOTE: this function does NOT apply any calibration,
         nor it updates the RCP. The reason for this is: this should
         be done before calling this function, so that a calibrated
         model map can eventually be subtracted from the data before
         running the processing.
         
         The succession of steps performed by this function is as follows:
           if model: flagSource(model,threshold=flag0)
           flagFractionRms(ratio=flag1)
           if corBlind: correctBlind(data)
           if corHe3: correctHe3(data)
           medianNoiseRemoval(chanRef=-2,factor=fac1,nbloop=nbl1)
           flagFractionRms(ratio=flag1p5)
           medianNoiseRemoval(chanRef=-1,factor=fac2,nbloop=nbl2)
           flagFractionRms(ratio=flag2)
           despike(above=spike1,below=-1*spike1)
           correlbox(data,factor=fac3,nbloop=nbl3)
           correlgroup(data,factor=fac4,nbloop=nbl4)
           flagFractionRms(ratio=flag3)
           despike(above=src,flag=8)
           data.flattenFreq(below=below,hiref=hiref)
           base(order=order)
           
    NOTE: this function should be called AFTER the data have been read
    """
    global myTiming
    myTiming.setTime()

    if model:
        data.flagSource(threshold=flag0,model=model)  # source model

    # Flag bad channels
    cross = getLabocaCross(data.ScanParam.MJD[0])
    data.flagChannels(cross)
    data.flagChannels(resistor)
    data.flagChannels(sealed_may07)

    if corBlind:
        # Correct for "blind" signal
        correctBlind(data)
    elif corHe3:
        # Correct for monitor readings of He3 T
        correctHe3(data)
        
    data.zeroEnds()
    data.flagFractionRms(ratio=flag1)

    # Flag stationary points and high acceleration
    data.flagSpeed (below=20.)
    data.flagSpeed(above=500.)
    data.flagAccel(above=800.)

    data.medianBaseline()

    # First correlated noise removal on all channels
    if fac1:
        data.medianNoiseRemoval(chanRef=-2,factor=fac1,nbloop=nbl1)
    if flag1p5:
        data.flagFractionRms(ratio=flag1p5,below=0)

    if fac2:
        data.medianNoiseRemoval(chanRef=-1,factor=fac2,nbloop=nbl2)
    data.flagFractionRms(ratio=flag2,below=0)

    data.despike(above=spike1,below=-1*spike1)

    # default: fac3=0.94, nbl3=3
    if fac3:
        # correlated noise removal by amplifier box
        correlbox(data,factor=fac3,nbloop=nbl3)
        data.flagFractionRms(ratio=flag3,below=0)

    # default: fac4=0.90, nbl4=2
    if fac4:
        # correlated noise removal by groups of channels
        # connected to the same cables
        correlgroup(data,factor=fac4,nbloop=nbl4)

    # Flag noisy channels
    data.flagFractionRms(ratio=flag3,below=0)

    if src:
        data.despike(above=src,flag=8)  # probably sources

    # flag NaNs - this should not hurt...
    data.flagNan()
    
    # Low frequency filtering
    if below:
        data.flattenFreq(below=below,hiref=hiref)

    # baseline on full scan
    data.polynomialBaseline(order=order)

    data.computeWeight()
    data.unflag(flag=8)  # sources back

    print myTiming

def jumps(ratio = 2.):
    """
    DES: This function flags time intervals corresponding to spikes (or jumps)
         seen by all bolometers, as inferred form the variations in the median
         of absolute signals of all channels. Timestamps where this quantity
         is above ratio x its median value are flagged.
    INP: (float) ratio = the threshold above which to consider a spike
    OUT: a list of time intervals is returned. Each element is a 2-element
         tuple: (time_start, time_end)
    NOTE: this is not specific to Laboca data
    """
    global data
    
    m = data._DataAna__computeMedianAbsSignal()
    t = data.ScanParam.get('mjd')
    # smooth them
    nb = len(t)  # number of timestamps
    nb2 = int(nb/4)
    m2 = zeros((nb2),'f')
    t2 = zeros((nb2),'f')
    for i in range(nb2):
        m2[i] = (m[4*i] + m[4*i+1] + m[4*i+2] + m[4*i+3])/4.
        t2[i] = (t[4*i] + t[4*i+1] + t[4*i+2] + t[4*i+3])/4.

    mm = fStat.f_mean(m2)
    toFlag = []
    i = 0

    start = -1
    while i < nb2-1:
        while m2[i] < mm * ratio and i < nb2-1:
            i += 1
        if i<nb2-1:
            # start of a jump to flag
            start = max([0,i-1])
            while m2[i] >= mm* ratio and i < nb2-1:
                i += 1
            # here it goes back below theshold
            end = min([i+1,nb2-1])
        
            if start >= 0:
                toFlag.append([t2[start],t2[end]])

    print "found %i jumps"%(len(toFlag))
    if toFlag:
        # if anything found
        for gap in toFlag:
            data.flagMJD(above=gap[0],below=gap[1])

    return toFlag

# -----------------------------------------------------------------------------
# Functions to get the correct RCP file and list of bad channels at a given MJD
#
def getLabocaRCP(mjdref):
    try:
        f = file(os.path.join(os.getenv('BOA_HOME_LABOCA'), 'Laboca-rcps.dat'))
    except IOError:
        self.__MessHand.error("could not open Laboca-rcps.dat")	
        return

    param = f.readlines()	
    f.close()

    mjdstart = []
    mjdstop = []
    rcp=[]
    for i in range(len(param)):
        if param[i][0] not in ['!','#','n']:    # skip comments
            tmp = string.split(param[i])
            mjdstart.append(float(tmp[0]))
            mjdstop.append(float(tmp[1]))
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

def getLabocaCross(mjdref):
    try:
        f = file(os.path.join(os.getenv('BOA_HOME_LABOCA'), 'Laboca-cross.dat'))
    except IOError:
        self.__MessHand.error("could not open Laboca-cross.dat")	
        return

    param = f.readlines()	
    f.close()

    mjdstart = []
    mjdstop = []
    cr=[]
    for i in range(len(param)):
        if param[i][0] not in ['!','#','n']:    # skip comments
            tmp = string.split(param[i])
            mjdstart.append(float(tmp[0]))
            mjdstop.append(float(tmp[1]))
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


