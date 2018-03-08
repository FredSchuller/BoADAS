#!/usr/bin/env python

# APEX - Atacama Pathfinder EXperiment Project

# artemis-util.py : collection of utilities and helper functions used
#                   for the reduction of Artemis data
#
# who         when        what
# --------    ----------  ----------------------------------------------
# fschuller   2016-11-27  Created: moved or redefined some functions
#                         from fred.py to here
# fschuller   2017-02-10  Function getPWVtoTau defined
#
"""
Name: artemis-util.py

This module defines a number of utilities and helper functions,
useful for the reduction of ArTeMiS data.

Content:
  - getArtemisRCP(mjdref,band)
  - getArtemisCalFactor(mjdref,band)
  - getPWVtoTau(mjdref,band)
  - shift(data,N)
  - applyFlatField(data,chanList)
  - center(radius)
  - scanPWV(data,rms)
  - plots2n(data)
  - computeNEFD(data)
  - get_model_tau(frontend, pwv)

Version: 2017-06-28

"""


# Here's the location where ancillary files should be located:
ARTEMIS_DIR = os.getenv('BOA_HOME_LABOCA')+'/../artemis/'

#------------------------------------------------------------------
#
# Helper functions to select proper RCP and calib. factor depending
# on observation date
#
#------------------------------------------------------------------
def getArtemisRCP(mjdref,band=350):
    """
    getArtemisRCP(mjdref,band=350):
      Returns the RCP file name for a given MJD date.
    INP:
      mjdref (float) : the MJD date of the observations
      band (int)     : select 350 or 450 microns
    """

    # Read the list of RCP files with MJD time intervals
    filename = ARTEMIS_DIR+'Artemis%i-rcps.dat'%(band)
    try:
        f = file(filename)
    except IOError:
        data.MessHand.error("Could not open %s."%(filename))
        return 0

    param = f.readlines()	
    f.close()

    # Look for the file corresponding to the input MJD
    rcpname = ''
    for i in range(len(param)):
        if param[i][0] not in ['!','#','n']:    # skip comments
            tmp = string.split(param[i])
            mjdstart = float(tmp[0])
            mjdstop  = float(tmp[1])
            if mjdstart <= mjdref and mjdref < mjdstop:
                rcpname = tmp[2]

    return rcpname

def getArtemisCalFactor(mjdref,band=350):
    """
    getArtemisCalFactor(mjdref,band=350):
      Returns the calibration factor (V-to-Jy) for a given MJD date.
    INP:
      mjdref (float) : the MJD date of the observations
      band (int)     : select 350 or 450 microns
    """

    # Read the list of cal. factors with MJD time intervals
    filename = ARTEMIS_DIR+'Artemis%i-calfactor.dat'%(band)
    try:
        f = file(filename)
    except IOError:
        data.MessHand.error("Could not open %s."%(filename))
        return 0

    param = f.readlines()	
    f.close()

    # Look for the line corresponding to the input MJD
    calfactor = 0.
    for i in range(len(param)):
        if param[i][0] not in ['!','#','n']:    # skip comments
            tmp = string.split(param[i])
            mjdstart = float(tmp[0])
            mjdstop  = float(tmp[1])
            if mjdstart <= mjdref and mjdref < mjdstop:
                calfactor = float(tmp[2])

    return calfactor

def getPWVtoTau(mjdref,band=350):
    """
    getPWVtoTau(mjdref,band=350):
      Returns the a,b coefficients of a linear approximation to convert
      PWV to zenith opacity, for a given MJD date.
    INP:
      mjdref (float) : the MJD date of the observations
      band (int)     : select 350 or 450 microns
    """
    # Read the list of coefficients with MJD time intervals
    filename = ARTEMIS_DIR+'Artemis%i-PWVtoTau.dat'%(band)
    try:
        f = file(filename)
    except IOError:
        data.MessHand.error("Could not open %s."%(filename))
        return 0,0

    param = f.readlines()	
    f.close()

    # Look for the line corresponding to the input MJD
    a,b = 0.,0.
    for i in range(len(param)):
        if param[i][0] not in ['!','#','n']:    # skip comments
            tmp = string.split(param[i])
            mjdstart = float(tmp[0])
            mjdstop  = float(tmp[1])
            if mjdstart <= mjdref and mjdref < mjdstop:
                a,b = float(tmp[2]),float(tmp[3])

    return a,b


#------------------------------------------------------------------
#
# Functions to perform basic tasks on the data
#
# -----------------------------------------------------------------
def shift(data,N=18):
    # Shift the data by N frames
    bak = copy.deepcopy(data.Data)
    dT  = fStat.f_median(data.ScanParam.get('deltaT'))
    if N > 0:
        data.Data[N:,:] = bak[:-N,:]
        data.flagMJD(below=N*dT+dT/2.)
    else:
        data.Data[:N,:] = bak[-N:,:]

# -----------------------------------------------------------------
def applyFlatField(data,chanList=[]):
    """
    Apply flat field by dividing signals by bolo gain, contained
    in data.BolometerArray.Gain
    This is a bug-fixed version of a buggy Boa method.
    """

    # get channel list
    chanList = data.BolometerArray.checkChanList(chanList)

    if len(chanList)<1: 
        data.MessHand.error("no valid channel")
        return

    for chan in chanList:
        num = data.BolometerArray.getChanIndex(chan)[0]
        numGain = chan-1
        data.Data[:,num] = data.Data[:,num] / array((data.BolometerArray.Gain[numGain]),'f')
        
    data._DataAna__resetStatistics()

#------------------------------------------------------------------
#
# Other utility functions
#
#------------------------------------------------------------------
def center(radius=40.):
    """
    center(radius):
      Return a list of pixels within a given radius of (0,0)
    """
    cc = data.BolometerArray.checkChanList([])
    result = []
    for c in cc:
        offX = data.BolometerArray.Offsets[0,c-1] - data.BolometerArray.RefOffX
        offY = data.BolometerArray.Offsets[1,c-1] - data.BolometerArray.RefOffY
        if sqrt(offX**2 + offY**2) < radius:
            result.append(c)
    return result

# -----------------------------------------------------------------
def scanPWV(data,rms=0):
    """
    scanPWV(data):
      Returns mean value of PWV in the current scan.
      Data must have been read with argument readPWV=1
    """
    pwv = data.ScanParam.PWV
    if not pwv:
        print " ***  No PWV available - returning 0 !!  ***"
        return 0.

    # Convert this list of list of arrays to a 1D array
    # pp = ravel(pwv)
    # but this fails if there is more than 1 value in a subscan - changed 2017-06-28
    result = []
    for p in pwv:
        result.extend(p)
    pp = ravel(result)
    
    # Compute mean value
    valmean = fStat.f_mean(pp)
    if rms:
        rmsPWV  = fStat.f_rms(pp,valmean)
        return valmean,rmsPWV
    else:
        return valmean

# -----------------------------------------------------------------
def plots2n(data):
    """
    plots2n(data):
      Build a map showing s2n (or 1/rms) per pixel
    """
    local_data = copy.deepcopy(data)
    if not local_data._DataAna__statisticsDone:
        local_data._DataAna__statistics()
    nbolo = local_data.BolometerArray.NUsedChannels
    local_data.ScanParam.AzOff = array([0.,0.],'f')
    local_data.ScanParam.ElOff = array([0.,0.],'f')
    local_data.Data = ones((2,nbolo),'f')
    local_data.DataWeights = ones((2,nbolo),'f')
    for i in range(nbolo):
        local_data.Data[0,i] = 1./local_data.ChanRms[i]
        local_data.Data[1,i] = 1./local_data.ChanRms[i]
    
    local_data.FlagHandler._aFlags = zeros((2,nbolo),'i')
    local_data.doMap(oversamp=2,showRms=0)

# -----------------------------------------------------------------
def computeNEFD(data):
    """
    Computes, print and return the NEFD.
    """
    chans = data.BolometerArray.checkChanList([])
    if not data._DataAna__statisticsDone:
        data._DataAna__statistics()
    rms = []
    for c in chans:
        rms.append(data.ChanRms[data.BolometerArray.getChanIndex(c)[0]])
    rms = array(rms,'f')    
    sampling = 1./fStat.f_median(data.ScanParam.get('deltaT'))
    # Weighted average NEFD
    w_nefd = sum(1./rms) / sum(1./rms**2) / sqrt(sampling)

    data.MessHand.info("--------------------------------------------")
    data.MessHand.info("Weighted average NEFD = %6.1f mJy sqrt(Hz)"%(w_nefd * 1.e3))
    data.MessHand.info("--------------------------------------------")
    return w_nefd

# -----------------------------------------------------------------
def get_model_tau(frontend='saboca', pwv=0.5):
    '''
    This function gives the value of tau_z calculated by the ATM model
    for the specified bolometer camera and given pwv
    INP:    (sting) frontend : it can be 'laboca' or 'saboca', must be lowercase
                 (float) pwv : the radiometer los water vapor in mm
    OUT:         (float) tau : model output
    '''
    outfile = os.path.join(os.getenv('HOME'),'watertau.dat')
    if frontend == 'laboca':
        cmd = 'curl -s www.apex-telescope.org/dbtest/atm/tau-water.php?water='+str(pwv)+'  > ' + outfile

    if frontend == 'saboca':
        cmd = 'curl -s www.apex-telescope.org/dbtest/atm/tau-water-saboca.php?water='+str(pwv)+'  > ' + outfile

    if (frontend <> 'saboca' and frontend <> 'laboca'):
       print 'Wrong frontend definition!'
       return

    os.system(cmd)
    f = file(os.path.join(os.getenv('HOME'),'watertau.dat'))
    s = f.read()
    f.close()
    t = s.split()
    tau = t[3]
    try:
        tau = float(tau[0:9])
    except ValueError:
        data.MessHand.warning("*************************************")
        data.MessHand.warning("Cannot use tau: %s - using tau = 0.3"%(str(tau)))
        data.MessHand.warning("*************************************")
        tau = 0.3
        
    return tau

# -----------------------------------------------------------------
def pwv_to_tau(data,band=350):
    '''
    Convert PWV to tau using a linear approximation.
    '''
    mjdref = data.ScanParam.MJD[0]
    a,b = getPWVtoTau(mjdref,band=band)
    if a == 0 and b == 0:
        print "No linear fit coefficients found for this date. Exiting."
        return 0.

    pwv = scanPWV(data)
    if pwv:
        return a * pwv + b
    else:
        return 0.
    
# -----------------------------------------------------------------
# -----------------------------------------------------------------
def face():
    """
    Perfectly useless function, just for fun.
    """
    print "  ^-^-^\n /_   _\ \n( O , O )\n |  ^  |\n \ <=> /\n  \_ _/\n    V"

def show():
    """
    Displays some basic information about the current data.
    """
    print "%i  %-9s %4.1f  %4.2f  %4.2f   %4.2f\n"%(data.ScanParam.ScanNum,
                                                    data.ScanParam.Object,
                                                    fStat.f_median(data.ScanParam.get('el')),
                                                    scanPWV(data),
                                                    pwv_to_tau(data,band=350),
                                                    pwv_to_tau(data,band=450))
