# This file defines which channels are connected to resistors,
# short or unconnected, and defines functions specific to Laboca
#
# Last update: 2007/8/17 - F. Schuller

import os

if not os.getenv('BOA_HOME_RCP'):
   raise 'Environment variable BOA_HOME_RCP undefined'
filename = os.path.join(os.getenv('BOA_HOME_RCP'), 'channels_full_set.20070509.txt')

# File format:
#BEch   box     group   amp     board   cable   pin     pixel   state
#1       1       1       3       A       I       11      206     B   

# Read "state" column
#data.BolometerArray.readAdditionnalIndexFile(filename,indexColumn=8,comment='#')
#resistor = data.BolometerArray.selectAdditionnalIndex('R')
resistor = [ 15,  16,  87,  88,  95, 139, 143, 145, 149, 150, 157, 178,
             181, 186, 187, 188, 190, 255, 310, 317]

# Some peculiar channels, as determined in May 2007
# cross-talk, as of 2007 May 10
#cross = [87,101,103,106,110,139,143,149,163,183,241,255] 
# update 2007 May 20 - the others are resistors, except 183
cross = [101,103,106,110,241]
blind = [4,163]
sealed_may07= blind
bad_may07=[86,90,96,136,173,183,191,192,201,202,204,205,218,220,309]

# Conversion factor for LABOCA - May 2007
VtoJy = 6.3E6

# The following function does correlated noise removal by some group of channels
def correlSelection(data,filename,refColumn,indexColumn,comment,listOfGroup,
                    factor=0.8,nbloop=1):
    """
    DES: Remove skynoise by grouping the bolometer with a certain index read in
         a file. Relative gains are computed with respect to median of signals.
    """
    data.BolometerArray.readAdditionnalIndexFile(filename,
                                                 refColumn,indexColumn,comment)
    for group in listOfGroup:
        chanList = data.BolometerArray.selectAdditionnalIndex(group)
        #data.medianNoiseRemoval(chanList,chanRef=-2,
        #                        factor=factor,nbloop=nbloop)
        # changed 2009/1/25: use mean instead of median for reference
        # AND other change: requires at least 6 channels, otherwise
        # don't subtract anything
        if len(chanList) > 5:
           data.medianNoiseRemoval(chanList,chanRef=-1,
                                   factor=factor,nbloop=nbloop)        
        else:
           data.MessHand.warning("Only %i channels in group %s"%(len(chanList),
                                                                 group))
           data.MessHand.warning("... no correlated noise subtracted")
           
# One group contains 26 channels, connected to the same cable
def correlgroup(data,nbloop=2,factor=0.9):
    """
    subtract correlated noise per cable
    """
    correlSelection(data,filename,0,5,'#',
                    ['A','B','C','D','E','F','G','H','I','K','L','M'],
                    factor=factor,nbloop=nbloop)

# One box contains 80 channels in total
def correlbox(data,nbloop=2,factor=0.9):
    """
    subtract correlated noise per cable
    """
    correlSelection(data,filename,0,1,'#',
                    [1,2,3,4],
                    factor=factor,nbloop=nbloop)

# Conversion raw data to Volts
def CntstoV(data):
   be = data.BolometerArray.BEGain
   data.Data *= array(be/270.,'f')
   data._DataAna__resetStatistics()
   
def CntstoJy(data):
   CntstoV(data)
   data.Data *= array(VtoJy,'f')
   data._DataAna__resetStatistics()
   
############################################################
# Function to correct for He3 temperature drifts
# first, read the channel numbers-gains dictionary
he3dict = restoreFile(os.getenv('BOA_HOME_LABOCA')+'/he3gains.dat')

def correctHe3(data):
   # get full timestream, also at flagged timestamps,
   # to be able to modify the data.Data array
   he = data.ScanParam.he3SmoothInterpolate(flag='None')
   he = 1000.*he        # convert to mK
   he = he.astype('f')  # and to single precision floats
   ok = data.BolometerArray.checkChanList([])
   nb = len(ok)
   for i in range(nb):
      num = data.BolometerArray.getChanIndex(ok[i])[0]
      correction = he3dict[str(ok[i])] * (he-he[0]) * 1e-6
      correction = correction.astype('f')
      data.Data[:,num] -= correction

############################################################
# Function to correct for He3 T drifts using blind bolos
# first, read the channel numbers-gains dictionary
blind_gains = restoreFile(os.getenv('BOA_HOME_LABOCA')+'/blind-gains.dat')

def correctBlind(data):
   # get full timestream, also at flagged timestamps,
   # to be able to modify the data.Data array
   blind1 = copy.copy(data.getChanData('flux',4,flag='None'))
   blind2 = copy.copy(data.getChanData('flux',163,flag='None'))
   ok = data.BolometerArray.checkChanList([])
   nb = len(ok)
   for i in range(nb):
      num = data.BolometerArray.getChanIndex(ok[i])[0]
      correc1 = blind_gains[str(ok[i])][0] * (blind1-blind1[0])
      correc2 = blind_gains[str(ok[i])][1] * (blind2-blind2[0])
      correction = (correc1 + correc2)/2.  # average of both blind bolos
      correction = correction.astype('f')
      data.Data[:,num] -= correction
