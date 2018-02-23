# Copyright (C) 2002-2007
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
NAM: BoaDataEntity.py (file)
DES: contains the BoA data entity class
"""                    


__version__=  '$Revision: 2811 $'
__date__=     '$Date: 2015-03-04 15:32:01 +0100 (Wed, 04 Mar 2015) $'
# __tag__=      '$Name:  $'

#----------------------------------------------------------------------------------
#----- Import ---------------------------------------------------------------------
#----------------------------------------------------------------------------------

import os, string, time, cPickle, gc, copy

from Numeric import *
from slalib import *   
import arrayfns

from boa           import BoaMBFits, BoaMBFitsReader, BoaFlagHandler
from boa           import BoaMessageHandler, BoaCommandHistory, BoaConfig, BoaDir
from boa.Bogli     import Plot, MultiPlot, Forms
from boa.fortran   import fUtilities,fStat
from boa.Utilities import gaussian, Timing, attrStr, compressNan, prettyPrintList
from boa.Utilities import fitGaussian, modelgauss,tolist_boa

#----------------------------------------------------------------------------------
#----- BolometerArray Class -------------------------------------------------------
#----------------------------------------------------------------------------------

class Telescope:
    """
    NAM: Telescope (class)
    DES: Define all the useful parameters of a telescope
    """

    def __init__(self):
        """
        DES: Instanciation of a Telescope object
        """

        self.Name = ""                    # Telescope name
        self.Diameter  = 0.0              # in m
        self.Longitude = 0.0              # telescope longitude (in deg)
        self.Latitude  = 0.0              # telescope latitude  (in deg)
        self.Elevation = 0.0              # telescope altitude  (in m)

    def set(self,name="",diameter=0.0,longitude=0.0,latitude=0.0,elevation=0.0):
        """
        DES: set all the parameters
        """

        self.Name      = name
        self.Diameter  = diameter
        self.Longitude = longitude
        self.Latitude  = latitude
        self.Elevation = elevation

    def __str__(self):
        """
        DES: Defines a string which is shown when the print
             instruction is used.
        """
        return self.Name + " (" + "%3i m)" % self.Diameter

#----------------------------------------------------------------------------------
#----- BolometerArray Class -------------------------------------------------------
#----------------------------------------------------------------------------------

class BolometerArray:
    """
    NAM: BolometerArray (class)
    DES: Define all the useful parameters of a bolometer array
    """
    
    def __init__(self):
        """
        DES: Instanciation of a BolometerArray object
        """
        
        # Add a MessHand attribute - new MessageHandler 20050303
        self.__MessHand=BoaMessageHandler.MessHand(self.__module__)
        
        self.FeBe = ""                             # Backend Frontend combinaison
        self.Telescope = Telescope()
        
        self.TransmitionCurve = array([],Float32)  # 2D (frequency vs transmition)
        self.EffectiveFrequency  = 0.0             # Hz

        self.NChannels     = 0                     # The total number of channels
        self.NUsedChannels = 0                     # Number of channels in use, i.e. size of data
        self.RefChannel    = 0                     # The reference Channel
        self.UsedChannels  = array([],Int)         # List of used Channels (to map to data array)
        self.CurrChanList  = array([],Int)         # Current list of used array
        
        self.Offsets    = array([[],[]], Float32)  # The so called RCP in arcsec
        self.FWHM       = array([[],[]], Float32)  # Corresponding FWHM of modelled gaussian (major-minor) (arcsec)
        self.Tilt       = array([],      Float32)  # with tilt in degree
        self.AddIndex   = []                       # an additionnal index list

        
        self.Gain       = array([], Float32)       # Normalized gains (for point sources)
        self.ExtGain    = array([], Float32)       # Normalized gains (for extended emission, skynoise)
        self.JyPerCount = 1.                       # Jy / count conversion factor
        self.BEGain     = 1.                       # Backend attenuation factor (new v. 1.61)
        self.FEGain     = 0.                       # Frontend amplifier gain (actually the gain
                                                   # is 2^FEGain) - new v. 1.61
        self.DCOff      = array([], Float32)       # DC offsets

        self.DewUser    = 0.                       # Dewar angle (user relative to coord. system)
        self.DewExtra   = 0.                       # Extra dewar rotation angle
        
        # Derived parameters
        self.BeamSize = 0.                         # size of the beam in arcsec
        self.ChannelSep = array([[],[]],Float)     # the separation between channels

        self.FeedType    = []                      # array defining the type of the used feed
        self.FeedCode    = {}                      # Dictionnary describing the feed type
                                                   # FeedCode of 1 will be sky feeds

        self.FlagHandler = BoaFlagHandler.createFlagHandler(array([], Int32))
        self.cflags       = {'NOT CONNECTED'  : 1,
                             'BAD SENSITIVITY': 2,
                             'LOW SENSITIVITY': 3,
                             'DARK BOLOMETER' : 4,
                             'TEMPORARY'      : 8 }

    #----------------------------------------------------------------------------
    def __str__(self):
        """
        DES: Defines a string which is shown when the print
             instruction is used.
        """
        out = self.FeBe
        out += " - %3i/%3i/%3i channels at %3i GHz (%s)"% \
               (self.FlagHandler.nUnset(), \
                len(self.UsedChannels), \
                self.NChannels, \
                int(self.EffectiveFrequency/1.e9), \
                self.printCurrChanList())
        
        #out += " on " + self.Telescope.Name + " (" + \
        #       "%3i m" % self.Telescope.Diameter + \
        #       " - %4.1f\" default beam)" % self.BeamSize

        out += " - FE/BE Gain: %2i/%3i"%\
               (2.**self.FEGain,self.BEGain)
        
        if BoaConfig.DEBUG > 2:
            out += "\n" + \
                   attrStr(self,['Telescope','_BolometerArray__MessHand']) + \
                   "\n"

        return out

    # ---------------------------------------------------------------------
    def __fillFromMBFits(self,reader,febe,baseband,subscan,flag=1):
        """
        DES: fill a BolometerArray object using the MBFitsReader object reader.
        INP:            reader : MBFitsReader object
                    febe (str) : frontend-backend name to select
                baseband (int) : baseband number to select
             subscans (i list) : list of subscans numbers to read in
                 flag (i list) : flag for not connected feeds (default: 1 'NOT CONNECTED')
        """
        try:

            # MBFit files can have several UseBand take the first
            # corresponding to the baseband
            useBand = reader.read("UseBand", febe=febe)
            indexBaseband = -1
            for iBand in xrange(len(useBand)):
                if useBand[iBand] == baseband:
                    indexBaseband = iBand

            # TODO: Why some parameters needs indexBaseband and some other not
            nChannels    = reader.read("FebeFeed", febe=febe)
            self.NChannels = nChannels
            usedChannels = reader.read("UseFeed",  febe=febe)[indexBaseband]
            allUsedChannels = usedChannels
            nbUsedChan   = reader.read("NUseFeeds",febe=febe)[indexBaseband]
            usedChannels = usedChannels[:nbUsedChan]
            refChannel   = reader.read("RefFeed",  febe=febe)
            origRefChannel = refChannel

            # FEEDTYPE: tells us which channels are real bolometers, which are "dark"
            feedType   = reader.read("FeedType",febe=febe)[indexBaseband]
            feedString = reader.read("FeedCode",febe=febe)
            
            # Convert FeedCode to dictionnary
            feedCode = {}
            listType = feedString.split(',')
            goodType = []
            for i in range(len(listType)):
                num_type = listType[i].split(':')
                if len(num_type) > 1:
                    feedCode[num_type[1]] = int(num_type[0])
                    if num_type[1] in BoaConfig.goodFeedList:
                        goodType.append(int(num_type[0]))

            if not feedCode:
                feedCode['AC'] = 1
                feedCode['DC'] = 2
                goodType = [1,2]

            self.Offsets       = zeros((2,nChannels),Float32)
            self.Gain          = ones((nChannels),Float32)
            self.ExtGain       = ones((nChannels),Float32)
            if reader.getType()=="ApexMBFitsReader":
                # Read ALL the offsets, non sky feeds (having -9999 offsets) will be set to 0
                offsets = (array([reader.read("FeedOffX",febe=febe),\
                                  reader.read("FeedOffY",febe=febe)])*3600.).astype(Float32)
                if (rank(offsets) == 1):  # when only one pixel
                    offsets = transpose(array([offsets]))

                # Store array pointing reference position explicitly
                try:
                    if (refChannel > 0):
                        refOffX = offsets[0,refChannel-1]
                        refOffY = offsets[1,refChannel-1]
                    # MBFits 1.65 and higher feature the array reference position keywords
                    # for arrays which do not have a feed at that position.
                    if (refChannel == -1):
                        refOffX = float(reader.read("RefOffX",febe=febe))*3600.
                        refOffY = float(reader.read("RefOffY",febe=febe))*3600.

                        # Calculate pseudo reference channel for display purposes

                        # Valid positions
                        usedFeedoffX = take(offsets[0], allUsedChannels-1)
                        usedFeedoffY = take(offsets[1], allUsedChannels-1)
                        # Distance to reference position
                        refDist = sqrt((usedFeedoffX-refOffX)**2 + (usedFeedoffY-refOffY)**2)
                        # Index of minimum
                        minIndex = argmin(refDist)
                        # Approximate reference feed
                        refChannel = allUsedChannels[minIndex]
                        origRefChannel = -1

                except Exception, e:
                    print 'Error determining array reference position. Assuming (0,0). Exception: %s' % (e)
                    refOffX = 0.0
                    refOffY = 0.0

                # FEBEPAR_FLATFIEL contains one array of gains
                self.Gain  = reader.read("FlatField",febe=febe)[indexBaseband]

                # DC offsets (reseted before scan start)
                self.DCOff = reader.read("DCoffset",febe=febe)[indexBaseband]
            else:  # IramMBFitsReader
                self.updateRCP('Mambo120.rcp')

            # FE and BE Gains
            gains = reader._readGains(febe=febe)
            self.FEGain   = gains[0]
            if type(gains[1]) == type(array([0])):
                self.BEGain   = gains[1][indexBaseband]
            else:
                self.BEGain   = gains[1]
                
            self.JyPerCount    = reader.read("BolCalFc", febe=febe)
                        
            self.UsedChannels  = usedChannels
            self.CurrChanList  = usedChannels
            self.NUsedChannels = len(usedChannels)
            self.RefChannel    = refChannel
            self.RefOffX       = refOffX
            self.RefOffY       = refOffY
            
            # remember 0 indexed numbering
            usedChannels=usedChannels-1

            self.FlagHandler = BoaFlagHandler.createFlagHandler(zeros((nChannels),Int32))
            
            # Setting the initial Flags
            # in two steps: flag all, and unflag the UsedChannels...

            #   ... first flag all with -1 meaning not used ...
            for num in arrayrange(nChannels):
                if 0 <= num < nChannels:
                    self.FlagHandler.setOnIndex(num, flag)

            #   ... second unflag the used channels
            if type(feedType[0]) == type(1):
                for i in arrayrange(self.NUsedChannels):
                    num = usedChannels[i]
                    if 0 <= num < nChannels and feedType[i] in goodType:
                        self.FlagHandler.unsetOnIndex(num, flag)
                    elif 0 <= num < nChannels:
                        self.FlagHandler.setOnIndex(num, flag)
                        # put the offsets of non sky feeds to 0
                        offsets[0,num] = 0
                        offsets[1,num] = 0
                    else:
                        self.__MessHand.warning(str(num+1)+' channel not in the range')
                self.Offsets  = offsets
            else:  # MAMBO data: feedType are already strings
                for i in arrayrange(self.NUsedChannels):
                    num = usedChannels[i]
                    if 0 <= num < nChannels and feedType[i] == 'A':
                        self.FlagHandler.unsetOnIndex(num, flag)
                
            # Check that the reference channel is used
            if self.FlagHandler.isSetOnIndex(self.RefChannel-1): 
                self.__MessHand.warning('Reference channel not used')
                
            self.FeedCode = feedCode
            self.FeedType = feedType

            # The Channel separations can be computed     
            self._BolometerArray__computeChanSep()
            self.FeBe = febe
            # Store the telescope properties
            self.Telescope.set(reader.read("Telescope",      febe=febe),\
                               float(reader.read("Diameter", febe=febe)),\
                               float(reader.read("SiteLong", febe=febe)),\
                               float(reader.read("SiteLat",  febe=febe)),\
                               float(reader.read("SiteElev", febe=febe)))

            # And observing frequency
            freq = float(reader.read("RestFreq",febe=febe,baseband=baseband,subsnum=subscan))
            if (freq == 0.):
                freq = 2.5e11   # assume 1.2mm if not available
                self.MessHand.warning('No frequency found in the file, assuming 1.2mm')
            self.EffectiveFrequency = freq
            
            # The telescope beamSize can be computed
            self.__computeBeamSize()

            # Dewar rotation angles
            try:
                # Old data (e.g. ExampleData) does not have these keywords
                self.DewUser  = reader.read("DewUser", febe=febe)
                self.DewExtra = reader.read("DewExtra", febe=febe,subsnum=subscan)
            except:
                print "No DEWAR angle - probably too old MB-FITS"

        except Exception, data:
            raise

    # ---------------------------------------------------------------------

    def get(self,dataType,flag=[],getFlagged=0):
        """
        DES: get bolometers offsets or gain according to flag
        INP: (string)   dataType : type of data
             (integer list) flag : retrieve data flagged or unflagged accordingly
             (log)    getFlagged : flag revers to flagged/unflagged data
                                   flag   | getFlagged | Retrieve..
                                   'None' |  0         | all data
                                   []     |  0         | unflagged data (default)
                                   []     |  1         | data with at least one flag set
                                   1      |  0         | data with flag 1 not set
                                   1      |  1         | data with flag 1 set
                                   [1,2]  |  0         | data with neither flag 1 nor flag 2 set
                                   [1,2]  |  1         | data with either flag 1 or flag 2 set
        OUT: (float array)       : the requested data
        """

        # retrieve the data... (offsets are in arcsec)
        if dataType in ['BolometerPositionX','BoloPosX','BoloX','bX']:
            dataArray = self.Offsets[0,::]
        elif dataType in ['BolometerPositionY','BoloPosY','BoloY','bY']:
            dataArray = self.Offsets[1,::]
        elif dataType in ['BolometerGain','Gain','gain']:
            dataArray = self.Gain

        # .. and only return the desired flag
        if flag in ['', 'None']:
            return dataArray
        else:
            if getFlagged:
                mask = self.FlagHandler.isSetMask(flag)
            else:
                mask = self.FlagHandler.isUnsetMask(flag)
            return compress(mask, dataArray)

    # ---------------------------------------------------------------------
    def flipOffsets(self):
        """
        DES: flips the sign in Az/El of channel offsets. Used to convert (old) APEX-SZ
             scans into the same convention as for LABOCA
        INP:         
        """
        self.Offsets=-self.Offsets

    # ---------------------------------------------------------------------
    def readAsciiRcp(self,filename='boa.rcp'):
        """
        DES: update receiver channel offsets from a simple ascii file
             channelNumber AzOffset ElOffset Major(FWHM) Minor(FWHN) Tilt Gain
                           with unit of arcsec and degree
        INP: (string) filename: the filename to read in
        """

        try:
            f = file(os.path.join(BoaConfig.rcpPath,filename))
        except IOError:
            self.__MessHand.error("could not open file %s"%(file))
            return

        # Read the file and put the values in 1 list and 2 arrays

        asciiFile = f.readlines()
        f.close()

        Number    = []
        AzOffset  = []
        ElOffset  = []
        MinorFWHM = []
        MajorFWHM = []
        nTilt     = []
        nGain     = []

        for i in range(0,len(asciiFile)):
            
            # Skip comment lines
            if asciiFile[i][0] == "!":
                continue

            if asciiFile[i].find("no fit") > 0:
                continue

            # Split lines and fill arrays
            tmp = asciiFile[i].split()
            Number.append(int(tmp[0]))
            AzOffset.append(float(tmp[1]))
            ElOffset.append(float(tmp[2]))
            MinorFWHM.append(float(tmp[3]))
            MajorFWHM.append(float(tmp[4]))
            nTilt.append(float(tmp[5]))
            nGain.append(float(tmp[6]))

        AzOffset  = array(AzOffset)
        ElOffset  = array(ElOffset)
        MinorFWHM = array(MinorFWHM)
        MajorFWHM = array(MajorFWHM)
        nTilt     = array(nTilt)
        nGain     = array(nGain)
    
        # Process the array
        refChan = self.RefChannel                              
        Offsets = self.Offsets
        Gain    = self.Gain

        #FB260307 usually Tilt and FWHM are not instantiated
        if shape(self.Tilt )[0] == 0:
            nc = shape(Gain)[0]
            Tilt = zeros(nc,Float32)
            FWHM = zeros((2,nc),Float32)
        else:
            FWHM    = self.FWHM
            Tilt    = self.Tilt
            
        # If the reference channel is in the list remove its offsets
        # from the other ones and replace them but the one already
        # existing, otherwise we have to assume that they are aligned
        if refChan in Number:
            refIndex = nonzero(equal(Number,refChan))
            print refIndex,AzOffset[refIndex],ElOffset[refIndex] 
            AzOffset = AzOffset - AzOffset[refIndex] + Offsets[0,refChan-1]
            ElOffset = ElOffset - ElOffset[refIndex] + Offsets[1,refChan-1]
            #FB260307 I dont understand this, so here assume 0,0 is reference channel
            #AzOffset = AzOffset - AzOffset[refIndex] 
            #ElOffset = ElOffset - ElOffset[refIndex] 

        # Replace the offsets
        for i in range(len(Number)):
            Offsets[0,Number[i]-1] = AzOffset[i]
            Offsets[1,Number[i]-1] = ElOffset[i]
            FWHM[:,Number[i]-1]    = array([MajorFWHM[i],MinorFWHM[i]],Float32)
            Tilt[Number[i]-1]      = nTilt[i]
            Gain[Number[i]-1]      = nGain[i]
            

        self.Offsets = Offsets
        self.FWHM    = FWHM
        self.Tilt    = Tilt
        self.Gain    = Gain

        # delete local variable
        del AzOffset, ElOffset, MinorFWHM, MajorFWHM, nTilt, nGain
        del Offsets, FWHM, Tilt, Gain

        # recompute separations between channels
        self._BolometerArray__computeChanSep()

    # ---------------------------------------------------------------------
    def writeAsciiRcp(self,rcpFile='boa.rcp'):
        """
        NAM: writeRCPfile (method)
        DES: store current Receiver Channel Parameters (Offsets,
             Gain) to a file with mopsi like format
        INP: (string) rcpFile: complete name of output file
        """
        
        try:
            f = file(os.path.join(BoaConfig.rcpPath,rcpFile),'w')
        except IOError:
            self.__MessHand.error("could not open file %s in write mode"%(rcpFile))
            return

        #FB260307 in case Tilt and FWHM are not dimensioned yet
        if shape(self.Tilt )[0] == 0:
            nc = shape(self.Gain)[0]
            self.Tilt = zeros(nc,Float32)
            self.FWHM = zeros((2,nc),Float32)

        # Write header
        f.write("! Chan # Az/EL offset Major/Minor FWHM Tilt Gain\n")
        # Write parameters for all channels
	for i in range(len(self.Gain)):
            f.write("%i %f %f %f %f %f %f\n"% \
                    (i+1, self.Offsets[0,i],self.Offsets[1,i],\
                     self.FWHM[0,i], self.FWHM[1,i],\
                     self.Tilt[i], self.Gain[i]))
        f.close()

    # ---------------------------------------------------------------------

    def readAszcaRCP_matlab(self,beamfile,calfile):
        """
        NAM: readRCPfile (method)
        DES: update Receiver Channel Parameters for Aszca (attributes Offsets,
                 Gain and ChannelSep) from the content of a file
                 USING CALIBRATION/BEAM PARAMETER FILES FROM THE MATLAB PIPELINE
        INP: (string) beamfile : complete name of file containing beam parameters
        (string) calfile  : complete name of file containing calibrations
        """

  
        b = file(beamfile)
        c = file(calfile)

        # read and process RCP file
        beam = b.readlines()
        b.close()
        cal = c.readlines()
        c.close()
        offX, offY, gain, flat = [], [], [], []   # local lists to store X and Y offsets and Gain
        fwhmx, fwhmy, tilt = [],[],[]
        timeconst = []
        
        for i in range(len(beam)):
            if beam[i][0] != '%':              # skip comments
                tmp = string.split(beam[i])
                offX.append((-1.)*3600.*string.atof(tmp[1]))
                offY.append((-1.)*3600.*string.atof(tmp[2]))
                fwhmx.append(3600.*string.atof(tmp[4]))
                fwhmy.append(3600.*string.atof(tmp[5]))
                tilt.append(string.atof(tmp[6]))

        for i in range(len(cal)):
            if cal[i][0] != '%':              # skip comments
                tmp = string.split(cal[i])
                gain.append(1./string.atof(tmp[2]))
                flat.append(1./string.atof(tmp[2]))        
                
        # set the attributes to default values
        nChannels = len(gain)

        #offX=(-1.)*offX*3600.
        #offY=(-1.)*offY*3600.
        #fwhmx=fwhmx*3600.
        #fwhmy=fwhmy*3600.
        
        self.NChannels     = nChannels
        self.NUsedChannels = nChannels
        self.Offsets       = zeros((nChannels,nChannels),Float32)
        self.Gain          = zeros((nChannels),Float32)
        self.ExtGain       = zeros((nChannels),Float32)
        self.FWHM          = zeros((nChannels,nChannels),Float32)
        self.Tilt          = zeros((nChannels),Float32)
        self.TimeConst     = zeros((nChannels), Float32)

        self.FlagHandler = BoaFlagHandler.createFlagHandler(zeros((nChannels),Int32))

        # By default use all the channels :
        self.UsedChannels = arrayrange(nChannels)+1
        self.CurrChanList = arrayrange(nChannels)+1

        # Remember RCP file are written in arcsec
        self.Offsets = array([offX,offY],Float32)
        self.Gain    = array(gain,Float32)
        self.ExtGain = array(flat,Float32)
        self.FWHM    = array([fwhmx,fwhmy],Float32)
        self.Tilt    = array(tilt,Float32)
        #self.TimeConst  = array(timeconst, Float32)
        
        # recompute separations between channels
        self._BolometerArray__computeChanSep()
        

    
    def readAszcaRCP(self,rcpFile):
        """
        NAM: readAszcaRCP (method)
        DES: update Receiver Channel Parameters for Aszca (attributes Offsets,
             Gain and ChannelSep) from the content of a file.
             Also read beam shape and time constant
        INP: (string) rcpFile: complete name of file to read in
        """

        try:
            f = file(os.path.join(BoaConfig.rcpPath,rcpFile))
        except IOError:
            self.__MessHand.error("could not open file %s"%(rcpFile))
            return

        # read and process RCP file
        param = f.readlines()
        f.close()
        offX, offY, gain, flat = [], [], [], []   # local lists to store X and Y offsets and Gain
        fwhmx, fwhmy, tilt = [],[],[]
        timeconst = []

        for i in range(len(param)):	        # -1: skip last line
            if param[i][0] not in ['!','\n','','#','%']:              # skip comments
                tmp = string.split(param[i])
                offX.append(string.atof(tmp[3]))
                offY.append(string.atof(tmp[4]))
                gain.append(string.atof(tmp[1]))
                flat.append(string.atof(tmp[2]))
                fwhmx.append(string.atof(tmp[5]))
                fwhmy.append(string.atof(tmp[6]))
                tilt.append(string.atof(tmp[7]))
                timeconst.append(string.atof(tmp[9]))
                
        # set the attributes to default values
        nChannels = len(gain)
        
        self.NChannels     = nChannels
        self.NUsedChannels = nChannels
        self.Offsets       = zeros((nChannels,nChannels),Float32)
        self.Gain          = zeros((nChannels),Float32)
        self.ExtGain       = zeros((nChannels),Float32)
        self.FWHM          = zeros((nChannels,nChannels),Float32)
        self.Tilt          = zeros((nChannels),Float32)
        self.TimeConst     = zeros((nChannels), Float32)

        self.FlagHandler = BoaFlagHandler.createFlagHandler(zeros((nChannels),Int32))

        # By default use all the channels :
        self.UsedChannels = arrayrange(nChannels)+1
        self.CurrChanList = arrayrange(nChannels)+1

        # Remember RCP file are written in arcsec
        self.Offsets = array([offX,offY],Float32)
        self.Gain    = array(gain,Float32)
        self.ExtGain = array(flat,Float32)
        self.FWHM    = array([fwhmx,fwhmy],Float32)
        self.Tilt    = array(tilt,Float32)
        self.TimeConst  = array(timeconst, Float32)
        
        # recompute separations between channels
        self._BolometerArray__computeChanSep()
        

    # ---------------------------------------------------------------------
    def readRCPfile(self,rcpFile):
        """
        NAM: readRCPfile (method)
        DES: update Receiver Channel Parameters (attributes Offsets,
             Gain and ChannelSep) from the content of a file.
             Also read beam shape if available
        INP: (string) rcpFile: complete name of file to read in
        """

        try:
            f = file(os.path.join(BoaConfig.rcpPath,rcpFile))
        except IOError:
            self.__MessHand.error("could not open file %s"%(rcpFile))
            return

        # read and process RCP file
        param = f.readlines()
        f.close()
        offX, offY, gain, flat = [], [], [], []   # local lists to store X and Y offsets and Gain
        fwhmx, fwhmy, tilt = [],[],[]
        
        for i in range(len(param)-1):	        # -1: skip last line
            if param[i][0] != '!':              # skip comments
                tmp = string.split(param[i])
                offX.append(string.atof(tmp[3]))
                offY.append(string.atof(tmp[4]))
                gain.append(string.atof(tmp[1]))
                flat.append(string.atof(tmp[2]))
                if len(tmp) > 5:
                    if tmp[7]:
                        fwhmx.append(string.atof(tmp[5]))
                        fwhmy.append(string.atof(tmp[6]))
                        tilt.append(string.atof(tmp[7]))
                
        # set the attributes to default values
        nChannels = len(gain)
        
        self.NChannels     = nChannels
        self.NUsedChannels = nChannels
        self.Offsets       = zeros((nChannels,nChannels),Float32)
        self.Gain          = zeros((nChannels),Float32)
        self.ExtGain       = zeros((nChannels),Float32)
        self.FWHM          = zeros((nChannels,nChannels),Float32)
        self.Tilt          = zeros((nChannels),Float32)

        self.FlagHandler = BoaFlagHandler.createFlagHandler(zeros((nChannels),Int32))


        # By default use all the channels :
        self.RefChannel   = 1
        self.UsedChannels = arrayrange(nChannels)+1
        self.CurrChanList = arrayrange(nChannels)+1

        # Remember RCP file are written in arcsec
        self.Offsets = array([offX,offY],Float32)
        self.Gain    = array(gain,Float32)
        self.ExtGain = array(flat,Float32)
        self.FWHM    = array([fwhmx,fwhmy],Float32)
        self.Tilt    = array(tilt,Float32)
        
        # recompute separations between channels
        self._BolometerArray__computeChanSep()
        
    # --------------------------------------------------------------------
    def updateRCP(self,rcpFile,scale=1.,readTimeConst=0):
        """
        NAM: updateRCP
        DES: update only offsets and gains from the content of a file
        INP: (string) rcpFile: complete name of file to read in
             (float)  scale:   scale factor to tune initial guess ASZCA rcp  FB20070324 
        """

        try:
            f = file(os.path.join(BoaConfig.rcpPath,rcpFile))
        except IOError:
            self.__MessHand.error("could not open file %s"%(rcpFile))
            return

        # read and process RCP file
        param = f.readlines()
        f.close()
        offsets = self.Offsets
        gain    = self.Gain
        flat    = self.ExtGain

        #FB260307 if Tilt and FWHM are not dimensioned
        if shape(self.Tilt )[0] == 0:
            nc = shape(gain)[0]
            Tilt = zeros(nc,Float32)
            FWHM = zeros((2,nc),Float32)
        else:
            FWHM    = self.FWHM
            Tilt    = self.Tilt

        if readTimeConst:
            nc = shape(gain)[0]
            TimeConst = zeros(nc, Int32)

        flagList = []
        for i in range(len(param)):
            if param[i][0] not in ['!','#','n']:    # skip comments
                tmp = string.split(param[i])
                num = string.atoi(tmp[0])
                offsets[0,num-1] = string.atof(tmp[3])
                offsets[1,num-1] = string.atof(tmp[4])
                gain[num-1]      = string.atof(tmp[1])
                flat[num-1]      = string.atof(tmp[2]) 
                if len(tmp) > 7:
                    FWHM[0,num-1]= string.atof(tmp[5])
                    FWHM[1,num-1]= string.atof(tmp[6])
                    Tilt[num-1]  = string.atof(tmp[7])
                    if len(tmp) > 10:
                        if string.atoi(tmp[10]):
                            flagList.append(num)
	            else:
		        if len(tmp) > 9:
                            if string.atoi(tmp[9]):
                                flagList.append(num)
                if readTimeConst:
                    TimeConst[num-1] = string.atoi(tmp[8])
                    

        self.Offsets   = offsets*scale
        self.Gain      = gain
        self.ExtGain   = flat
        self.FWHM      = FWHM
        self.Tilt      = Tilt
        if readTimeConst:
            self.TimeConst = TimeConst
        
        # recompute separations between channels
        self._BolometerArray__computeChanSep()
        #if flagList:
        #    self._flagChannels(chanList=flagList,flag=1)
        return flagList
    
    # ---------------------------------------------------------------------
    def writeRCPfile(self,rcpFile='rcpBoa.rcp'):
        """
        NAM: writeRCPfile (method)
        DES: store current Receiver Channel Parameters (Offsets,
             Gains, Beam shape) to a file with mopsi like format
        INP: (string) rcpFile: complete name of output file
        """
        
        try:
            f = file(os.path.join(BoaConfig.rcpPath,rcpFile),'w')
        except IOError:
            self.__MessHand.error("could not open file %s in write mode"%(rcpFile))
            return

        #FB260307 in case Tilt and FWHM are not dimensioned yet
        if shape(self.Tilt )[0] == 0:
            nc = shape(self.Gain)[0]
            self.Tilt = zeros(nc,Float32)
            self.FWHM = zeros((2,nc),Float32)

	f.write("!"+rcpFile+"\n")
        # Write parameters for all channels
	for i in range(len(self.Gain)):
            f.write("%i %f %f %f %f %f %f %f \n"%(i+1,\
		self.Gain[i],self.ExtGain[i],\
		self.Offsets[0,i],self.Offsets[1,i],\
                self.FWHM[0,i],self.FWHM[1,i],self.Tilt[i] ))
 
	f.write("! nb bolo...")
	# TO DO: write last line as in rcp files
        f.close()

    #----------------------------------------------------------------------------
    def readAdditionnalIndexFile(self,indexFile='match.dat', refColumn=0, indexColumn=1, comment='!'):
        """
        DES: Read a list of additional index from an ASCII file, to be used with selectAdditionnalIndex()
        INP:   indexFile : the name of the file to read the ...
               refColumn : the column of channel number and ... (default 0, the first column)
             indexColumn : the column to match the channel with (default 1, the second column)
                 comment : comment character (default '!')
        """
        
        try:
            f = file(indexFile)
        except IOError:
            self.__MessHand.error("could not open file %s"%(indexFile))
            return

        # read and process RCP file
        fileContent = f.readlines()
        f.close()

        boloIndex, addIndex = [], []
        
        for i in range(len(fileContent)):
            # skip comments
            if fileContent[i][0] == comment:
                continue
            
            tmp = string.split(fileContent[i])
            boloIndex.append(string.atoi(tmp[refColumn]))
            addIndex.append(tmp[indexColumn])

        s_addIndex = [None]*self.NChannels
        
        # Reorder if necessary 
        for i in range(len(boloIndex)):
            s_addIndex[boloIndex[i]-1] = addIndex[i]
            
        self.AddIndex = s_addIndex
        
    #----------------------------------------------------------------------------
    def selectAdditionnalIndex(self,value=None):
        """
        DES: Select according to the additionnal Index
        INP:  (s) value : the value to test 
        """
        
        addIndex = self.AddIndex
        List = []
        for index in range(len(addIndex)):
            if addIndex[index] == str(value):
                List.append(index+1)        

        # Check for observed channels
        if List != []:
            List = self.checkChanList(List)
        
        self.__MessHand.info("selected %i channel(s) with index eq %s"%(len(List),value)) 

        return List

    #----------------------------------------------------------------------------
    def checkChanList(self,inList,flag=[],getFlagged=0):
        """
        DES: Return a list of valid channels
        INP: (int list/string) inList: list of channel numbers to get, or
             empty list to get the complete list of unflagged channels, or
             'all' or 'al' or 'a' to get the complete list of channels
             (integer list) flag : retrieve data flagged or unflagged accordingly
             (log)    getFlagged : flag revers to flagged/unflagged data
                                   flag   | getFlagged | Retrieve..
                                   'None' |  0         | all data
                                   []     |  0         | unflagged data (default)
                                   []     |  1         | data with at least one flag set
                                   1      |  0         | data with flag 1 not set
                                   1      |  1         | data with flag 1 set
                                   [1,2]  |  0         | data with neither flag 1 nor flag 2 set
                                   [1,2]  |  1         | data with either flag 1 or flag 2 set
        OUT: (int list) list of channel numbers
        """

        chanList = []
        if inList in ['all','al','a']:
             inList = self.UsedChannels
        elif inList == []:
             inList = self.CurrChanList
        elif type(inList) == int:
             inList = [inList]

        UsedChannels = self.UsedChannels

        if flag in ['','None']:
            for num in inList:
                if num in UsedChannels:
                    chanList.append(num)
        else:
            if getFlagged:
                for num in inList:
                    if (num in UsedChannels) and \
                       self.FlagHandler.isSetOnIndex(num-1, flag):
                        chanList.append(num)
            else:
                for num in inList:
                    if num in UsedChannels and \
                       self.FlagHandler.isUnsetOnIndex(num-1, flag):
                        chanList.append(num)
        return array(chanList)

    #----------------------------------------------------------------------------
    def setCurrChanList(self,chanList='?'):
        """
        DES: set list of channels to be treated
        INP: (int list/string) chanList = list of channels, or string '?'
             to get current list of channels, or string 'a' or 'al' or 'all'
             to set current list to all possible channels. Default: '?'
        """
        if chanList in ["all","al","a"]:
            self.CurrChanList = self.checkChanList('all')
        elif type(chanList) == list:
            self.CurrChanList = self.checkChanList(chanList)
        else:
            self.__MessHand.error("channel list required ")

        self.__MessHand.info("selected channels = " + self.printCurrChanList())
        
    #----------------------------------------------------------------------------
    def printCurrChanList(self):
        """
        DES: print the current channel list in somehow "clever" way 
        OUT: a string representing the current channel list
        """

        return prettyPrintList(self.CurrChanList)

    #----------------------------------------------------------------------------
    def plotArray(self,overplot=0,num=0,limitsX=[],limitsY=[],ci=3):
        """
        DES: plot the receiver parameters
        INP: (optional) overplot (logical) = overplot?
             (optional) num (logical) = indicate chan numbers?
        """

        offsetsX, offsetsY = self.Offsets
        refChannel         = self.RefChannel
        refOffX            = self.RefOffX
        refOffY            = self.RefOffY
        halfBeamSize       = self.BeamSize/2
        offsetsX           = offsetsX - refOffX
        offsetsY           = offsetsY - refOffY

        #TODO plot circles which diameter is propto gain
        #TODO plot ellipses instead of circle
        
        if not overplot:
            if not limitsX:
                limitsX = [min(offsetsX),max(offsetsX)]
            if not limitsY:
                limitsY = [min(offsetsY),max(offsetsY)]

            # empty plot
            Plot.plot(offsetsX,offsetsY, \
                      limitsX=limitsX, limitsY=limitsY, \
                      caption = self.FeBe, \
                      labelX ="\gD Az ['']", labelY="\gD El ['']",\
                      aspect=1,nodata=1)        

        # plot the flagged channel in grey
        mask = self.FlagHandler.isSetMask()
        offsetFlagged = tolist_boa(compress(mask, array([offsetsX,offsetsY])))
        if len(offsetFlagged[0]) != 0:
            for i in range(len(offsetFlagged[0])):
                Forms.circle(offsetFlagged[0][i],\
                             offsetFlagged[1][i],\
                             halfBeamSize,ci=14)

        # overplot used channels
        for i in self.UsedChannels:
            if self.FlagHandler.isUnsetOnIndex(i-1):
                # RefChannel in red
                if i == refChannel:
                    Forms.circle(offsetsX[i-1],offsetsY[i-1], halfBeamSize, ci=2)
                else:
                    Forms.circle(offsetsX[i-1],offsetsY[i-1], halfBeamSize, ci=ci)

        # overplot the channels number for the UsedChannels only
        if num:
            for i in self.UsedChannels:
                if self.FlagHandler.isUnsetOnIndex(i-1):
                    Plot.xyout(offsetsX[i-1],offsetsY[i-1],str(i))

    #----------------------------------------------------------------------------
    def plotGain(self, style='idl4'):
        """
        DES: plot the gain of the Array
        INP: (str) style : the style to be used (default idl4)
        WAR: the bolometer without know offsets should be flagged 
        """

        # This function needs a special treatment.
        # Do not delete unless you know what you are doing.

        import interp
        NCP = 4
        bolometerPositionX = self.get('BolometerPositionX')
        bolometerPositionY = self.get('BolometerPositionY')
        bolometerGain      = self.get('BolometerGain')
        nX = 100
        nY = 100
        X = arange(nX,typecode=Float)/(nX)*\
            (max(bolometerPositionX)-min(bolometerPositionX))\
            +min(bolometerPositionX)
        Y = arange(nY,typecode=Float)/(nY)*\
            (max(bolometerPositionY)-min(bolometerPositionY))\
            +min(bolometerPositionY)
        image = interp.interpsf(NCP,\
                                bolometerPositionX,\
                                bolometerPositionY,\
                                bolometerGain,\
                                X,Y)

        # Flag data outside of the array, assume the array is embeded in a circle
        image = ravel(image)
        imageX = transpose(array([X]*nY))
        imageY = array([Y]*nX)
        dist = ravel(sqrt(imageX**2+imageY**2))
        maxDist = max(sqrt(bolometerPositionX**2+bolometerPositionY**2))
        image = where(greater_equal(maxDist,dist),image,float('Nan'))
        image = reshape(image,[nX,nY])
        

        Plot.draw(image,sizeX=[min(X),max(X)],\
                  sizeY=[min(Y),max(Y)],\
                  aspect=1,wedge=1,nan=1,caption = self.FeBe, \
                  labelX ="\gD Az ['']", labelY="\gD El ['']", style=style)

        self.plotArray(overplot=1)
        
    #----------------------------------------------------------------------------
    def __computeBeamSize(self):
        """
        DES: Compute the beam size in arcsec
        """
        
        c0 = 299792458.0	                        #speed of light m/s
        
        self.BeamSize = 1.22*c0/(self.EffectiveFrequency)/ \
                        self.Telescope.Diameter/ \
                        (pi/180.0/3600)                 # airy pattern in arcsec

    #----------------------------------------------------------------------------
    def getChanSep(self,chanList=[]):
        """
        DES: return the channel separation in both direction from the reference channel
        """
        # cast to be sure
        chanList = array(chanList)-1
        
        boloAz = take(self.Offsets[0,::],chanList)-self.RefOffX
        boloEl = take(self.Offsets[1,::],chanList)-self.RefOffY
        return boloAz, boloEl
        
    #----------------------------------------------------------------------------
    def __computeChanSep(self):
        """
        DES: Compute separation between pixels (in arcsec)
        """
        
        nc = self.NChannels
        ChannelSep = zeros((nc,nc),Float)
        Parameters = self.Offsets
        if (nc > 1): # to avoid crash for single-pixel receivers
            for i in range(nc):
                dx = Parameters[0,i] - Parameters[0,::]
                dy = Parameters[1,i] - Parameters[1,::]
                ChannelSep[i,::] = sqrt(dx**2 + dy**2)
                
        self.ChannelSep = ChannelSep

    #----------------------------------------------------------------------------
    def fourpixels(self):
        """
        DES: returns a list of 4 non-flagged channel numbers, selected as follows:
             - the reference channel
             - the two closest neighbours to the ref
             - the furthest one
        """
        ok = self.checkChanList([])  # non-flagged channels
        ok = ok - 1       # python numbering
        ref = self.RefChannel - 1
        refSep = self.ChannelSep[ref,::] # distances to ref chan
        refSep = take(refSep,ok)     # keep only unflagged ones
        result = [self.RefChannel]
        # now find the nearest neighbour
        minsep = max(refSep)
        minNum = 0
        for i in range(len(ok)):
            if refSep[i] < minsep and refSep[i] > 0:
                minsep = refSep[i]
                minNum = ok[i]+1
        result.append(minNum)
        # do it again to find the next neighbour and the furthest away
        minsep2 = max(refSep)
        minNum2 = 0
        for i in range(len(ok)):
            if refSep[i] < minsep2 and refSep[i] > minsep:
                minsep2 = refSep[i]
                minNum2 = ok[i]+1
            if refSep[i] == max(refSep):
                maxNum = ok[i]+1
        result.append(minNum2)
        result.append(maxNum)
        return result
        
    #----------------------------------------------------------------------------
    def flag(self,chanList=[],flag=1):
        """
        DES: assign flags to a list of channels
        INP: (integer array) chanList : list the channels to be flaged
             (integer list)      flag : flag values (default 1)
        """

        # cast to be sure
        chanList = array(chanList)-1
        
        nChannels = self.NChannels

        if self.RefChannel-1 in chanList and flag !=0 :
            self.__MessHand.warning('Reference channel flagged')

        countFlag = 0
        flaggedChan = []

        for num in chanList:
            if 0 <= num < nChannels and self.FlagHandler.isUnsetOnIndex(num, flag):
                self.FlagHandler.setOnIndex(num, flag)
                countFlag += 1
                flaggedChan.append(num+1)
        self.__MessHand.info(' %i channel(s) flagged (%s) with flag %s'%\
                             (countFlag,str(flaggedChan),str(flag)))

        self.__MessHand.debug("Flags="+str(self.FlagHandler.getFlags()))

        return countFlag
    #----------------------------------------------------------------------------

    def unflag(self,chanList=[],flag=[]):
        """
        DES: unflags a list of channels
        INP: (integer array) chanList : list of channels to be unflaged (default all)
             (integer list)      flag : flag values (default []: unset all flags)
        """
        if chanList == []:
            chanList = arrayrange(self.NChannels)
        else:
            chanList = array(chanList)-1

        nChannels = self.NChannels

        countFlag = 0
        flaggedChan = []

        for num in chanList:
            if 0 <= num < nChannels and self.FlagHandler.isSetOnIndex(num, flag):
                self.FlagHandler.unsetOnIndex(num, flag)
                countFlag += 1
                flaggedChan.append(num+1)
        self.__MessHand.info(' %i channel(s) unflagged (%s) with flag %s'%\
                             (countFlag,str(flaggedChan),str(flag)))

        self.__MessHand.debug("Flags="+str(self.FlagHandler.getFlags()))

        return countFlag

    #----------------------------------------------------------------------------

    def rotateArray(self,angle):
        """
        DES: rotate array offsets by a given angle
        INP: (float) angle (in degree)
        """
        # TODO: rotate w.r.t optical axis
        angle = angle*pi/180.
        rotMatrix = array([[cos(angle),-1.*sin(angle)],[sin(angle),cos(angle)]],'f')
        nc = self.NChannels
        for i in range(nc):
            self.Offsets[:,i] = dot(rotMatrix,self.Offsets[:,i])

    def rotateDewar(self):
        """
        DES: rotate array using dewar rotation angle
        """
        self.rotateArray(self.DewUser+self.DewExtra)
        
    #----------------------------------------------------------------------------
    def getChanIndex(self,chanList=[]):
        """
        DES: convert from physical channel number to index in UsedChannel
        INP: (i list) chanList : the physical channel number
        OUT: (i list )           the corresponding index (-1 if failed)
        """

        if type(chanList) == int:
            chanList = [chanList]

        indexing      = range(self.NUsedChannels)
        usedChannels  = self.UsedChannels

        outIndexes = []
        notUsedChannel = []
        
        for chan in chanList:
            out, n = fUtilities.icompress(indexing,usedChannels,chan)
            if n == 1: 
                outIndexes.append(out[0])
            else:
                outIndexes.append(-1)
                notUsedChannel.append(chan)
                
        if notUsedChannel:
            self.__MessHand.error("channel %s not used"%str(notUsedChannel))

        return array(outIndexes)


#----------------------------------------------------------------------------------
#----- ScanParameter Class -------------------------------------------------------
#----------------------------------------------------------------------------------

class ScanParameter:
    """
    NAM: ScanParameter (class)
    DES: Define all parameters (coordinates, time) for a scan
    """

    def __init__(self):
        """
        DES: Instanciation of a new ScanParameter object
        """

        self.__MessHand=BoaMessageHandler.MessHand(self.__module__)

        # General parameters about the SCAN type and geometry
        self.ScanNum  = 0
        self.DateObs  = ""
        self.ScanType = ""
        self.ScanMode = ""
        self.ScanDir  = []  # will be a list of strings
        self.LineLen  = 0.0
        self.LineYsp  = 0.0
        self.AzVel    = 0.0

        # Object
        self.Object   = ""        # Object name
        self.Equinox  = 2000.0    # Default Equinox J2000
        self.Basis    = ("", "")  # Astronomical Basis Frame should be ('ALON-SFL', 'ALAT-SFL')
        self.Coord    = (0.0,0.0) # object coordinates in the basis frame
        self.Frames   = ""        # basis + user frames, e.g. "EQEQHO"

        # Wobbler
        self.WobUsed    = 0
        self.WobMode    = ""
        self.WobThrow   = 0.0
        self.WobCycle   = 0.0
        self.WobblerSta = []     # LIST of strings (strings don't support arrays)
        self.WobblerPos = array([],Float32)
        # this will contain pairs of On-Off integration numbers, if wobbler used
        self.OnOffPairs  = []

        # Pointing and focus status at scan start
        self.Nula     = 0.0
        self.Nule     = 0.0
        self.Colstart = 0.0

        # Total number of integrations:
        self.NInt     = 0
        
        # Focus positions: X, Y, Z, XTILT, YTILT
        self.FocX    = array([],Float32)
        self.FocY    = array([],Float32)
        self.FocZ    = array([],Float32)
        self.PhiX    = array([],Float32)
        self.PhiY    = array([],Float32)
    
        # Telescope coordinates in User and Basis systems
        self.LST      = array([],Float32)
        self.MJD      = array([],Float)
        self.UT       = array([],Float32)
        self.Az       = array([],Float32)   # This is always absolute Az
        self.El       = array([],Float32)   # This is always absolute El

        # added 20080707 MN
        self.MeanRa   = array([],Float32)
        self.MeanDec  = array([],Float32)

        self.LonOff   = array([],Float32)   # Offsets in User native frame
        self.LatOff   = array([],Float32)
        self.BasLon   = array([],Float32)   # Absolute positions in Astron. basis frame
        self.BasLat   = array([],Float32)
        self.Rot      = array([],Float32)
        self.LonPole  = array([],Float32)
        self.LatPole  = array([],Float32)
        self.NoddingState = array([],Int0)  # array of integers
        self.AddLonWT = 0       # add wobbler throw to get right azimuth offset?
        self.AddLatWT = 0       # add wobbler throw to get right elev. offset?
        # Galactic coordinates
        self.GLon     = array([],Float32)
        self.GLat     = array([],Float32)
        self.GalAngle = array([],Float32)
        
        # Offsets in horizontal and equatorial systems
        self.AzOff    = array([],Float32)
        self.ElOff    = array([],Float32)
        self.RAOff    = array([],Float32)
        self.DecOff   = array([],Float32)
        # Source position in equatorial system
        self.RA0      = 0.
        self.Dec0     = 0.
        # Telescope positions in equatorial system (source + offsets)
        self.RA       = array([],Float32)
        self.Dec      = array([],Float32)
        # Parallactic angle = rotation between HO and EQ
        self.ParAngle = array([],Float32)

        # Source coordinates in basis frame for a moving target
        self.Mcrval1  = array([],Float32)
        self.Mcrval2  = array([],Float32)

        # Arrays representing 'subscans'
        self.NObs         = 0                 # Number of subscans in dataset
        self.SubscanNum   = []                # list of integers
        self.SubscanIndex = array([],Int)     # Start/end indices of subscans
        self.SubscanType  = []                # List of strings, ('REF', 'ON', 'OFF'...)
        self.SubscanTime  = []

        # Refraction: one value per subscan
        self.Refraction   = []

        # Ambient temperature - use only value at start
        self.T_amb        = 273.
        
        # He3 temperature: values and timestamps
        self.TimeHe3      = []  # timestamps in Monitor
        self.TempHe3      = []  # corresponding values
        self.He3Temp      = []  # interpolated to MJD(data)

        # Wind speed and direction:
        self.WindSpeed    = [] # m/s
        self.WindDir      = []
        self.TimeWind     = []

        #pwv
        self.PWV          = [] # mm

        # Bias aplitude, QDAC amplitude, and Bias potsetting to calculate Bias voltage:
        self.BiasAmplitude  = [] # Numeric array per subscan
        self.QdacAmplitude  = [] # Numeric array per subscan
        self.BiasPotsetting = [] # Numeric array per subscan
        self.GainSetting    = [] 
        
        # Flags by integration (independant of channels)
        self.FlagHandler = BoaFlagHandler.createFlagHandler(array([], Int32))
        self.iflags = {'TURNAROUND'                   : 1,
                       'ACCELERATION THRESHOLD'       : 2,
                       'ELEVATION VELOCITY THRESHOLD' : 3,
                       'SUBSCAN FLAGGED'              : 7,
                       'TEMPORARY'                    : 8,
                       'BLANK DATA'                   : 9 }
        
    #----------------------------------------------------------------------------
    def __str__(self):
        """
        DES: Defines a string, shown when the print instruction is used.
        """

        out = "Scan number %s : %s  %s on %s contains %3i subscan(s), %5i/%5i records"%\
              (str(self.ScanNum),\
               self.ScanMode,\
               self.ScanType,\
               self.Object,\
               len(self.SubscanNum),\
               self.FlagHandler.nUnset(), \
               len(self.LST))
    
        if BoaConfig.DEBUG > 2:
            out += "\n" + \
                   attrStr(self,['_ScanParameter__MessHand']) + \
                   "\n"
                    
        return out

    # ---------------------------------------------------------------------
    # Overload addition operator: used to combine two datasets
    # ---------------------------------------------------------------------
    def __add__(self,other):
        result = copy.deepcopy(self)
        result._coadd(other)
        return result
    
    def _coadd(self,other):
        # TODO: check that it makes sense to co-add these two datasets
        # TODO: check for missing keywords/attributes!!!
        #
        self.LST           = concatenate((self.LST,         other.LST))
        self.MJD           = concatenate((self.MJD,         other.MJD))
        self.Az            = concatenate((self.Az,          other.Az))
        self.El            = concatenate((self.El,          other.El))
        self.RA            = concatenate((self.RA,          other.RA))
        self.Dec           = concatenate((self.Dec,         other.Dec))
        self.MeanRa        = concatenate((self.MeanRa,      other.MeanRa))
        self.MeanDec       = concatenate((self.MeanDec,     other.MeanDec))
        self.AzOff         = concatenate((self.AzOff,       other.AzOff))
        self.ElOff         = concatenate((self.ElOff,       other.ElOff))
        self.RAOff         = concatenate((self.RAOff,       other.RAOff))
        self.DecOff        = concatenate((self.DecOff,      other.DecOff))
        self.LonOff        = concatenate((self.LonOff,      other.LonOff))
        self.LatOff        = concatenate((self.LatOff,      other.LatOff))
        self.BasLon        = concatenate((self.BasLon,      other.BasLon))
        self.BasLat        = concatenate((self.BasLat,      other.BasLat))
        self.LonPole       = concatenate((self.LonPole,     other.LonPole))
        self.LatPole       = concatenate((self.LatPole,     other.LatPole))
        self.Rot           = concatenate((self.Rot,         other.Rot))
        self.ParAngle      = concatenate((self.ParAngle,    other.ParAngle))
        self.NoddingState  = concatenate((self.NoddingState,other.NoddingState))
        self.UT            = concatenate((self.UT,          other.UT))
        self.WobblerPos    = concatenate((self.WobblerPos,  other.WobblerPos))
        self.FocX          = concatenate((self.FocX,        other.FocX))
        self.FocY          = concatenate((self.FocY,        other.FocY))
        self.FocZ          = concatenate((self.FocZ,        other.FocZ))
        self.PhiX          = concatenate((self.PhiX,        other.PhiX))
        self.PhiY          = concatenate((self.PhiY,        other.PhiY))

        slfFlags  = concatenate((self.FlagHandler.getFlags(), other.FlagHandler.getFlags()))
        self.FlagHandler = BoaFlagHandler.createFlagHandler(slfFlags)

        # Update subscans-related infos - specific case, some attributes are lists
        self.SubscanNum.extend(other.SubscanNum)
        self.SubscanType.extend(other.SubscanType)
        # SubscanIndex: take into account that we merged two datasets
        # integ number of the last point in dataset 1
        max1 = self.SubscanIndex[1,-1]
        self.SubscanIndex = transpose(concatenate((transpose(self.SubscanIndex),
                                                   transpose(other.SubscanIndex+max1))))
        self.NInt = self.NInt + other.NInt
        self.NObs = self.NObs + other.NObs

    #----------------------------------------------------------------------------
    def caption(self):
        """
        DES: Return a short caption of the scan
        """
        out = "Scan: %i (%s) -* %s *- %s"%\
              (self.ScanNum, prettyPrintList(self.SubscanNum), self.Object, self.DateObs)
        return out

    # ---------------------------------------------------------------------
    def __fillFromMBFits(self,reader,febe,baseband,subscans,flag=9,\
                       readHe=0,readAzEl0=0,readT=0,readWind=0,readBias=0,readPWV=0):
        """
        DES: fill a ScanParam object using the MBFitsReader object reader.
        INP:            reader : MBFitsReader object
                    febe (str) : frontend-backend name to select
                baseband (int) : baseband number to select
             subscans (i list) : list of subscans numbers to read in
                 flag (i list) : flag for blanked integrations (default: 9 'BLANK DATA')
              (logical) readHe, readAzEl0, readT, readWind, readBias : see DataEntity.read
        """

        self.__MessHand.debug('start of ScanParam.fillFromMBfits')
        try:
            # Some infos about the SCAN - do not change from one obs to the other
            self.ScanNum  = reader.read("ScanNum")
            self.DateObs  = reader.read("DateObs")
            self.ScanType = reader.read("ScanType")
            self.ScanMode = reader.read("ScanMode",subsnum=subscans[0])

            self.Object   = reader.read("Object")
            self.Equinox  = reader.read("Equinox")
            self.Basis    = reader.read("Basis")
            self.Coord    = reader.read("Coord")

            # Description of Basis and User frames
            self.Frames   = reader.read("UsrFrame",subsnum=subscans[0])

            WobUsed = reader.read("WobUsed")
            if (not WobUsed) or (WobUsed == 'F'):
                self.WobUsed  = 0
                self.WobCycle = 0.
                self.WobThrow = 0.
                self.WobMode  = ""
            else:
                self.WobUsed  = 1
                self.WobCycle = float(reader.read("WobCycle"))
                self.WobThrow = reader.read("WobThrow")
                self.WobMode  = reader.read("WobMode")

            #self.Nula     = reader.read("Nula",febe=febe)
            #self.Nule     = reader.read("Nule",febe=febe)
            #self.Colstart = reader.read("Colstart",febe=febe)

            self.DeltaCA  = float(reader.read("DeltaCA"))
            self.DeltaIE  = float(reader.read("DeltaIE"))

            # Time differences between UT1, UTC and TAI
            self.TAIUTC = float(reader.read("TAIUTC"))
            self.UTCUT1 = float(reader.read("UTCUT1"))
            
            nIntegSubscan = {}
            nIntegTotal = 0
            for subscan in subscans:
                subscanWasOpened = reader.openSubscan(subsnum=subscan)
                
                nIntegSubscan[subscan] = reader.read("NInteg", \
                                                     subsnum=subscan, \
                                                     febe=febe, \
                                                     baseband=baseband)
                nIntegTotal += nIntegSubscan[subscan]

                if subscanWasOpened:
                    reader.closeSubscan(subsnum=subscan)
            
            self.NInt = nIntegTotal
            self.NObs = reader.read("NObs")
            
            LST     = zeros((nIntegTotal),Float32)
            MJD     = zeros((nIntegTotal),Float)   # Use Float64 for MJD
            Az      = zeros((nIntegTotal),Float32)
            El      = zeros((nIntegTotal),Float32)
            MeanRa  = zeros((nIntegTotal),Float32)
            MeanDec = zeros((nIntegTotal),Float32)
            LonOff  = zeros((nIntegTotal),Float32)
            LatOff  = zeros((nIntegTotal),Float32)
            BasLon  = zeros((nIntegTotal),Float32)
            BasLat  = zeros((nIntegTotal),Float32)
            LonPole = zeros((nIntegTotal),Float32)
            LatPole = zeros((nIntegTotal),Float32)
            Rot     = zeros((nIntegTotal),Float32)
            Mcrval1 = zeros((nIntegTotal),Float32)
            Mcrval2 = zeros((nIntegTotal),Float32)
            # Focus positions
            FocX    = zeros((nIntegTotal),Float32)
            FocY    = zeros((nIntegTotal),Float32)
            FocZ    = zeros((nIntegTotal),Float32)
            PhiX    = zeros((nIntegTotal),Float32)
            PhiY    = zeros((nIntegTotal),Float32)

            # Subscan info:
            SubscanNum  = []
            SubscanType = []
            SubscanDir  = []

            # compute index for 1st subscan start and end using python-style numbering:
            # starting at zero, and ending at nb_elements since last element is excluded
            # WARNING: in fortran modules, will have to use [start+1, end]
            subscan_start = [0]
            subscan_end = [nIntegSubscan[subscans[0]]]

            # Now fill local arrays one subscan after the other
            tmpLen = 0

            for subscan in subscans:
                nbData = nIntegSubscan[subscan]
                if nbData:
                    subscanWasOpened = reader.openSubscan(subsnum=subscan)
                    
                    if tmpLen:  # means it's not the 1st subscan
                        subscan_start.append(subscan_end[-1])
                        subscan_end.append(subscan_end[-1]+nbData)

                    MJD[tmpLen:tmpLen+nbData] = reader.read("MJD",
                                                            subsnum=subscan,febe=febe).astype(Float)

                    LST[tmpLen:tmpLen+nbData]     = reader.read("LST",
                                                                subsnum=subscan,
                                                                febe=febe).astype(Float32)
                    Az[tmpLen:tmpLen+nbData]      = reader.read("Az",
                                                                subsnum=subscan,
                                                                febe=febe).astype(Float32)
                    El[tmpLen:tmpLen+nbData]      = reader.read("El",
                                                                subsnum=subscan,
                                                                febe=febe).astype(Float32)

                    try:
                        MeanRa[tmpLen:tmpLen+nbData]  = reader.read("Ra",
                                                                    subsnum=subscan,
                                                                    febe=febe).astype(Float32)
                        MeanDec[tmpLen:tmpLen+nbData] = reader.read("Dec",
                                                                    subsnum=subscan,
                                                                    febe=febe).astype(Float32)
                    except:
                        MeanRa=[]
                        MeanDec=[]
                    
                    LonOff[tmpLen:tmpLen+nbData]  = reader.read("LonOff",
                                                                subsnum=subscan,
                                                                febe=febe).astype(Float32)
                    LatOff[tmpLen:tmpLen+nbData]  = reader.read("LatOff",
                                                                subsnum=subscan,
                                                                febe=febe).astype(Float32)
                    BasLon[tmpLen:tmpLen+nbData]  = reader.read("BasLon",
                                                                subsnum=subscan,
                                                                febe=febe).astype(Float32)
                    BasLat[tmpLen:tmpLen+nbData]  = reader.read("BasLat",
                                                                subsnum=subscan,
                                                                febe=febe).astype(Float32)
                    Rot[tmpLen:tmpLen+nbData]     = reader.read("RotAngle",
                                                                subsnum=subscan,
                                                                febe=febe).astype(Float32)

                    if self.WobUsed:
                        self.WobblerSta.extend(reader.read("WobblerSta",subsnum=subscan,febe=febe))
                    
                    # Refraction
                    if self.Frames[:2] == 'HO':
                        self.Refraction.append(reader.read("Refract",subsnum=subscan))
                    else:
                        self.Refraction.append([0])

                    # Ambient temperature at scan start
                    if readT and subscan == subscans[0]:
                        self.T_amb = reader.read("T_amb",subsnum=subscan)
                        self.T_amb = 273.15 + self.T_amb[0]
                    if reader.getType() == "ApexMBFitsReader":
                        # Read data not present in IRAM MB-FITS:
                        LonPole[tmpLen:tmpLen+nbData] = reader.read("LonPole",
                                                                    subsnum=subscan,
                                                                    febe=febe).astype(Float32)
                        LatPole[tmpLen:tmpLen+nbData] = reader.read("LatPole",
                                                                    subsnum=subscan,
                                                                    febe=febe).astype(Float32)

                        # Coordinates for a moving target
                        Mcrval1[tmpLen:tmpLen+nbData] = reader.read("MVAL1",
                                                                    subsnum=subscan,
                                                                    febe=febe).astype(Float32)
                        Mcrval2[tmpLen:tmpLen+nbData] = reader.read("MVAL2",
                                                                    subsnum=subscan,
                                                                    febe=febe).astype(Float32)

                        # Focus positions
                        FocX[tmpLen:tmpLen+nbData]    = reader.read("FocX",
                                                                    subsnum=subscan,
                                                                    febe=febe).astype(Float32)
                        FocY[tmpLen:tmpLen+nbData]    = reader.read("FocY",
                                                                    subsnum=subscan,
                                                                    febe=febe).astype(Float32)
                        FocZ[tmpLen:tmpLen+nbData]    = reader.read("FocZ",
                                                                    subsnum=subscan,
                                                                    febe=febe).astype(Float32)
                        PhiX[tmpLen:tmpLen+nbData]    = reader.read("PhiX",
                                                                    subsnum=subscan,
                                                                    febe=febe).astype(Float32)
                        PhiY[tmpLen:tmpLen+nbData]    = reader.read("PhiY",
                                                                    subsnum=subscan,
                                                                    febe=febe).astype(Float32)

                        if readHe:
                            # He3 temperature - WARNING: specific to LABOCA!
                            tmpHe3, timeHe3 = reader.read("TempHe",subsnum=subscan)
                            # Interpolate to the data timestamps
                            mjd     = 86400.*(MJD[tmpLen:tmpLen+nbData] - MJD[0])
                            timeHe3 = 86400.*(timeHe3 - MJD[0])
                            tmpInterp = arrayfns.interp(tmpHe3,timeHe3,mjd)

                            self.TimeHe3.extend(timeHe3)
                            self.TempHe3.extend(tmpHe3)
                            self.He3Temp.extend(tmpInterp)

                        if readBias:
                            biasAmplitude  = reader.read("BiasAmplitude", subsnum=subscan)
                            qdacAmplitude  = reader.read("QdacAmplitude", subsnum=subscan)
                            biasPotsetting = reader.read("BiasPotsetting", subsnum=subscan)
                            gainSetting    = reader.read("GainSetting", subsnum=subscan)
                            self.BiasAmplitude.append(biasAmplitude)
                            self.QdacAmplitude.append(qdacAmplitude)
                            self.BiasPotsetting.append(biasPotsetting)
                            self.GainSetting.append(gainSetting)

                        if readWind:
                            windSpeed, windDir, timeWind = reader.read("WindSpeedDir",
                                                                       subsnum=subscan)
                            self.WindSpeed.append(windSpeed)
                            self.WindDir.append(windDir)
                            self.TimeWind.append(timeWind)

                        if readPWV:
                            pwv=reader.read("PWV",subsnum=subscan)
                            self.PWV.append(pwv)

                    ######################################################
                    # Fill the subscans related fields
                    SubscanNum.append(subscan)
                    # Subscan type (REF, ON, OFF)
                    SubscanType.append(reader.read("ObsType",subsnum=subscan,febe=febe))
                    # Subscan direction
                    SubscanDir.append(reader.read("ScanDir",subsnum=subscan,febe=febe))
                    # LST time at subscan start
                    self.SubscanTime.append(LST[tmpLen])
                
                    # Compute float value of subscan start time, using SLALIB
                    date_tmp = reader.read("SubsStart",
                                           subsnum=subscan,
                                           febe=febe) # string, e.g. 2003-06-08T20:41:21
                    if date_tmp:
                        year, day, status = sla_clyd(string.atof(date_tmp[0:4]),
                                                     string.atof(date_tmp[5:7]),
                                                     string.atof(date_tmp[8:10]))
                        date_flt = year + (day-1.)/365.
                    else:
                        date_flt = 0.

                    tmpLen += nbData

                    # The following is needed for pointing scans, to generate a .dat
                    # file that is used to compute pointing model
                    if readAzEl0 and subscan == subscans[0]:
                        antenna, encoder = reader.read("AzEl0",subsnum=subscans[0])
                        self.EncAz0 = encoder[0]
                        self.EncEl0 = encoder[1]
                        self.AntAz0 = antenna[0]
                        self.AntEl0 = antenna[1]
                        # when these are requested, we also need the following:
                        self.PDeltaCA = reader.read("PDeltaCA")
                        self.PDeltaIE = reader.read("PDeltaIE")
                        self.FDeltaCA = reader.read("FDeltaCA")
                        self.FDeltaIE = reader.read("FDeltaIE")

                    # Everything has been read, close the file
                    if subscanWasOpened:
                        reader.closeSubscan(subsnum=subscan)
                
                else:
                    # nbData = 0: subscan not readable (no data)
                    self.__MessHand.warning("Subscan %i not readable (no data)"%(subscan))

            # The following attributes are for the full scan
            NoddingState = zeros(nIntegTotal,Int0)
            UT           = zeros(nIntegTotal,Float32)
            WobblerPos   = zeros(nIntegTotal,Float32)
            SubscanIndex = array([subscan_start,subscan_end],Int)

            # At reading flag the data with blanked LST
            self.FlagHandler = BoaFlagHandler.createFlagHandler(zeros(nIntegTotal,Int32))
            flagValue = reader.getBlankFloat()
            mask = equal(LST,flagValue)
            self.FlagHandler.setOnMask(mask, self.iflags['BLANK DATA'])
            
            self.LST          = fUtilities.as_column_major_storage(LST)
            self.MJD          = fUtilities.as_column_major_storage(MJD)
            self.Az           = fUtilities.as_column_major_storage(Az)
            if MeanRa:
                self.MeanRa       = fUtilities.as_column_major_storage(MeanRa)
                self.MeanDec      = fUtilities.as_column_major_storage(MeanDec)
            self.El           = fUtilities.as_column_major_storage(El)
            self.LonOff       = fUtilities.as_column_major_storage(LonOff)
            self.LatOff       = fUtilities.as_column_major_storage(LatOff)
            self.BasLon       = fUtilities.as_column_major_storage(BasLon)
            self.BasLat       = fUtilities.as_column_major_storage(BasLat)
            self.LonPole      = fUtilities.as_column_major_storage(LonPole)
            self.LatPole      = fUtilities.as_column_major_storage(LatPole)
            self.Rot          = fUtilities.as_column_major_storage(Rot)
            self.NoddingState = fUtilities.as_column_major_storage(NoddingState)
            self.UT           = fUtilities.as_column_major_storage(UT)
            self.WobblerPos   = fUtilities.as_column_major_storage(WobblerPos)
            self.SubscanIndex = fUtilities.as_column_major_storage(SubscanIndex)
            self.FocX         = fUtilities.as_column_major_storage(FocX)
            self.FocY         = fUtilities.as_column_major_storage(FocY)
            self.FocZ         = fUtilities.as_column_major_storage(FocZ)
            self.PhiX         = fUtilities.as_column_major_storage(PhiX)
            self.PhiY         = fUtilities.as_column_major_storage(PhiY)
            self.SubscanNum   = SubscanNum
            self.SubscanType  = SubscanType
            self.ScanDir      = SubscanDir

            self.Mcrval1      = fUtilities.as_column_major_storage(Mcrval1)
            self.Mcrval2      = fUtilities.as_column_major_storage(Mcrval2)

            del LST, MJD, Az, El, LonOff, LatOff, BasLon, BasLat, \
                LonPole, LatPole, Rot, NoddingState, \
                UT, WobblerPos, SubscanIndex, \
                FocX, FocY, FocZ, PhiX, PhiY

            if self.NInt:
                # at least one subscan successfully read
                # Compute the missing coordinates / offsets
                self.computeRa0De0()
                if self.Frames[-2::] == 'HO':            # user frame = horizontal
                    self.AzOff = copy.copy(self.LonOff) # then we already have the HO offsets
                    self.ElOff = copy.copy(self.LatOff)
                else: 
                    self.computeAzElOffsets()
                if self.Frames[:2] == 'EQ':   # basis frame = equatorial
                    if fStat.f_rms(self.BasLon,fStat.f_mean(self.BasLon)):
                        # if not constant
                        self.RA  = copy.copy(self.BasLon) # RA, Dec already available
                        self.Dec = copy.copy(self.BasLat)
                else:
                    # changed 20080717 MN
                    # if basis fram is horizontal, try to read ra, dec from
                    # meanRa, meanDec (not all MBFits versions)
                    if (self.MeanRa) and (self.MeanDec):
                        self.RA  = copy.copy(self.MeanRa)
                        self.Dec = copy.copy(self.MeanDec)
                    else:
                        self.__MessHand.warning('No equatorial coordinte information')
                        self.__MessHand.warning('Computing RA/Dec from telescope coordinates')
                        self.computeRaDec()
                if self.Frames[-2::] == 'EQ':             # user frame = equatorial
                    self.RAOff = copy.copy(self.LonOff)  # then we already have the EQ offsets
                    self.DecOff = copy.copy(self.LatOff)
                else:
                    self.computeRaDecOffsets()
                self.computeParAngle()

            else:
                self.__MessHand.warning("No subscans readable, no data")

        except Exception, data:
            raise

    #----------------------------------------------------------------------------
    def he3SmoothInterpolate(self,flag=[], getFlagged=0):
        """
        DES: this is a *function* which *returns* an array with He3 temperatures
             interpolated to the data timestamps, with a smoothing (boxcar window
             applied) before interpolating
        INP: (integer list) flag : retrieve data flagged or unflagged accordingly
             (log)    getFlagged : flag revers to flagged/unflagged data
                                   flag   | getFlagged | Retrieve..
                                   'None' |  0         | all data
                                   []     |  0         | unflagged data (default)
                                   []     |  1         | data with at least one flag set
                                   1      |  0         | data with flag 1 not set
                                   1      |  1         | data with flag 1 set
                                   [1,2]  |  0         | data with neither flag 1 nor flag 2 set
                                   [1,2]  |  1         | data with either flag 1 or flag 2 set
        OUT: (f array) interpolated He3 temperatures are returned
        """
        if not len(self.TimeHe3):
            self.__MessHand.error("No He3 temperature available - returning")
            self.__MessHand.error("You should use 'read(<scan>,readHe=1)'")
            return
        
        he3time = array(self.TimeHe3)  # timestamps in Monitor table
        he3temp = array(self.TempHe3)  # corresponding values
        mjd = self.get('mjd',flag=flag,getFlagged=getFlagged)
        
        # smooth he3temp monitor points
        nb = int(max(he3time))+1
        newx = array(range(nb),'f')
        newy = zeros((nb),'f')

        # compute max. of delta time
        tmptime = he3time[1::]-he3time[:-1]
        tmptime = tmptime.tolist()
        tmptime.append(2.*he3time[0])  # delta time between start and 1st datapoint
        deltatime = int(max(tmptime)/2.+1)

        # smoothing
        for i in range(nb):
            mask = nonzero(less(abs(he3time - newx[i]),deltatime))
            newy[i] = fStat.f_mean(take(he3temp,mask))

        # interpolate he3temp to a regular time grid
        dt = fStat.f_median (mjd[1::]-mjd[:-1])  # delta time
        n = int(max(mjd)/dt)+1
        tt = dt*array(range(n),'f')  # regular timestream covering the scan
        yy = arrayfns.interp(newy,newx,tt)

        # finally interpolate he3temp to original time stamps
        he3T = arrayfns.interp(yy,tt,mjd)
        return he3T
    
    #----------------------------------------------------------------------------
    #----------------------------------------------------------------------------
    def computeRa0De0(self):
        """
        DES: compute source coordinates in equatorial system
        """
        self.__MessHand.debug('start of ScanParam.computeRa0De0')
        if self.Frames[:2] == 'EQ':   # basis frame = equatorial
            self.RA0  = self.Coord[0]
            self.Dec0 = self.Coord[1]
        elif self.Frames[:2] == 'HO': # basis frame = horizontal
            az0 = self.Coord[0] * pi/180.
            #el0 = self.Coord[1] * pi/180.
            el0  = self.Coord[1] * pi/180. - self.Refraction[0][0]*pi/180.
            phi = self.Telescope.Latitude * pi/180.
            ha, dec = sla_h2e(az0,el0,phi)
            # last = gmst + data.BolometerArray.Telescope.Longitude * pi/180. + sla_eqeqx(mjd)
            # + UT1-TAI...
            lst0 = self.LST[0] * 360./86400. * pi/180.
            ra = lst0 - ha # in equinox = obs. date, precess to J2000
            date = self.DateObs
            year,day,status = sla_calyd(int(date[:4]),int(date[5:7]),int(date[8:10]))
            ra_0,de_0 = sla_preces ('FK5',float(year)+float(day)/365.25,2000.,ra,dec)
            self.RA0  = ra_0 * 180./pi
            self.Dec0 = de_0 * 180./pi
                         
        else:   # other systems not supported yet
            self.__MessHand.warning('Unsupported astronomical basis frame: '+self.Frames[:2])
            self.__MessHand.warning('RA, Dec of the source not computed')
            
    #----------------------------------------------------------------------------
    def computeAzElOffsets(self):
        """
        DES: compute telescope Az, El offsets w.r.t. the source, using antenna
             Az, El and RA, Dec of the source
        """
        self.__MessHand.debug('start of ScanParam.computeAzElOffsets')
        nb = self.NInt
        AzOff = zeros((nb),Float32)
        ElOff = zeros((nb),Float32)
        ra0 = self.RA0 * pi/180.
        de0 = self.Dec0 * pi/180.
        # Precess to date of obs.
        mjd0 = self.MJD[0]
        epoch = sla_epj(mjd0)
        ra0,de0 = sla_map(ra0,de0,0,0,0,0,2000,mjd0)  # precess + mean to apparent

        TAI_TT  = 32.184
        TAI_UTC = -1.*self.TAIUTC
        UTC_UT1 = self.UTCUT1
        TAI_UT1 = TAI_UTC + UTC_UT1
        beta= self.Telescope.Longitude * pi / 180.
        phi = self.Telescope.Latitude * pi / 180.

        # convert RA0, Dec0 to Az0, El0 for each time stamp
        for i in xrange(nb):
            mjd = self.MJD[i]   # this is MJD TAI
            mjd_ut1 = mjd + TAI_UT1/86400.  # MJD UT1
            mjd_tt = mjd + TAI_TT/86400.    # MJD TT, = TDB within a few ms
            gmst = sla_gmsta(int(mjd_ut1), mjd_ut1-int(mjd_ut1))
            lst = gmst + beta + sla_eqeqx(mjd_tt)
            ha = lst - ra0
            az0, el0 = sla_e2h(ha,de0,phi)
            AzOff[i] = self.Az[i]*pi/180. - az0
            AzOff[i] = sla_drange(AzOff[i])
            AzOff[i] = AzOff[i] * 180./pi * cos(self.El[i]*pi/180.)
            ElOff[i] = self.El[i]*pi/180. - el0
            ElOff[i] = sla_drange(ElOff[i]) * 180./pi
        self.AzOff = fUtilities.as_column_major_storage(AzOff)
        self.ElOff = fUtilities.as_column_major_storage(ElOff)
                    
    #----------------------------------------------------------------------------
    def computeRaDec(self):
        """
        DES: compute telescope RA, Dec positions from Az, El
        """
        self.__MessHand.debug('start of ScanParam.computeRaDec')
        nb = self.NInt
        ra  = zeros((nb),Float32)
        dec = zeros((nb),Float32)
        az  = self.Az * pi/180.
        el  = self.El * pi/180. - self.Refraction[0][0]*pi/180.
        phi = self.Telescope.Latitude*pi/180.
        mjd0 = self.MJD[0]

        # convert Az, El to RA, Dec for each time stamp
        for i in xrange(nb):
            lst = self.LST[i] * 360./86400. * pi/180. # lst in radians
            ha, dec0 = sla_dh2e(az[i],el[i],phi) # radians, radians
            ra0 = lst - ha # in equinox = obs. date, precess to J2000
            ra_0,de_0 = sla_amp(ra0,dec0,mjd0,2000.)  # precess + apparent to mean
            ra[i]  = ra_0 * 180./pi
            dec[i] = de_0 * 180./pi
        self.RA  = fUtilities.as_column_major_storage(ra)
        self.Dec = fUtilities.as_column_major_storage(dec)
                
    #----------------------------------------------------------------------------
    def computeRaDecOffsets(self):
        """
        DES: compute telescope RA, Dec offsets w.r.t. the source
        """
        self.__MessHand.debug('start of ScanParam.computeRaDecOffsets')
        self.RAOff = (self.RA - self.RA0)*cos(self.Dec*pi/180)
        self.DecOff = self.Dec - self.Dec0
        
    #----------------------------------------------------------------------------
    def computeParAngle(self):
        """
        DES: compute parallactic angle
        """
        self.__MessHand.debug('start of ScanParam.computeParAngle')
        nb = self.NInt
        parAng = zeros((nb),Float32)
        lst = self.LST * 360./86400. * pi/180.
        ra  = self.RA * pi/180.
        ha = lst - ra
        dec = self.Dec * pi / 180.
        phi = self.Telescope.Latitude*pi/180.
        for i in xrange(nb):
            ang = sla_pa (ha[i], dec[i], phi)
            parAng[i] = ang * 180./pi  # in degrees
        self.ParAngle = fUtilities.as_column_major_storage(parAng)
        
    #----------------------------------------------------------------------------
    def computeGal(self):
        """
        DES: compute telescope GLon, GLat positions from RA, Dec
        """
        self.__MessHand.debug('start of ScanParam.computeGal')
        nb = self.NInt
        gl  = zeros((nb),Float32)
        gb = zeros((nb),Float32)
        ra  = self.RA * pi/180.
        dec = self.Dec * pi/180.
        # convert RA,Dec to gl,gb for each time stamp
        for i in xrange(nb):
            tmpgl,tmpgb = sla_eqgal(ra[i],dec[i])
            gl[i] = tmpgl * 180./pi
            if gl[i] > 180.:
                gl[i] = gl[i]-360.
            gb[i] = tmpgb * 180./pi
        self.GLon = fUtilities.as_column_major_storage(gl)
        self.GLat = fUtilities.as_column_major_storage(gb)
                
    #----------------------------------------------------------------------------
    def computeGalAngle(self):
        """
        DES: compute angle EQ to GAL
        """
        self.__MessHand.debug('start of ScanParam.computeGalAngle')
        nb = self.NInt
        galAng = zeros((nb),Float32)
        ra  = self.RA * pi/180.
        dec = (self.Dec+1.) * pi / 180.
        gl  = self.GLon * pi / 180.
        gb  = self.GLat * pi / 180.
        for i in xrange(nb):
            gl1,gb1 = sla_eqgal(ra[i],dec[i])
            if gl1 > pi:
                gl1 = gl1-2.*pi
            ang = arctan((gl1-gl[i])/(gb1-gb[i]))
            galAng[i] = ang * 180./pi  # in degrees
        self.GalAngle = fUtilities.as_column_major_storage(galAng)
        
    #----------------------------------------------------------------------------
    #----------------------------------------------------------------------------
    def flipOffsets(self,system='eq'):
        """
        DES: change sign of telescope offsets w.r.t. reference position
        INP: (string) system = 'eq' or 'ho', to flip RA/Dec offsets or Az/El
                               offsets (default: 'eq')
        """
        if string.upper(system) == 'EQ':
            self.RA = 2.*self.RA0 - self.RA
            self.Dec = 2.*self.Dec0 - self.Dec
            self.RAOff = -1.*self.RAOff
            self.DecOff = -1.*self.DecOff
        elif string.upper(system) == 'HO':
            self.AzOff = -1.*self.AzOff
            self.ElOff = -1.*self.ElOff
        else:
            self.__MessHand.error('Unkown coordinate system')

    #############################################################################
    # Methods related with chopped data
    #----------------------------------------------------------------------------
    def computeOnOff(self):
        """
        DES: determine ON-OFF pairs from content of WobblerSta, and fill
             OnOffPairs attribute with pairs of integration numbers.
             The result is a 2 x Nb_Integ. array of integers.
        """

        # Get number of integrations
        nd = self.NInt
        num1,num2 = 0,0

        if self.OnOffPairs:
            # if already exists, may be of type array
            # => convert to list for appending new pairs
            self.OnOffPairs = tolist_boa(self.OnOffPairs)
        while ((num1 < nd-1) & (num2 < nd-1)):
            # Initialise flags, set when one ON-OFF pair is found
            okOn, okOff = 0,0
            # Find the next (ON,OFF) pair
            while ((okOn == 0) & (num1 < nd-1)):
                if (self.WobblerSta[num1] not in ['ON',1]):
                    num1 += 1
                else:
                    okOn = 1
            num2 = num1+1
            while ((okOff == 0) & (num2 < nd)):
                if (self.WobblerSta[num2] in ['ON',1]):
                    num1 = num2
                    num2 += 1
                else:
                    okOff = 1

            if (okOn & okOff):
                self.OnOffPairs.append([num1,num2])
                num1 += 1
        
        if self.OnOffPairs:
            self.OnOffPairs = array(self.OnOffPairs)
            self.NInt = len(self.OnOffPairs[:,0])
            
        
    #----------------------------------------------------------------------------
    def _phaseDiffParam(self):
        """
        NAM: phaseDiffParam (method)
        DES: Compute the phase differences for data associated parameters.
             Times are average of ON and OFF, coordinates are ON positions.
        """

        # compute times: (time(on)+time(off)) / 2.
        self.LST = (take(self.LST,self.OnOffPairs[:,0]) + \
        			take(self.LST,self.OnOffPairs[:,1]))/2.
        self.MJD = (take(self.MJD,self.OnOffPairs[:,0]) + \
                                take(self.MJD,self.OnOffPairs[:,1]))/2.
        self.UT = (take(self.UT,self.OnOffPairs[:,0]) + \
                                take(self.UT,self.OnOffPairs[:,1]))/2.

        # compute positions = position(on)
        self.Az      = take(self.Az,      self.OnOffPairs[:,0])
        self.El      = take(self.El,      self.OnOffPairs[:,0])
        self.AzOff   = take(self.AzOff,   self.OnOffPairs[:,0])
        self.ElOff   = take(self.ElOff,   self.OnOffPairs[:,0])
        self.RA      = take(self.RA,      self.OnOffPairs[:,0])
        self.Dec     = take(self.Dec,     self.OnOffPairs[:,0])
        self.MeanRa  = take(self.MeanRa,  self.OnOffPairs[:,0])
        self.MeanDec = take(self.MeanDec, self.OnOffPairs[:,0])
        self.RAOff   = take(self.RAOff,   self.OnOffPairs[:,0])
        self.DecOff  = take(self.DecOff,  self.OnOffPairs[:,0])
        self.LonOff  = take(self.LonOff,  self.OnOffPairs[:,0])
        self.LatOff  = take(self.LatOff,  self.OnOffPairs[:,0])
        self.BasLon  = take(self.BasLon,  self.OnOffPairs[:,0])
        self.BasLat  = take(self.BasLat,  self.OnOffPairs[:,0])
        self.LonPole = take(self.LonPole, self.OnOffPairs[:,0])
        self.LatPole = take(self.LatPole, self.OnOffPairs[:,0])
        self.Rot     = take(self.Rot,     self.OnOffPairs[:,0])
        self.ParAngle= take(self.ParAngle,self.OnOffPairs[:,0])
        
        # Focus positions: use average of ON and OFF positions
        self.FocX   = (take(self.FocX,self.OnOffPairs[:,0]) + \
                        take(self.FocX,self.OnOffPairs[:,1]))/2.
        self.FocY   = (take(self.FocY,self.OnOffPairs[:,0]) + \
                        take(self.FocY,self.OnOffPairs[:,1]))/2.
        self.FocZ   = (take(self.FocZ,self.OnOffPairs[:,0]) + \
                        take(self.FocZ,self.OnOffPairs[:,1]))/2.
        self.PhiX   = (take(self.PhiX,self.OnOffPairs[:,0]) + \
                        take(self.PhiX,self.OnOffPairs[:,1]))/2.
        self.PhiY   = (take(self.PhiY,self.OnOffPairs[:,0]) + \
                        take(self.PhiY,self.OnOffPairs[:,1]))/2.

        # Flag: could use the max of flag1 and flag2 (if one is flagged,
        # then the phase diff should also be flagged) - for now, use binary_or:

        slfFlags = self.FlagHandler.getFlags()
        slfFlags = take(slfFlags,self.OnOffPairs[:,0]) | \
                   take(slfFlags,self.OnOffPairs[:,1])
        self.FlagHandler = BoaFlagHandler.createFlagHandler(slfFlags)

    #----------------------------------------------------------------------------
    def selectPhase(self,phase):
        """
        NAM: selectPhase (method)
        DES: Keep only parameters (times, positions) associated with
             Data(ON) or Data(OFF)
        INP: (int) phase: phase to keep, 1=ON, 2=OFF
        """

        ph = phase-1  # index in OnOffPairs: 0 = ON, 1 = OFF
        self.LST     = take(self.LST,     self.OnOffPairs[:,ph])
        self.MJD     = take(self.MJD,     self.OnOffPairs[:,ph])
        self.UT      = take(self.UT,      self.OnOffPairs[:,ph])
        self.Az      = take(self.Az,      self.OnOffPairs[:,ph])
        self.El      = take(self.El,      self.OnOffPairs[:,ph])
        self.MeanRa  = take(self.MeanRa,  self.OnOffPairs[:,ph])
        self.MeanDec = take(self.MeanDec, self.OnOffPairs[:,ph])
        self.LonOff  = take(self.LonOff,  self.OnOffPairs[:,ph])
        self.LatOff  = take(self.LatOff,  self.OnOffPairs[:,ph])
        self.BasLon  = take(self.BasLon,  self.OnOffPairs[:,ph])
        self.BasLat  = take(self.BasLat,  self.OnOffPairs[:,ph])
        self.LonPole = take(self.LonPole, self.OnOffPairs[:,ph])
        self.LatPole = take(self.LatPole, self.OnOffPairs[:,ph])
        self.Rot     = take(self.Rot,     self.OnOffPairs[:,ph])
        self.FocX    = take(self.FocX,    self.OnOffPairs[:,ph])
        self.FocY    = take(self.FocY,    self.OnOffPairs[:,ph])
        self.FocZ    = take(self.FocZ,    self.OnOffPairs[:,ph])
        self.PhiX    = take(self.PhiX,    self.OnOffPairs[:,ph])
        self.PhiY    = take(self.PhiY,    self.OnOffPairs[:,ph])
        #self.Flags   = take(self.Flags,   self.OnOffPairs[:,ph])

        slfFlags = self.FlagHandler.getFlags()
        slfFlags = take(slfFlags, self.OnOffPairs[:,ph])
        self.FlagHandler = BoaFlagHandler.createFlagHandler(slfFlags)


    def __computeSubIndex(self):
        """
        DES: Compute start and end indices per subscan
        """

        self.__MessHand.debug('computeSubIndex start...')
        for i in range(len(self.SubscanTime)):
            # build array of indices where LST > LST_0
            after = nonzero(greater(self.LST,self.SubscanTime[i]))
            self.SubscanIndex[0,i] = after[0]

        self.SubscanIndex[1,0:-1] = self.SubscanIndex[0,1::]
        self.SubscanIndex[1,-1] = len(self.LST)
        self.__MessHand.debug('... computeSubIndex end')


    #----------------------------------------------------------------------------
    def get(self,dataType=' ', flag=[], getFlagged=0, subscans=[]):
        """
        DES: get data of the ScanParam class
        INP: (string) dataType : type of data
                                 LST MJD Az El AzOff ElOff focX focY focZ 
             (integer list) flag : retrieve data flagged or unflagged accordingly
             (log)    getFlagged : flag revers to flagged/unflagged data
                                   flag   | getFlagged | Retrieve..
                                   'None' |  0         | all data
                                   []     |  0         | unflagged data (default)
                                   []     |  1         | data with at least one flag set
                                   1      |  0         | data with flag 1 not set
                                   1      |  1         | data with flag 1 set
                                   [1,2]  |  0         | data with neither flag 1 nor flag 2 set
                                   [1,2]  |  1         | data with either flag 1 or flag 2 set
             (i list) subscans : optionnally select subscan(s)
        OUT: (float array)     : the requested data

        returned data are in the stored unit except for offsets which are
        converted to arcsec
        """
        
        # retrieve the data... (offsets are in arcsec)
        if dataType in ['LST','lst']:
            dataArray = self.LST
        elif dataType in ['MJD','mjd']:
            dataArray = (self.MJD - self.MJD[0])*86400. # in seconds since the beginning of the scan
        elif dataType in ['Azimuth','azimuth','Az','az']:
            dataArray = self.Az
        elif dataType in ['Elevation','elevation','El','el']:
            dataArray = self.El
        elif dataType in ['AzimuthOffset', 'azimuthoffset', 'AzOff', 'azoff','azo']:
            dataArray = self.AzOff*3600.
            if self.AddLatWT!=0:
                dataArray=dataArray+self.AddLatWT*self.WobThrow
        elif dataType in ['ElevationOffset', 'elevationoffset', 'ElOff','eloff','elo']:
            dataArray = self.ElOff*3600.
            if self.AddLonWT!=0:
                dataArray=dataArray+self.AddLonWT*self.WobThrow
        elif dataType in ['LonOff','Lonoff','lonoff']:
            dataArray = self.LonOff*3600.
        elif dataType in ['LatOff','Latoff','latoff']:
            dataArray = self.LatOff*3600.
        elif dataType in ['UT']:
            dataArray = self.UT
        elif dataType in ['FOCUS-X','FocusX','FocX','focX','focx']:
            dataArray = self.FocX
        elif dataType in ['FOCUS-Y','FocusY','FocY','focY','focy']:
            dataArray = self.FocY
        elif dataType in ['FOCUS-Z','FocusZ','FocZ','focZ','focz']:
            dataArray = self.FocZ
        elif dataType in ['FOCUS-XTILT','Focus-XTILT']:
            dataArray = self.PhiX
        elif dataType in ['FOCUS-YTILT','Focus-YTILT']:
            dataArray = self.PhiY
        elif dataType in ['BasLon','Baslon', 'baslon']:
            dataArray = self.BasLon
        elif dataType in ['BasLat','Baslat', 'baslat']:
            dataArray = self.BasLat
        elif dataType in ['ra','Ra','RA']:
            dataArray = self.RA
        elif dataType in ['dec','Dec','DEC']:
            dataArray = self.Dec
        elif dataType in ['meanRa','MeanRa','meanRA','MeanRA']:
            dataArray = self.MeanRa
        elif dataType in ['meanDec','MeanDec','meanDEC','MeanDEC']:
            dataArray = self.MeanDec
        elif dataType in ['raoff','RaOff','RAOff','raOffset','RaOffset','RAOffset']:
            dataArray = self.RAOff * 3600.
        elif dataType in ['decoff','DecOff','DECOff','decOffset','DecOffset','DECOffset']:
            dataArray = self.DecOff * 3600.
        elif dataType in ['GLon','glon','Glon']:
            dataArray = self.GLon
        elif dataType in ['GLat','glat','Glat']:
            dataArray = self.GLat
        elif dataType in ['HeT','HeTemp']:
            dataArray = self.He3Temp
        elif dataType in ['azspeed','Azspeed','AzSpeed']:
            azimuth   = self.Az
            dt        = self.get('deltat',flag='None')
            elev      = self.El
            dataArray = (azimuth[1::]-azimuth[:-1])*cos(elev[1::]*pi/180.)/dt[:-1]
            dataArray = 3600.*concatenate((dataArray,[dataArray[-1]]))
        elif dataType in ['elspeed','Elspeed','ElSpeed']:
            elevation = self.El
            dt        = self.get('deltat',flag='None')
            dataArray = (elevation[1::]-elevation[:-1])/dt[:-1]
            dataArray = 3600.*concatenate((dataArray,[dataArray[-1]]))
        elif dataType in ['azacc','Azacc','AzAcc']:
            dt        = self.get('deltat',flag='None')
            azspeed   = self.get('azspeed',flag='None')
            dataArray = (azspeed[1::]-azspeed[:-1])/dt[:-1]
            dataArray = concatenate((dataArray,[dataArray[-1]]))
        elif dataType in ['elacc','Elacc','ElAcc']:
            dt        = self.get('deltat',flag='None')
            elspeed   = self.get('elspeed',flag='None')
            dataArray = (elspeed[1::]-elspeed[:-1])/dt[:-1]
            dataArray = concatenate((dataArray,[dataArray[-1]]))
        elif dataType in ['speed','Speed']:
            azimuthSpeed   = self.get('azspeed',flag='None')
            elevationSpeed = self.get('elspeed',flag='None')
            dataArray = sqrt(azimuthSpeed**2+elevationSpeed**2)
        elif dataType in ['acc','Acc','accel','Accel']:
            azimuthAcc   = self.get('azacc',flag='None')
            elevationAcc = self.get('elacc',flag='None')
            dataArray = sqrt(azimuthAcc**2+elevationAcc**2)
        elif dataType in ['deltat','deltaT','DeltaT']:
            MJD  = (self.MJD - self.MJD[0])*86400.
            dt   = MJD[1::]-MJD[:-1]
            bad  = nonzero(equal(dt,0.))
            good = nonzero(not_equal(dt,0))
            if bad:
                for nbad in bad:
                    dt[nbad] = min(take(dt,good))
            dt = concatenate((dt,[min(dt)]))
            dataArray = dt
        elif dataType in ['phase','Phase','wob','wobbler']:
            dataArray = self.WobblerSta
            
        # Check if subscans was asked
        if subscans:
            # Create a mask representing the asked subscans
            subscanMask = zeros(shape(dataArray))
            
            SubscanNum   = self.SubscanNum
            SubscanIndex = self.SubscanIndex
            
            for subscan in subscans:
                if subscan in SubscanNum:
                    isub = SubscanNum.index(subscan)
                    subscanMask[SubscanIndex[0,isub]:SubscanIndex[1,isub]] = 1
                else:
                    self.MessHand.error("subscan "+subscan+" does not exist")
                    return
        else:
            subscanMask = ones(shape(dataArray))
                
        # Extract that mask from the data
        dataArray = compress(subscanMask, dataArray)

        # .. and only return the desired flag 
        if flag in ['', 'None']:
            return dataArray
        else:
            if getFlagged:
                mask = self.FlagHandler.isSetMask(flag)
            else:
                mask = self.FlagHandler.isUnsetMask(flag)
            mask = compress(subscanMask, mask)

            return compress(mask, dataArray)
        
    #----------------------------------------------------------------------------

    def flag(self,dataType='', below='?', above='?', flag=1):
        """
        DES: flag data based on dataType
        INP: (float) below : flag dataType < below (default max)
             (float) above : flag dataType > above (default min)
             (integer list) flag : flag values (default 1)

             below and above should be in unit of the flagged data,
             except for 'Lon' and 'Lat' where they should be in arcsec
        """
        self.__MessHand.debug("flag start...")

        # flag on dataType
        dataTest = self.get(dataType=dataType,flag='None')

        # default inputs
        if above == '?':
            above = min(dataTest)
        if below == '?':
            below = max(dataTest)

        mask = where(bitwise_and(dataTest >= above,dataTest <= below),1,0)

        if len(nonzero(mask)) > 0:
            n0 = self.FlagHandler.nSet(flag)
            self.FlagHandler.setOnMask(mask, flag)
            n1 = self.FlagHandler.nSet(flag)
            if (n1-n0):
                self.__MessHand.info("%5i timestamps flagged (%5.2f %%) with flag %s"
                                     % ((n1-n0) ,100.*float(n1-n0)/self.NInt, str(flag)))
            else:
                self.__MessHand.warning("Nothing flagged")
        else:
            self.__MessHand.warning("Nothing flagged")

        self.__MessHand.debug("flag end")

    #----------------------------------------------------------------------------

    def unflag(self,dataType='', below='?', above='?', flag=[]):
        """
        DES: unflag data based on dataType
        INP: (float) below : unflag dataType < below (default max)
             (float) above : unflag dataType > above (default min)
             (integer list) flag : flag values (default []: unset all flags)

             below and above should be in unit of the flagged data,
             except for 'Lon' and 'Lat' where they should be in arcsec
        """
        self.__MessHand.debug("unflag start...")

        # flag on dataType
        dataTest = self.get(dataType=dataType,flag='None')

        # default inputs
        if above == '?':
            above = min(dataTest)
        if below == '?':
            below = max(dataTest)

        mask = where(bitwise_and(dataTest >= above,dataTest <= below),1,0)

        if len(nonzero(mask)) > 0:
            n0 = self.FlagHandler.nUnset(flag)
            self.FlagHandler.unsetOnMask(mask, flag)
            n1 = self.FlagHandler.nUnset(flag)
            if (n1-n0):
                self.__MessHand.info("%5i timestamps unflagged (%5.2f %%) with flag %s"
                                     % ((n1-n0) ,100.*float(n1-n0)/self.NInt, str(flag)))
            else:
                self.__MessHand.warning("Nothing unflagged")
        else:
            self.__MessHand.warning("Nothing unflagged")

        self.__MessHand.debug("unflag end")

    #----------------------------------------------------------------------------
    def _computeSubscanEfficiency(self):
        """
        DES: compute the scan efficiencies wrt subscans dead time
        OUT: return % of time spent IN subscans wrt total scan length
        """

        # The following variables should contain all the needed information
        # We do not use flagging info, we just compute the SubScan Efficiency
        
        MJD          = self.MJD             # Use MJD because LST could be flagged
        Nobs         = self.NObs
        SubscanIndex = self.SubscanIndex

        ScanDuration = MJD[-1]-MJD[0]
        SubscanDuration = 0
        for i in arange(Nobs):
            SubscanDuration += MJD[SubscanIndex[1,i]-1]-MJD[SubscanIndex[0,i]]

        return SubscanDuration/ScanDuration*100

    #----------------------------------------------------------------------------
    def plotAzimuth(self,flag=[],plotFlagged=0, \
                    limitsX=[],limitsY=[], \
                    style='l',ci=1,overplot=0,aspect=1):
        """
        DES: plot time series of azimuth
        INP: (int list)   flag : plot data flagged or unflagged accordingly
             (log) plotFlagged : flag revers to flagged/unflagged data
                                   flag   | plotFlagged | Plot..
                                   'None' |  0          | all data
                                   []     |  0          | unflagged data (default)
                                   []     |  1          | data with at least one flag set
                                   1      |  0          | data with flag 1 not set
                                   1      |  1          | data with flag 1 set
                                   [1,2]  |  0          | data with neither flag 1 nor flag 2 set
                                   [1,2]  |  1          | data with either flag 1 or flag 2 set
        """

        dataX = self.get('MJD',flag=flag,getFlagged=plotFlagged)
        dataY = self.get('Azimuth',flag=flag,getFlagged=plotFlagged)

        xLabel = "MJD [sec]"
        yLabel = "Az [Deg]"

        Plot.plot(dataX,dataY,\
                  limitsX = limitsX, limitsY = limitsY, \
                  labelX = xLabel, labelY = yLabel, caption=self.caption(), \
                  style=style,ci=ci,overplot=overplot,aspect=aspect)
        
    #----------------------------------------------------------------------------
    def plotElevation(self,flag=[],plotFlagged=0,limitsX=[],limitsY=[],
                      style='l',ci=1,overplot=0,aspect=1):
        """
        DES: plot time series of elevation
        INP: (int list)   flag : plot data flagged or unflagged accordingly
             (log) plotFlagged : flag revers to flagged/unflagged data
                                   flag   | plotFlagged | Plot..
                                   'None' |  0          | all data
                                   []     |  0          | unflagged data (default)
                                   []     |  1          | data with at least one flag set
                                   1      |  0          | data with flag 1 not set
                                   1      |  1          | data with flag 1 set
                                   [1,2]  |  0          | data with neither flag 1 nor flag 2 set
                                   [1,2]  |  1          | data with either flag 1 or flag 2 set
        """

        dataX = self.get('MJD',flag=flag,getFlagged=plotFlagged)
        dataY = self.get('Elevation',flag=flag,getFlagged=plotFlagged)

        xLabel = "MJD [sec]"
        yLabel = "El [Deg]"

        Plot.plot(dataX,dataY,\
                  limitsX = limitsX, limitsY = limitsY, \
                  labelX = xLabel, labelY = yLabel, caption=self.caption(), \
                  style=style,ci=ci,overplot=overplot,aspect=aspect)
        
    #----------------------------------------------------------------------------
    def plotAzEl(self,flag=[],plotFlagged=0,limitsX=[],limitsY=[],style='l',
                 ci=1,overplot=0,aspect=1):
        """
        DES: plot azimuth vs. elevation
        INP: (int list)   flag : plot data flagged or unflagged accordingly
             (log) plotFlagged : flag revers to flagged/unflagged data
                                   flag   | plotFlagged | Plot..
                                   'None' |  0          | all data
                                   []     |  0          | unflagged data (default)
                                   []     |  1          | data with at least one flag set
                                   1      |  0          | data with flag 1 not set
                                   1      |  1          | data with flag 1 set
                                   [1,2]  |  0          | data with neither flag 1 nor flag 2 set
                                   [1,2]  |  1          | data with either flag 1 or flag 2 set
        """

        dataX = self.get('Azimuth',flag=flag,getFlagged=plotFlagged)
        dataY = self.get('Elevation',flag=flag,getFlagged=plotFlagged)
        
        xLabel = "Az [Deg]"
        yLabel = "El [Deg]"

        Plot.plot(dataX,dataY,\
                  limitsX = limitsX, limitsY = limitsY, \
                  labelX = xLabel, labelY = yLabel, caption=self.caption(), \
                  style=style,ci=ci,overplot=overplot,aspect=aspect)
        
    #----------------------------------------------------------------------------
    def plotElevationOffset(self,flag=[],plotFlagged=0,limitsX=[],limitsY=[],
                            style='l',ci=1,overplot=0,aspect=1):
        """
        DES: plot time series of elevation offset
        INP: (int list)   flag : plot data flagged or unflagged accordingly
             (log) plotFlagged : flag revers to flagged/unflagged data
                                   flag   | plotFlagged | Plot..
                                   'None' |  0          | all data
                                   []     |  0          | unflagged data (default)
                                   []     |  1          | data with at least one flag set
                                   1      |  0          | data with flag 1 not set
                                   1      |  1          | data with flag 1 set
                                   [1,2]  |  0          | data with neither flag 1 nor flag 2 set
                                   [1,2]  |  1          | data with either flag 1 or flag 2 set
        """
        
        dataX = self.get('MJD',flag=flag,getFlagged=plotFlagged)
        dataY = self.get('ElevationOffset',flag=flag,getFlagged=plotFlagged)
        
        xLabel = "MJD [sec]"
        yLabel = "\GD El ['']"

        Plot.plot(dataX,dataY,\
                  limitsX = limitsX, limitsY = limitsY, \
                  labelX = xLabel, labelY = yLabel, caption=self.caption(), \
                  style=style,ci=ci,overplot=overplot,aspect=aspect)
    #----------------------------------------------------------------------------
    def plotAzimuthOffset(self,flag=[],plotFlagged=0,limitsX=[],limitsY=[],
                          style='l',ci=1,overplot=0,aspect=1):
        """
        DES: plot time series of azimuth offset
        INP: (int list)   flag : plot data flagged or unflagged accordingly
             (log) plotFlagged : flag revers to flagged/unflagged data
                                   flag   | plotFlagged | Plot..
                                   'None' |  0          | all data
                                   []     |  0          | unflagged data (default)
                                   []     |  1          | data with at least one flag set
                                   1      |  0          | data with flag 1 not set
                                   1      |  1          | data with flag 1 set
                                   [1,2]  |  0          | data with neither flag 1 nor flag 2 set
                                   [1,2]  |  1          | data with either flag 1 or flag 2 set
        """

        dataX = self.get('MJD',flag=flag,getFlagged=plotFlagged)
        dataY = self.get('AzimuthOffset',flag=flag,getFlagged=plotFlagged)
        
        xLabel = "MJD [sec]"
        yLabel = "\gD Az ['']"

        Plot.plot(dataX,dataY,\
                  limitsX = limitsX, limitsY = limitsY, \
                  labelX = xLabel, labelY = yLabel, caption=self.caption(), \
                  style=style,ci=ci,overplot=overplot,aspect=aspect)

    #----------------------------------------------------------------------------
    def plotAzElOffset(self,flag=[],plotFlagged=0,limitsX=[],limitsY=[],
                       style='l',ci=1,overplot=0,aspect=0,caption='',num=1):
        """
        DES: plot elevation offset versus azimuth offset
        INP: (int list)   flag : plot data flagged or unflagged accordingly
             (log) plotFlagged : flag revers to flagged/unflagged data
                                   flag   | plotFlagged | Plot..
                                   'None' |  0          | all data
                                   []     |  0          | unflagged data (default)
                                   []     |  1          | data with at least one flag set
                                   1      |  0          | data with flag 1 not set
                                   1      |  1          | data with flag 1 set
                                   [1,2]  |  0          | data with neither flag 1 nor flag 2 set
                                   [1,2]  |  1          | data with either flag 1 or flag 2 set
        """

        NObs         = self.NObs
        SubscanNum   = self.SubscanNum
        SubscanIndex = self.SubscanIndex
        dataX        = self.get('AzimuthOffset',flag=flag,getFlagged=plotFlagged)
        dataY        = self.get('ElevationOffset',flag=flag,getFlagged=plotFlagged)

        
        if not overplot: 
            xLabel = "\gD Az ['']"
            yLabel = "\gD El ['']"
            if not caption:
                caption = self.caption()
            Plot.plot(dataX,dataY,\
                      limitsX = limitsX, limitsY = limitsY, \
                      labelX = xLabel, labelY = yLabel, caption=caption, \
                      style=style,ci=ci,aspect=aspect,nodata=1)

        #for i in arange(NObs):
        for obs in SubscanNum:
            dataX = self.get('AzimuthOffset',flag=flag,getFlagged=plotFlagged,subscans=[obs])
            dataY = self.get('ElevationOffset',flag=flag,getFlagged=plotFlagged,subscans=[obs])
            if (shape(dataX)[0] > 0):
                Plot.plot(dataX,dataY,style=style,ci=ci,overplot=1)
                if num:
                    Plot.xyout(dataX[0],dataY[0],str(obs))
                      
    #----------------------------------------------------------------------------
    def plotAzElSpeed(self,flag=[],plotFlagged=0,limitsX=[],limitsY=[],
                      style='l',ci=1,overplot=0,aspect=1):
        """
        DES: plot azimuth and elevation speed
        INP: (int list)   flag : plot data flagged or unflagged accordingly
             (log) plotFlagged : flag revers to flagged/unflagged data
                                   flag   | plotFlagged | Plot..
                                   'None' |  0          | all data
                                   []     |  0          | unflagged data (default)
                                   []     |  1          | data with at least one flag set
                                   1      |  0          | data with flag 1 not set
                                   1      |  1          | data with flag 1 set
                                   [1,2]  |  0          | data with neither flag 1 nor flag 2 set
                                   [1,2]  |  1          | data with either flag 1 or flag 2 set
        """
        MJD            = self.get('MJD',flag=flag,getFlagged=plotFlagged)
        azimuthSpeed   = self.get('azspeed',flag=flag,getFlagged=plotFlagged)
        elevationSpeed = self.get('elspeed',flag=flag,getFlagged=plotFlagged)
        
        xLabel = "MJD (sec)"
        yLabel = "|Speed| [''/s]"
        
        MultiPlot.plot(['Azimuth', 'Elevation','Combined'],\
                       [MJD[:-1],MJD[:-1],MJD[:-1]],\
                       [abs(azimuthSpeed[:-1]), abs(elevationSpeed[:-1]),\
                        sqrt(azimuthSpeed[:-1]**2+elevationSpeed[:-1]**2)],\
                       limitsX = limitsX, limitsY = limitsY, \
                       labelX = xLabel, labelY = yLabel, caption=self.caption(), \
                       style=style,ci=ci,overplot=overplot)
    
    #----------------------------------------------------------------------------
    def plotAzElAcceleration(self,flag=[],plotFlagged=0,limitsX=[],limitsY=[],
                             style='l',ci=1,overplot=0,aspect=1):
        """
        DES: plot azimuth and elevation acceleration
        INP: (int list)   flag : plot data flagged or unflagged accordingly
             (log) plotFlagged : flag revers to flagged/unflagged data
                                   flag   | plotFlagged | Plot..
                                   'None' |  0          | all data
                                   []     |  0          | unflagged data (default)
                                   []     |  1          | data with at least one flag set
                                   1      |  0          | data with flag 1 not set
                                   1      |  1          | data with flag 1 set
                                   [1,2]  |  0          | data with neither flag 1 nor flag 2 set
                                   [1,2]  |  1          | data with either flag 1 or flag 2 set
        INP: (int) flag : flag to be plot (default 0 : valid data, -1 plot all)
        """

        MJD       = self.get('MJD',flag=flag,getFlagged=plotFlagged)
        azimuthAccel   = self.get('azacc',flag=flag,getFlagged=plotFlagged)
        elevationAccel = self.get('elacc',flag=flag,getFlagged=plotFlagged)
        
        xLabel = "MJD (sec)"
        yLabel = "|Acceleration| [''/s^2]"
        
        MultiPlot.plot(['Azimuth', 'Elevation','Combined'],\
                       [MJD[:-2],MJD[:-2],MJD[:-2]],\
                       [abs(azimuthAccel), abs(elevationAccel),
                        sqrt(azimuthAccel**2+elevationAccel**2)],\
                       limitsX = limitsX, limitsY = limitsY, \
                       labelX = xLabel, labelY = yLabel, caption=self.caption(), \
                       style=style,ci=ci,overplot=overplot)
    

    #----------------------------------------------------------------------------
    # Subscan related methods
    #----------------------------------------------------------------------------

    def findSubscanByOffset(self,off=60.,combine=10):
        """
        DES: compute subscan indices by looking for sufficient spatial offset
        (in any direction, but in the az/el system)
        INP: (float) off = minimum spatial offset between subscans,
                           in az/el system, in arcseconds
        """

        # get az/el offsets
        azoff = self.get('AzimuthOffset',flag='None')  # in arcsec
        eloff = self.get('ElevationOffset',flag='None')

        SubIndex = [0]
        sizeAzoff = size(azoff)
        
        oldpos=[azoff[0],eloff[0]]

        for i in range(sizeAzoff):
            newpos=[azoff[i],eloff[i]]
            dist=sqrt((oldpos[0]-newpos[0])**2+(oldpos[1]-newpos[1])**2)
            if (dist > off):
                SubIndex.append(i)
                oldpos=copy.deepcopy(newpos)

        # combine subscans found in this way
        NewSubIndex=[0]
        cycle=0
        for index in SubIndex:
            cycle+=1
            if (cycle == combine):
                NewSubIndex.append(index)
                cycle=0

        NewSubIndex.append(size(azoff))
    
        # Reinitialise Subscan related attributes
        self.SubscanIndex = []
        self.SubscanNum   = []
        self.SubscanType  = []

        nb = len(NewSubIndex)-1
        for i in range(nb):
            self.SubscanIndex.append([NewSubIndex[i],NewSubIndex[i+1]]) 
            self.SubscanNum.append(i+1)
            self.SubscanType.append('ON')

        self.SubscanIndex = transpose(self.SubscanIndex)
        self.NObs = nb
        self.__MessHand.info(str("Found %i subscans"%(nb)))


    #----------------------------------------------------
                
    def findSubscanFB(self,azMax=1000.,eq=0):
        """
        DES: compute subscan indices from steps in az, el
        INP: (float) azMax = azimuth offset where subscans are marked
             (logical) eq  - for EQ scan patterns 
        """
        # Retrieve all Az, El offsets (even if some are flagged) to
        # indices refering to the complete dataset
        if eq:
            azOff = self.get('raoff',flag='None')  # in arcsec
            elOff = self.get('decoff',flag='None')
        else:
            azOff = self.get('AzimuthOffset',flag='None')  # in arcsec
            elOff = self.get('ElevationOffset',flag='None')

        SubIndex = [0]
        
        if (abs(azOff[0]) <= azMax): SubPos=[0]
        else:
            if (azOff[0] > 0): SubPos = [1]
            else: SubPos = [-1]

        sizeAzoff = size(azOff)
        inside = 0
        if abs(azOff[0])<azMax:
            inside = 1

        for i in range(sizeAzoff):
            if inside:
                if (abs(azOff[i]) > azMax):
                    inside = 0
                    SubIndex.append(i)
                    if (azOff[i] > 0):
                        SubPos.append(1)
                    else:
                        SubPos.append(-1)
            else:
                if (abs(azOff[i]) <= azMax):
                    inside = 1
                    SubIndex.append(i)
                    SubPos.append(0)
                    
        SubIndex.append(sizeAzoff)
        
        self.SubscanIndex = array([SubIndex[:-1],SubIndex[1::]])
        self.SubscanNum   = []
        self.SubscanType  = []
        nb = shape(self.SubscanIndex)[1]
        self.NObs = nb

        for i in range(nb):
            self.SubscanNum.append(i+1)
            self.SubscanType.append('ON')

        self.SubscanPos = SubPos
        
        self.__MessHand.info(str("Found %i subscans"%(nb)))



    #----------------------------------------------------------------------------
    def findSubscan(self,direction='El',combine=1):
        """
        DES: compute subscan indices for circular scans by looking for sign change in az/el speed
        INP: (string) direction = 'Az' or 'El' - direction in which to look for stationary points
             (int)    combine - number of found subscans to combine into one
                                (useful for irregular scan patterns)
        """

        # scan speed in az/el
        if (direction=='Az' or direction=='az'):
            azsp=self.get('AzSpeed',flag='None')
        else:
            azsp=self.get('ElSpeed',flag='None')

        # replace zeroes with a small epsilon
        eps=0.00001
        azsp0=(where(azsp==0,eps,0))
        azsp=azsp+azsp0
        
        # find where sign changes
        sgn=azsp[1::]*azsp[:-1]
        sgn=sgn/(abs(sgn))
                
        #mask=(nonzero(where(less(sgn,0),1,0)))
        mask=nonzero(where(sgn<0,1,0))

        # make sure there are at least smin=20 timestamps in each subscan
        SubIndex = []
        smin=20
        last=0
        for i in range(len(mask)-1):
            if (mask[i]-last > smin): 
                SubIndex.append(mask[i])
            last=mask[i]

        # combine subscans found in this way
        NewSubIndex=[0]
        cycle=0
        for index in SubIndex:
            cycle+=1
            if (cycle == combine):
                NewSubIndex.append(index)
                cycle=0

        NewSubIndex.append(size(azsp))
    
        # Reinitialise Subscan related attributes
        self.SubscanIndex = []
        self.SubscanNum   = []
        self.SubscanType  = []

        nb = len(NewSubIndex)-1
        for i in range(nb):
            self.SubscanIndex.append([NewSubIndex[i],NewSubIndex[i+1]]) 
            self.SubscanNum.append(i+1)
            self.SubscanType.append('ON')

        self.SubscanIndex = transpose(self.SubscanIndex)
        self.NObs = nb
        self.__MessHand.info(str("Found %i subscans"%(nb)))


    #----------------------------------------------------------------------------

    def findSubscanSpiral(self,threshold=1500.,combine=1):
        """
        DES: compute subscan indices for spiral scans by looking for large acceleration
        INP: (float)  threshold - mark new subscan where acceleration exceeds this value
             (int)    combine - number of found subscans to combine into one
                                (useful for somewhat irregular scan patterns)
        """

        # scan speed in az/el
        acc=self.get('acc',flag='None')

        mask=nonzero(where(acc>threshold,1,0))

        # make sure there are at least smin=400 timestamps in each subscan
        SubIndex = []
        smin=400
        last=0
        for i in range(len(mask)-1):
            if (mask[i]-last > smin): 
                SubIndex.append(mask[i])
            last=mask[i]

        # combine subscans found in this way
        NewSubIndex=[0]
        cycle=0
        for index in SubIndex:
            cycle+=1
            if (cycle == combine):
                NewSubIndex.append(index)
                cycle=0

        NewSubIndex.append(size(acc))
    
        # Reinitialise Subscan related attributes
        self.SubscanIndex = []
        self.SubscanNum   = []
        self.SubscanType  = []

        nb = len(NewSubIndex)-1
        for i in range(nb):
            self.SubscanIndex.append([NewSubIndex[i],NewSubIndex[i+1]]) 
            self.SubscanNum.append(i+1)
            self.SubscanType.append('ON')

        self.SubscanIndex = transpose(self.SubscanIndex)
        self.NObs = nb
        self.__MessHand.info(str("Found %i subscans"%(nb)))


    #----------------------------------------------------------------------------

    def findSubscanByFlagArray(self,gflags):

        gflags=where(gflags > 0, 1, 0)

        flagVal=gflags[0]
        SubFlags=[flagVal]
        SubIndex=[0]
        for i in range(len(gflags)):
            if (gflags[i] != flagVal):
                flagVal=gflags[i]
                SubIndex.append(i)
                SubFlags.append(flagVal)
        SubIndex.append(size(gflags))

        # Reinitialise Subscan related attributes
        self.SubscanIndex = []
        self.SubscanNum   = []
        self.SubscanType  = []

        nb = len(SubIndex)-1
        for i in range(nb):
            self.SubscanIndex.append([SubIndex[i],SubIndex[i+1]]) 
            self.SubscanNum.append(i+1)
            if (SubFlags[i] > 0):
                self.SubscanType.append('OFF')
            else:
                self.SubscanType.append('ON')

        self.SubscanIndex = transpose(self.SubscanIndex)
        self.NObs = nb
        self.__MessHand.info(str("Found %i subscans"%(nb)))
        
    #----------------------------------------------------------------------------    

    def findSubscanCircle(self,combine=1,minLen=400):
        """
        DES: compute subscan indices for ''families of circles''
        INP: 
             (int)    combine - number of found subscans to combine into one
                                (useful for somewhat irregular scan patterns)
             (int)    minLen  - minimum length of one subscan (in time steps)
        """

        # scan speed
        azSpeed   = self.get('azspeed',flag='None',getFlagged=0)
        elSpeed   = self.get('elspeed',flag='None',getFlagged=0)
        speed     = sqrt(azSpeed**2+elSpeed**2)
        medSpeed  = fStat.f_median(speed)
        offSpeed  = abs(speed-medSpeed)
        kern=array(range(100))*0.0+1.0
        a=convolve(offSpeed,kern,mode=1)
        medOff=fStat.f_median(array(a))
        threshold=medOff*7.

        # offset from circle in az/el offset coord ("offset from offset")
        AzEl = array([self.get('AzimuthOffset',flag='None'), \
                          self.get('ElevationOffset',flag='None')])
        azeloff=AzEl[0,::]**2 + AzEl[1,::]**2
        # this only works for (well-behaved) drift circles!!
        scanradius=fStat.f_median(azeloff)
        offcircle = where(bitwise_or((azeloff > scanradius*1.5),(azeloff < scanradius*0.5)),1,0)
        #offcircle = where((azeloff > scanradius*1.2),1,0)
        o=convolve(offcircle,kern,mode=1)

        lrg=where(bitwise_and((a > threshold),(offcircle > 0)),1,0)
        mask=nonzero(lrg)

        # make sure there are at least smin=400 timestamps in each subscan
        SubIndex = []
        smin=minLen
        last=0
        for i in range(len(mask)-1):
            if (mask[i]-last > smin): 
                SubIndex.append(mask[i])
            last=mask[i]

        # combine subscans found in this way
        NewSubIndex=[0]
        cycle=0
        for index in SubIndex:
            cycle+=1
            if (cycle == combine):
                NewSubIndex.append(index)
                cycle=0

        NewSubIndex.append(size(speed))
    
        # Reinitialise Subscan related attributes
        self.SubscanIndex = []
        self.SubscanNum   = []
        self.SubscanType  = []

        nb = len(NewSubIndex)-1
        for i in range(nb):
            self.SubscanIndex.append([NewSubIndex[i],NewSubIndex[i+1]]) 
            self.SubscanNum.append(i+1)
            self.SubscanType.append('ON')

        self.SubscanIndex = transpose(self.SubscanIndex)
        self.NObs = nb
        self.__MessHand.info(str("Found %i subscans"%(nb)))

        ## SET FLAGS
        # put the flags on the ScanParam flag array
        #self.FlagHandler.setOnMask(mask,iFlags=flag)
        # put the flags on the main flag array
        #for chan in chanListIndices:
        #    self.FlagHandler.setOnMask(mask, self.rflags['INTEGRATION FLAGGED'], \
        #                               dim=1, index=chan)



    #----------------------------------------------------------------------------

    
    def findSubscanOld(self,threshold=1.):
        """
        DES: compute subscan indices from steps in az, el
        INP: (float) threshold = value (in arcsec^2) of (d_az^2 + d_el^2) step
                                 used to detect turnovers / stationnary points
        """
        # Retrieve all Az, El offsets (even if some are flagged) to
        # indices refering to the complete dataset
        azOff = self.get('AzimuthOffset',flag='None')  # in arcsec
        elOff = self.get('ElevationOffset',flag='None')

        # 1st derivative
        d_az = azOff[1::] - azOff[:-1]
        d_el = elOff[1::] - elOff[:-1]

        # compute squared radius
        radius = (d_az*d_az) + (d_el*d_el)
        # radius go  below 10. at turnovers and stationary positions
        # TODO: the latter could be used to estimate noise
        # compute difference with threshold to look at its sign
        radius_1 = radius - threshold
        diff_radius_1 = radius_1[1::] * radius_1[:-1]
        # this is negative when the sign of radius_1 has changed
        mask = nonzero(where(less(diff_radius_1,0),1,0))
    
        # Reinitialise Subscan related attributes
        self.SubscanIndex = []
        self.SubscanNum   = []
        self.SubscanType  = []

        nb = len(mask)-1
        for i in range(nb):
          self.SubscanIndex.append([mask[i],mask[i+1]-1])  # -1: the gaps
                                  # between subscans are not in any subscan
          self.SubscanNum.append(i+1)
          self.SubscanType.append('ON')

        self.SubscanIndex = transpose(self.SubscanIndex)
        self.NObs = nb
        self.__MessHand.info(str("Found %i subscans"%(nb)))

    #----------------------------------------------------------------------------

    def plotSubscan(self):
        """
        DES: generate a plot showing starting and ending times of subscans
        """
        index = self.SubscanIndex
        nbSub = len(self.SubscanNum)
        mjd = self.get('mjd',flag='None')
        mini = mjd[index[0,0]]
        maxi = mjd[index[1,-1]-1]
        maxSub = max(self.SubscanNum)
        Plot.plot([mini,maxi],[0,0],limitsY=[-1,maxSub],style='l',\
                  labelX='MJD - MJD[0] (s)',labelY='Subscan number')
        for num in range(nbSub):
            Plot.plot([mjd[index[0,num]],mjd[index[1,num]-1]],\
                      [self.SubscanNum[num],self.SubscanNum[num]],\
                      style='l',overplot=1)


    def plotSubscanOffsets(self,overplot=0):
        """
        DES: Use four colours to show subscans on the Az, El pattern
        INP: (logical) overplot : if set, do not plot AzElOffset - assume
                                  these have been plotted already
        """
        if not overplot:
            self.plotAzElOffset()
        index = self.SubscanIndex
        nbSub = len(self.SubscanNum)
        azOff = 3600. * self.LonOff
        elOff = 3600. * self.LatOff

        for num in range(nbSub):
            Plot.plot(azOff[index[0,num]:index[1,num]-1],\
                      elOff[index[0,num]:index[1,num]-1],\
                      style='p',overplot=1,ci=num-4*int(num/4.)+2)


#----------------------------------------------------------------------------------
#----- BoA Data Entity Class ------------------------------------------------------
#----------------------------------------------------------------------------------
class DataEntity: 
    """
    NAM: DataEntity (class)
    DES: Objects of this class store the data and associated
         parameters of a scan, which can contain several observations
         (or subscans).
         They also contain additional arrays in which the current
         results of the data reduction are stored. 
         This class also provides the interface between the MB-FITS
         files and BoA, by the means of the fillFromMBFits() method.
    """

    def __init__(self):
        """
        DES: Instanciation of a new DataEntity object.
             All attributes are defined and set to default values.
        """

        # Add a MessHand attribute - new MessageHandler 20050303
        self.MessHand = BoaMessageHandler.MessHand(self.__module__)
   
        # Dictionary containing status informations
        self.Status_Dic = {'Gain_Ele_Cor_Done':0}
        self.Status_Dic['Baseline_Cor_Ord'] = 0
        self.Status_Dic['Noi_Cor_Done'] = 0
        self.Status_Dic['Opa_Cor_Done'] = 0
        self.Status_Dic['Flux_Cal_Done'] = 0
        self.Status_Dic['Pha_Dif_Done'] = 0
        # ...to be completed...

        # Initialisation of all other attributes
        # All are arrays of floats, except when stated
        self.Data        = array([],Float32)
        self.DataBackup  = array([],Float32)   # The backup copy of data
        self.CorrelatedNoise    = array([],Float32)   # The associated noise array
        self.CorMatrix   = array([],Float32)
        self.FFCF_CN     = array([],Float32)
        self.FFCF_Gain   = array([],Float32)
        self.DataWeights = array([],Float32)
        self.Weight      = array([],Float32)   # TODO: move to BolometerArray

        self.BolometerArray = BolometerArray() # Contains all the array parameters
        self.ScanParam      = ScanParameter()  # contains all coordinates and times
        
        # Statistics 
        self.ChanRms  =  array([],Float32)
        self.ChanMean =  array([],Float32)
        self.ChanMed  =  array([],Float32)
        self.ChanRms_s  =  array([],Float32)
        self.ChanMean_s =  array([],Float32)
        self.ChanMed_s  =  array([],Float32)

        # File name
        self.FileName = ''

        # Flags:
        dataFlags = array([],Int8)
        dataFlags.shape = (0,0)
        self.FlagHandler = BoaFlagHandler.createFlagHandler(dataFlags)

        self.dflags = {'SPIKE'              : 1,
                       'GLITCH TYPE 1'      : 2,
                       'GLITCH TYPE 2'      : 3,
                       'SQUID FLUX JUMP'    : 4,
                       'DEMOD RAILED'       : 5,
                       'TEMPORARY'          : 8 }
        # Reserved flag values:
        self.rflags = {'INTEGRATION FLAGGED': 6,
                       'CHANNEL FLAGGED'    : 7 }

    def reset(self):
        """
        DES: Reset all attributes - useful before reading a new file
        """
        # Keep the same BolometerArray and ScanParam object --
        # otherwise shortcut will be lost !
        BolometerArray = self.BolometerArray
        ScanParam = self.ScanParam
 
        BolometerArray.__init__()
        ScanParam.__init__()

        # Also keep Message handler to keep trace of max. weight
        MessHand = self.MessHand
        
        self.__init__()
        self.BolometerArray = BolometerArray
        self.ScanParam      = ScanParam
        self.MessHand       = MessHand
        
        gc.collect()
        
        
    def __str__(self):
        """
        DES: Defines a string which is shown when the print instruction is
             used. It contains the sizes and typecodes of all attributes.
        """

        if not self._existData():
            return "No data read in yet\n"
        
        out = self.ScanParam.__str__()+"\n"
        out += self.BolometerArray.__str__()+"\n"
        
        if BoaConfig.DEBUG > 2:
            out += "\n" + \
                   attrStr(self,['MessHand','BolometerArray','Map','ScanParam'])+ \
                   "\n"
            
        return out

    # ---------------------------------------------------------------------
    # Overload addition operator: used to combine two datasets
    # ---------------------------------------------------------------------
    def __add__(self,other):
        # TODO: check that it makes sense to co-add these two datasets
        # (e.g. same bolo array, same source if pointing...)
        result = copy.deepcopy(self)
        result._coadd(other)
        return result

    def _coadd(self,other):
        # Special case of addition, applied to 'self' rather than returning a result
        self.ScanParam._coadd(other.ScanParam)

        self.Data        = concatenate((self.Data,        other.Data))
        self.DataWeights = concatenate((self.DataWeights, other.DataWeights))
        self.DataBackup  = concatenate((self.DataBackup,  other.DataBackup))
        #self.CorrelatedNoise    = concatenate((self.CorrelatedNoise,    other.CorrelatedNoise))

        slfFlags  = concatenate((self.FlagHandler.getFlags(), other.FlagHandler.getFlags()))
        self.FlagHandler = BoaFlagHandler.createFlagHandler(slfFlags)

                
    # ---------------------------------------------------------------------
    # General Input/Output methods
    # ---------------------------------------------------------------------
    def read(self,inFile='',febe='',baseband=0,subscans=[],update=0,phase=0, \
             channelFlag=1, integrationFlag=9, \
             readHe=0,readAzEl0=0,readT=0,readWind=0,readBias=0,readPWV=0):
        """
        DES: fill a data entity object
        INP: (int/string)  inFile: scan number / path to the dataset to be read
             (int list) subscans : subscan numbers to read (default: all)
                (logical) update : if true, do not reset previous entity object
                     (int) phase : phase to be stored (default: phase diff)
            channelFlag (i list) : flag for not connected feeds (default: 1 'NOT CONNECTED')
        integrationFlag (i list) : flag for blanked integrations (default: 9 'BLANK DATA')
                (logical) readHe : do we read LABOCA He3 tempe? (def: no)
             (logical) readAzEl0 : do we read monitor Az, El(0)? (def: no)
             (logical)     readT : do we read T_amb from monitor? (def: no)
             (logical)  readWind : do we read wind speed, dir...? (def: no)
             (logical)  readBias : do we need ASZCa bias settings? (def: no)
             (logical)  readPWV  : do we read pwv? (def: no)
        OUT: (int)        status : 0 if reading ok, <> 0 if an error occured
                   (see BoaDataAnalyser.read for error codes description)
        """
        t0 = time.clock()

        if type(inFile) == type(1):
            inFile = str(inFile)
            
        status = 0
        if update:
            newData = DataEntity()
            status = newData.read(inFile,febe,baseband,subscans,update=0,phase=phase)
            self._coadd(newData)
            # free memory
            newData = 0
            gc.collect()
            
        else:
            BoaCommandHistory.tagHistory()
            
            if not inFile == '':
                BoaDir.setInFile(inFile)

            # Open dataset and create MBFitsReader:
            datasetWasOpen = 0
            try:
                dataset = BoaMBFits.importDataset(BoaConfig.inFile)
            except:
                self.MessHand.error(" could not open dataset %s"%(inFile))
                status = -1
                return status

            reader = BoaMBFitsReader.createReader(dataset)
                
            # Reset the currData if data already present
            if self._existData():
                self.reset()

            reader.openSubscan(subsnum=None)

            # Determine febe to be used:
            febesDataset = reader.read("Febes")
            if febe:
                if febe in febesDataset:
                    useFebe = febe
                else:
                    self.MessHand.error(" no data for Febe %s in dataset %s" \
                                        %(febe,inFile))
                    status = -2
                    return status
            else:
                if len(febesDataset)==1:
                    # Only one febe in dataset: ok
                    useFebe = febesDataset[0]
                elif type(febesDataset) == type('test'):
                    # only one string - IRAM case
                    useFebe = febesDataset
                elif len(febesDataset)==0:
                    # No febe in dataset: not ok
                    self.MessHand.error(" no Febe data in dataset %s" \
                                        %(inFile))
                    status = -2
                    return status
                else:
                    # More than one febe in dataset: not ok
                    print "found febes = ",febesDataset
                    self.MessHand.error(" must specify Febe for dataset %s" \
                                        %(inFile))
                    status = -2
                    return status

            # Determine baseband to be used:
            basebandsDataset = reader.read("UseBand", febe=useFebe)
            if baseband:
                if baseband in basebandsDataset:
                    useBaseband = baseband
                else:
                    self.MessHand.error(" no data for baseband %s in dataset %s" \
                                        %(baseband,inFile))
                    status = -3
                    return status
            else:
                if len(basebandsDataset)==1:
                    # Only one baseband in dataset: ok
                    useBaseband = basebandsDataset[0]
                elif len(basebandsDataset)==0:
                    # No baseband in dataset: not ok
                    self.MessHand.error(" no baseband data in dataset %s" \
                                        %(inFile))
                    status = -3
                    return status
                else:
                    # More than one baseband in dataset: not ok
                    self.MessHand.error(" must specify baseband for dataset %s" \
                                        %(inFile))
                    status = -3
                    return status

            # Determine subscans to be used:
            subscansDataset = reader.read("Subscans")
            if not subscansDataset:
                self.MessHand.error(" no subscan data in dataset %s" \
                                    %(inFile))
                status = -4
                return status
            if not subscans:
                useSubscans = subscansDataset
            else:
                # Check whether all specified subscans are present:
                for subscan in subscans:
                    if not subscan in subscansDataset:
                        self.MessHand.error(" no data for subscan %d in dataset %s" \
                                            %(subscan,inFile))
                        status = -4
                        return status
                useSubscans = subscans

            # fillFromMBFits: here is most of the work done
            # ---------------
            # Use 1st subscan to define bolometer array
            self.BolometerArray._BolometerArray__fillFromMBFits(reader=reader, \
                                                                febe=useFebe, \
                                                                baseband=useBaseband, \
                                                                subscan=useSubscans[0], \
                                                                flag=channelFlag)
            
            # Copy the Telescope object to ScanParam: we need the latitude there
            self.ScanParam.Telescope = self.BolometerArray.Telescope
            
            # Fill ScanParam related attributes
            self.ScanParam._ScanParameter__fillFromMBFits(reader=reader,
                                                          febe=useFebe,
                                                          baseband=useBaseband,
                                                          subscans=useSubscans,
                                                          flag=integrationFlag, \
                                                          readHe=readHe,readAzEl0=readAzEl0,
                                                          readT=readT,readWind=readWind,
                                                          readBias=readBias,readPWV=readPWV)
            
            self._DataEntity__fillFromMBFits(reader=reader, \
                                             febe=useFebe, \
                                             baseband=useBaseband, \
                                             subscans=useSubscans)

            # Now report channels flagged at read
            usedChannels = self.BolometerArray.UsedChannels
            arrayFlagHandler = self.BolometerArray.FlagHandler
            timeFlagHandler  = self.ScanParam.FlagHandler
            dataFlagHandler  = self.FlagHandler

            if timeFlagHandler.nSet()>0:
                timeMask = timeFlagHandler.isSetMask()
            else:
                timeMask= None
            for chan in usedChannels:
                index = self.BolometerArray.getChanIndex(chan)[0]
                if arrayFlagHandler.isSetOnIndex(chan-1):
                    dataFlagHandler.setAll(self.rflags['CHANNEL FLAGGED'], \
                                           dim=1, index=index)
                if timeMask:
                    dataFlagHandler.setOnMask(timeMask, \
                                              self.rflags['INTEGRATION FLAGGED'], \
                                              dim=1, index=index)

            # now choose phase to be stored in Boa
            # Compute integration numbers for both phases
#             if self.ScanParam.WobUsed:
#                 self.ScanParam.computeOnOff()
#             if phase == 0:
#                 # compute phase diff, if needed
#                 if self.ScanParam.OnOffPairs:
#                     self._phaseDiff()
#             else:
#                 # return appropriate phase, if exists
#                 if not (self.ScanParam.OnOffPairs):
#                     self.MessHand.warning("Only one phase in these data -" +\
#                                           " returning complete dataset")
#                 else:
#                     self.selectPhase(phase)

            dataset.close()

            t1 = time.clock()
            self.MessHand.debug(str(len(subscans))+" subscan(s) read  "+str(t1-t0))
            t0 = t1

            # store file name in the DataEntity object
            self.FileName = inFile
            # Display some general infos about the file
            self.MessHand.info(self.__str__())

            # At read time make the current channel selection equal the total channel list
            self.BolometerArray.setCurrChanList('all')

            self.MessHand.debug(" Dataset "+inFile+" has been read")
            
            gc.collect()
            return status
                    
    # ---------------------------------------------------------------------
    def __fillFromMBFits(self,reader,febe,baseband,subscans):
        """
        NAM: fillFromMBFits()
        DES: fill a DataEntity object using the MBFitsReader object reader.
        INP:            reader : MBFitsReader object
                baseband (int) : baseband number to select
             subscans (i list) : list of subscans numbers to read in
        """

        self.MessHand.debug('start of fillFromMBfits')
        t0 = time.clock()

        # Now get the data
        subscan = subscans[0]
        nInt = self.ScanParam.NInt
        nUseFeed = reader.read("NUseFeed", \
                               subsnum=subscan, \
                               febe=febe, \
                               baseband=baseband)
        subsIndex = 0
        Data = zeros((nInt,nUseFeed),Float32)
        for subscan in subscans:
            if subscan in self.ScanParam.SubscanNum:
                # this means ScanParam.__fillFromMBFits worked fine
                subscanWasOpened = reader.openSubscan(subsnum=subscan)
                subscanStart = self.ScanParam.SubscanIndex[0,subsIndex]
                subscanStop  = self.ScanParam.SubscanIndex[1,subsIndex]
                subsIndex += 1

                tmpData = reader.read("Data",
                                      subsnum=subscan,
                                      febe=febe,
                                      baseband=baseband)
                if shape(tmpData)[0] > subscanStop - subscanStart:
                    # Happens when more rows in ARRAYDATA than DATAPAR
                    tmpData = tmpData[:subscanStop - subscanStart,::]
                Data[subscanStart:subscanStop] = tmpData.astype(Float32)
                                
                if subscanWasOpened:
                    reader.closeSubscan(subsnum=subscan)

        # Add the DC offsets to the signals
        #for i in range(self.BolometerArray.NUsedChannels):
            # Note: for LABOCA-ABBA the ordering of DCOff may be wrong
            # (see Dirk's e-mail 2006/09/09)
            #Data[:,i] = Data[:,i] + array(self.BolometerArray.DCOff[i]).astype(Float32)
            # Commented out 2006/9/17 - no meaningful unit in LABOCA data

        # Apply gain factor - Note: to convert to Volts, one would have to divide
        # by 32768. in the specific case of LABOCA-ABBA - We need to know the
        # ADC dynamic range, not in MB-FITS format definition
        Data = Data / array(self.BolometerArray.BEGain).astype(Float32)
        # Also correct for frontend gain
        Data = Data / (array(2.**self.BolometerArray.FEGain).astype(Float32))
        self.Data = fUtilities.as_column_major_storage(Data)
        
        # Initialise Data Weights
        self.DataWeights = fUtilities.as_column_major_storage(ones((nInt,nUseFeed),Float32))
        
        # The other arrays can be initalised to appropriate sizes
        # Extract the numbers of channels (= pixels) and datapoints
        dataShape = shape(self.Data)
        self.FlagHandler = BoaFlagHandler.createFlagHandler(zeros(dataShape,Int8))

        t1 = time.clock()
        self.MessHand.debug(" Complementary information filled "+str(t1-t0)) 
        self.MessHand.debug('end of fillFromMBfits')
            
    # -------------------------------------------------------------------
    def dumpData(self,fileName='BoaData.sav'):
        """
        DES: save the current DataEntity object to a file
        INP: (string) fileName: name of the output file
             optional - default value = 'BoaData.sav'
        """
        #fileName = self.outDir+fileName
        try:
            f = file(os.path.join(BoaConfig.outDir,fileName),'w')
        except IOError:
            self.MessHand.error(" permission denied, please change outdir")
            return
        cPickle.dump(self,f,2)
        f.close()
        self.MessHand.info(" current data successfully written to %s"%fileName)
    # -------------------------------------------------------------------
    def restoreData(self,fileName='BoaData.sav'):
        """ 
        DES: restore a DataEntity object previously saved in a file, and
             set it as the currData attribute of BoaB
        INP: (string) fileName: name of the input file
             optional - default value = 'BoaData.sav'
        """
        #fileName = self.outDir+fileName
        try:
            f = file(os.path.join(BoaConfig.outDir,fileName))
        except IOError:
            self.MessHand.error(" could not open file %s"%(fileName))
            return
        self = cPickle.load(f)
        f.close()
        
    # -------------------------------------------------------------------
    def backup(self):
        """ 
        DES: backup the data
        """
        self.DataBackup = copy.copy(self.Data)

    # -------------------------------------------------------------------
    def restore(self):
        """ 
        DES: backup the data
        """
        self.Data = copy.copy(self.DataBackup)
        
    # -------------------------------------------------------------------
    def saveMambo(self,inName='',outName=''):
        """
        DES: convert an MB-Fits file to the MAMBO FITS format, readable
             by MOPSIC
        INP: (str) inName: name of the MB-Fits file (optional)
             (str) outName: name of the MAMBO output file (optional)

        """
        # if input name is not given, use that of the current file
        if inName == '':
            inName = self.FileName
        # if output name is not given, simply use the scan number
        if outName == '':
            #outName = self.outDir+str(self.scanNum)
            outName = str(self.ScanParam.ScanNum)
        # Now use the MamboMBFits.py module
        m = MamboMBFits.MamboMBFits(outName,inName)
        try:
            m.convertMB2MamboFits()
            self.MessHand.setMess(3," MB-FITS file "+inName)
            self.MessHand.setMess(3," successfully converted to "+outName)
        except:
            self.MessHand.setMess(1," something got wrong in MB to MAMBO conversion")
            raise

    # -------------------------------------------------------------------
    def saveExchange(self, fileName="", overwrite=0):
        """
        DES: save information from the DataEntity object to a
             Fits file for exchange with other reduction packages
        INP: (str) fileName: name of the Fits file (optional)
             (log) overwrite: Overwrite existing file (optional)
        """
        # if file name is not given, simply use the scan number
        if fileName == '':
            fileName = "%s.fits" % str(self.ScanParam.ScanNum)
        fileName = os.path.join(BoaConfig.outDir,fileName)

        if overwrite:
            if os.path.exists(fileName):
                os.remove(fileName)
        
        # Keywords for the PrimaryHeader:
        keywords = []
        keywords.append(BoaMBFits.Keyword(name="CREATOR", \
                                          value="BOA", \
                                          datatype="10A", \
                                          comment="ID of program that created this file"))
        keywords.append(BoaMBFits.Keyword(name="TIME", \
                                          value=time.strftime("%Y-%m-%d:%H:%M:%S"), \
                                          datatype="20A", \
                                          comment="Time of creation"))
        keywords.append(BoaMBFits.Keyword(name="MBFITS", \
                                          value=self.FileName, \
                                          datatype="30A", \
                                          comment="Name of original MBFits dataset"))
        keywords.append(BoaMBFits.Keyword(name="SCAN", \
                                          value=str(self.ScanParam.ScanNum), \
                                          datatype="10A", \
                                          comment="Scan number"))
        keywords.append(BoaMBFits.Keyword(name="OBJECT", \
                                          value=self.ScanParam.Object, \
                                          datatype="30A", \
                                          comment="Object observed"))
        keywords.append(BoaMBFits.Keyword(name="FEBE", \
                                          value=self.BolometerArray.FeBe, \
                                          datatype="17A", \
                                          comment="Frontend-backend ID"))
        
        ds = BoaMBFits.createDataset(fileName, fileName, keywords, "")

        # DATA Table:
        keywords = []
        keywords.append(BoaMBFits.Keyword(name="EXTNAME", \
                                          value="DATA", \
                                          datatype="20A"))
        
        colinfos = []
        colinfos.append(BoaMBFits.ColumnInfo(name="MJD", \
                                             datatype="D", \
                                             description="Mjd"))
        colinfos.append(BoaMBFits.ColumnInfo(name="DATA", \
                                             datatype="E", \
                                             repeat=self.Data.shape[1], \
                                             description="Processed signal"))
        tData = ds.addTable(keywords=keywords, colinfos=colinfos)

        colMJD = tData.getColumn("MJD")
        colData = tData.getColumn("DATA")

        colMJD.write(1, self.ScanParam.MJD)
        colData.write(1, self.Data)

        # FLAGS Table:
        keywords = []
        keywords.append(BoaMBFits.Keyword(name="EXTNAME", \
                                          value="FLAGS", \
                                          datatype="20A"))
        keywords.append(BoaMBFits.Keyword(name="DFLAG1", \
                                          value="SPIKE", \
                                          datatype="20A", \
                                          comment="Description of flag"))
        keywords.append(BoaMBFits.Keyword(name="DFLAG2", \
                                          value="GLITCH TYPE 1", \
                                          datatype="20A", \
                                          comment="Description of flag"))
        keywords.append(BoaMBFits.Keyword(name="DFLAG3", \
                                          value="GLITCH TYPE 2", \
                                          datatype="20A", \
                                          comment="Description of flag"))
        keywords.append(BoaMBFits.Keyword(name="DFLAG4", \
                                          value="SQUID FLUX JUMP", \
                                          datatype="20A", \
                                          comment="Description of flag"))
        keywords.append(BoaMBFits.Keyword(name="DFLAG5", \
                                          value="DEMOD RAILED", \
                                          datatype="20A", \
                                          comment="Description of flag"))
        keywords.append(BoaMBFits.Keyword(name="DFLAG6", \
                                          value="", \
                                          datatype="20A", \
                                          comment="Description of flag"))
        keywords.append(BoaMBFits.Keyword(name="DFLAG7", \
                                          value="", \
                                          datatype="20A", \
                                          comment="Description of flag"))
        keywords.append(BoaMBFits.Keyword(name="DFLAG8", \
                                          value="", \
                                          datatype="20A", \
                                          comment="Description of flag"))

        keywords.append(BoaMBFits.Keyword(name="IFLAG1", \
                                          value="TURNAROUND", \
                                          datatype="20A", \
                                          comment="Description of flag"))
        keywords.append(BoaMBFits.Keyword(name="IFLAG2", \
                                          value="ACCELERATION THRESHOLD", \
                                          datatype="20A", \
                                          comment="Description of flag"))
        keywords.append(BoaMBFits.Keyword(name="IFLAG3", \
                                          value="ELEVATION VELOCITY THRESHOLD", \
                                          datatype="20A", \
                                          comment="Description of flag"))
        keywords.append(BoaMBFits.Keyword(name="IFLAG4", \
                                          value="", \
                                          datatype="20A", \
                                          comment="Description of flag"))
        keywords.append(BoaMBFits.Keyword(name="IFLAG5", \
                                          value="", \
                                          datatype="20A", \
                                          comment="Description of flag"))
        keywords.append(BoaMBFits.Keyword(name="IFLAG6", \
                                          value="", \
                                          datatype="20A", \
                                          comment="Description of flag"))
        keywords.append(BoaMBFits.Keyword(name="IFLAG7", \
                                          value="SUBSCAN FLAGGED", \
                                          datatype="20A", \
                                          comment="Description of flag"))
        keywords.append(BoaMBFits.Keyword(name="IFLAG8", \
                                          value="", \
                                          datatype="20A", \
                                          comment="Description of flag"))
        for iflag in xrange(9,33):
            keywords.append(BoaMBFits.Keyword(name="IFLAG%d"%iflag, \
                                              value="", \
                                              datatype="20A", \
                                          comment="Description of flag"))

        
        colinfos = []
        colinfos.append(BoaMBFits.ColumnInfo(name="MJD", \
                                             datatype="D", \
                                             description="Mjd"))
        colinfos.append(BoaMBFits.ColumnInfo(name="DATAFLAG", \
                                             datatype="B", \
                                             repeat=self.FlagHandler.getFlags().shape[1], \
                                             description="Flag per channel and Mjd"))
        colinfos.append(BoaMBFits.ColumnInfo(name="INTEGFLAG", \
                                             datatype="J", \
                                             description="Flag per Mjd"))
        tFlags = ds.addTable(keywords=keywords, colinfos=colinfos)

        colMJD = tFlags.getColumn("MJD")
        colDataflag = tFlags.getColumn("DATAFLAG")
        colIntegflag = tFlags.getColumn("INTEGFLAG")

        colMJD.write(1, self.ScanParam.MJD)

        # Clear the Boa specific flags from DataFlags:
        dFlags = copy.copy(self.FlagHandler.getFlags())
        flagHandler = BoaFlagHandler.createFlagHandler(dFlags)
        
        flagHandler.unsetAll([self.rflags['INTEGRATION FLAGGED'], \
                              self.rflags['CHANNEL FLAGGED'], \
                              self.dflags['TEMPORARY'] ])
        colDataflag.write(1, flagHandler.getFlags())
        colIntegflag.write(1, self.ScanParam.FlagHandler.getFlags())


        # CHANNELFLAGS Table:
        keywords = []
        keywords.append(BoaMBFits.Keyword(name="EXTNAME", \
                                          value="CHANNELFLAGS", \
                                          datatype="20A"))

        keywords.append(BoaMBFits.Keyword(name="CFLAG1", \
                                          value="NOT CONNECTED", \
                                          datatype="20A", \
                                          comment="Description of flag"))
        keywords.append(BoaMBFits.Keyword(name="CFLAG2", \
                                          value="BAD SENSITIVITY", \
                                          datatype="20A", \
                                          comment="Description of flag"))
        keywords.append(BoaMBFits.Keyword(name="CFLAG3", \
                                          value="LOW SENSITIVITY", \
                                          datatype="20A", \
                                          comment="Description of flag"))
        keywords.append(BoaMBFits.Keyword(name="CFLAG4", \
                                          value="DARK BOLOMETER", \
                                          datatype="20A", \
                                          comment="Description of flag"))
        for iflag in xrange(5,33):
            keywords.append(BoaMBFits.Keyword(name="CFLAG%d"%iflag, \
                                              value="", \
                                              datatype="20A", \
                                          comment="Description of flag"))

        
        colinfos = []
        colinfos.append(BoaMBFits.ColumnInfo(name="CHANNEL", \
                                             datatype="J", \
                                             description="Channel number"))
        colinfos.append(BoaMBFits.ColumnInfo(name="FLAG", \
                                             datatype="J", \
                                             description="Flag per channel"))
        tCflags = ds.addTable(keywords=keywords, colinfos=colinfos)

        colChannel = tCflags.getColumn("CHANNEL")
        colFlag = tCflags.getColumn("FLAG")

        colChannel.write(1, arrayrange(self.BolometerArray.NChannels)+1)

        colFlag.write(1, self.BolometerArray.FlagHandler.getFlags())


        # SUBSCANS Table:
        keywords = []
        keywords.append(BoaMBFits.Keyword(name="EXTNAME", \
                                          value="SUBSCANS", \
                                          datatype="20A"))
        
        colinfos = []
        colinfos.append(BoaMBFits.ColumnInfo(name="SUBSNUM", \
                                             datatype="J", \
                                             description="Subscan number"))
        colinfos.append(BoaMBFits.ColumnInfo(name="ISTART", \
                                             datatype="J", \
                                             description="Start index of subscan"))
        colinfos.append(BoaMBFits.ColumnInfo(name="IEND", \
                                             datatype="J", \
                                             description="Stop index of subscan"))
        colinfos.append(BoaMBFits.ColumnInfo(name="TYPE", \
                                             datatype="A", \
                                             repeat=10, \
                                             description="Subscan type"))
        colinfos.append(BoaMBFits.ColumnInfo(name="TIME", \
                                             datatype="E", \
                                             description="Subscan duration"))


        tSubs = ds.addTable(keywords=keywords, colinfos=colinfos)

        colNum = tSubs.getColumn("SUBSNUM")
        colStart = tSubs.getColumn("ISTART")
        colEnd = tSubs.getColumn("IEND")
        colType = tSubs.getColumn("TYPE")
        colTime = tSubs.getColumn("TIME")

        colNum.write(1, self.ScanParam.SubscanNum)
        colStart.write(1, self.ScanParam.SubscanIndex[0])
        colEnd.write(1, self.ScanParam.SubscanIndex[1]-1)
        colType.write(1, self.ScanParam.SubscanType)
        colTime.write(1, self.ScanParam.SubscanTime)


        # HISTORY Table:
        keywords = []
        keywords.append(BoaMBFits.Keyword(name="EXTNAME", \
                                          value="HISTORY", \
                                          datatype="20A"))
        
        colinfos = []
        colinfos.append(BoaMBFits.ColumnInfo(name="PROGID", \
                                             datatype="A", \
                                             repeat=10, \
                                             description="ID of program"))
        colinfos.append(BoaMBFits.ColumnInfo(name="COMMAND", \
                                             datatype="A", \
                                             repeat=100, \
                                             description="Command"))
        tHist = ds.addTable(keywords=keywords, colinfos=colinfos)

        commands = BoaCommandHistory.getHistory(-1)
        progIDs = ["BOA"]*len(commands)
        if commands:
            colProgID = tHist.getColumn("PROGID")
            colCommand = tHist.getColumn("COMMAND")

            colProgID.write(1, progIDs)
            colCommand.write(1, commands)
        
        

        ds.close()

    # -------------------------------------------------------------------
    def loadExchange(self, fileName=""):
        """
        DES: read information from a Fits file for exchange with other
             reduction packages into the DataEntity object
        INP: (str) fileName: name of the Fits file
        """
        fileName = os.path.join(BoaConfig.outDir,fileName)

        try:
            ds = BoaMBFits.importDataset(fileName)
        except:
            self.MessHand.setMess(1," could not load file %s" % fileName)
            raise

        # Check information in Primary Header:
        mbfitsFile = ds.getKeyword("MBFITS").getValue()
        if mbfitsFile != self.FileName:
            self.MessHand.setMess(1," MBFitsFile in file %s does not match" % fileName)
            raise

        febe = ds.getKeyword("FEBE").getValue()
        if febe != self.BolometerArray.FeBe:
            self.MessHand.setMess(1," FeBe in file %s does not match" % fileName)
            raise

        # DATA Table:
        tData = ds.getTables(EXTNAME="DATA")[0]
        tData.open()
        
        colMJD = tData.getColumn("MJD")
        colData = tData.getColumn("DATA")

        self.ScanParam.MJD = colMJD.read().astype(Float)
        self.Data = colData.read().astype(Float32)

        tData.close()

        # FLAGS Table:
        tFlags = ds.getTables(EXTNAME="FLAGS")[0]
        tFlags.open()

        colDataflag = tFlags.getColumn("DATAFLAG")
        colIntegflag = tFlags.getColumn("INTEGFLAG")

        dataFlags = colDataflag.read().astype(Int8)
        timeFlags = colIntegflag.read().astype(Int32)
        
        tFlags.close()

        # CHANNELFLAGS Table:
        tCflags = ds.getTables(EXTNAME="CHANNELFLAGS")[0]
        tCflags.open()

        colFlag = tCflags.getColumn("FLAG")

        arrayFlags = colFlag.read().astype(Int32)
        
        tCflags.close()

        # Adjust the flags:
        usedChannels = self.BolometerArray.UsedChannels
        dataFlagHandler = BoaFlagHandler.createFlagHandler(dataFlags)
        timeFlagHandler = BoaFlagHandler.createFlagHandler(timeFlags)
        arrayFlagHandler = BoaFlagHandler.createFlagHandler(arrayFlags)
        
        if timeFlagHandler.nSet()>0:
            timeMask = timeFlagHandler.isSetMask()
        else:
            timeMask= None
        for chan in usedChannels:
            index = self.BolometerArray.getChanIndex(chan)[0]
            if arrayFlagHandler.isSetOnIndex(chan-1):
                dataFlagHandler.setAll(self.rflags['CHANNEL FLAGGED'], dim=1, index=index)
            if timeMask:
                dataFlagHandler.setOnMask(timeMask, \
                                          self.rflags['INTEGRATION FLAGGED'], \
                                          dim=1, index=index)

        self.FlagHandler = dataFlagHandler
        self.ScanParam.FlagHandler = timeFlagHandler
        self.BolometerArray.FlagHandler = arrayFlagHandler


        # SUBSCANS Table:
        tSubs = ds.getTables(EXTNAME="SUBSCANS")[0]
        tSubs.open()

        colNum = tSubs.getColumn("SUBSNUM")
        colStart = tSubs.getColumn("ISTART")
        colEnd = tSubs.getColumn("IEND")
        colType = tSubs.getColumn("TYPE")
        colTime = tSubs.getColumn("TIME")

        self.ScanParam.SubscanNum = colNum.read().tolist()
        self.ScanParam.NObs = len(self.ScanParam.SubscanNum)
        iStarts = colStart.read().astype(Int)
        iEnds = colEnd.read().astype(Int) + 1
        self.ScanParam.SubscanIndex = array([iStarts, iEnds], Int)
        self.ScanParam.SubscanType = colType.read()
        self.ScanParam.SubscanTime = colTime.read().tolist()
        
        tSubs.close()


        ds.close()


    # ---------------------------------------------------------------------
    # Methods for computing phase differences
    # ---------------------------------------------------------------------

    def _phaseDiff(self):
        """
        NAM: phaseDiff (method)
        DES: Compute phase differences: call ScanParam.phaseDiffParam for
             coordinates and times, and compute Data(ON) - Data(OFF)
        """

        # phase diff on positions and times
        t0 = time.clock()
        self.MessHand.debug("Entering ScanParam.phaseDiff")
        self.ScanParam._phaseDiffParam()
        t1 = time.clock()
        self.MessHand.debug("end of ScanParam.phaseDiff() "+str(t1-t0)) 
                
        # phase diff on data: ON - OFF
        self.Data = take(self.Data,self.ScanParam.OnOffPairs[:,0]) - \
                    take(self.Data,self.ScanParam.OnOffPairs[:,1])
        # Do not initialize backup on phasediff
        #self.DataBackup = copy.deepcopy(self.Data)
        # CorrelatedNoise also not initialised
        #self.CorrelatedNoise = take(self.CorrelatedNoise,self.ScanParam.OnOffPairs[:,0]) - \
        #                        take(self.CorrelatedNoise,self.ScanParam.OnOffPairs[:,1])
        # Flags: use binary_or of ON and OFF flags
        slfFlags = self.FlagHandler.getFlags()
        slfFlags = take(slfFlags,self.ScanParam.OnOffPairs[:,0]) | \
                   take(slfFlags,self.ScanParam.OnOffPairs[:,1])
        self.FlagHandler = BoaFlagHandler.createFlagHandler(slfFlags)

        # Weights: use average of ON and OFF
        self.DataWeights = (take(self.DataWeights,self.ScanParam.OnOffPairs[:,0]) + \
                                take(self.DataWeights,self.ScanParam.OnOffPairs[:,1]))/2.
       
        self.ScanParam._ScanParameter__computeSubIndex()

    def selectPhase(self,phase):
        """
        NAM: selectPhase (method)
        DES: Keep only Data(ON) or Data(OFF)
        INP: (int) phase: phase to keep, 1=ON, 2=OFF
        """

        #BROKEN : because of Subscan Index calculation

        # Select times and positions
        self.ScanParam.selectPhase(phase)
        
        # select data
        ph = phase-1  # index in OnOffPairs: 0 = ON, 1 = OFF
        self.Data       = take(self.Data,self.ScanParam.OnOffPairs[:,ph])
        # Do not initialize backup on phase
        #self.DataBackup = copy.deepcopy(self.Data)
        self.CorrelatedNoise    = take(self.CorrelatedNoise,self.ScanParam.OnOffPairs[:,ph])
        self.DataWeights = take(self.DataWeights,self.ScanParam.OnOffPairs[:,ph])

        slfFlags = self.FlagHandler.getFlags()
        slfFlags = take(slfFlags,self.ScanParam.OnOffPairs[:,ph])
        self.FlagHandler = BoaFlagHandler.createFlagHandler(slfFlags)

        self.ScanParam._ScanParameter__computeSubIndex()

    #----------------------------------------------------------------------------
    def _existData(self):
        """
        DES: check if the DataEntity object has been filled with data
        OUT: (int) result: 0 if no data, 1 otherwise
        """
        if (len(self.Data) > 0):
            return 1
        else:
            return 0

    #----------------------------------------------------------------------------
    #----- Data extraction ------------------------------------------------------
    #----------------------------------------------------------------------------
    def getChanData(self,dataType=' ',chan='None', flag=[], getFlagged=0,
                    flag2=None, subscans=[]):
        """
        DES: get data for one channel 
        INP: (string)   dataType : type of data
             (int)          chan : channel number 
             (integer list) flag : retrieve data flagged or unflagged accordingly
             (log)    getFlagged : flag revers to flagged/unflagged data
                                   flag   | getFlagged | Retrieve..
                                   'None' |  0         | all data
                                   []     |  0         | unflagged data (default)
                                   []     |  1         | data with at least one flag set
                                   1      |  0         | data with flag 1 not set
                                   1      |  1         | data with flag 1 set
                                   [1,2]  |  0         | data with neither flag 1 nor flag 2 set
                                   [1,2]  |  1         | data with either flag 1 or flag 2 set
             (int list) subscans : list of wanted subscan (default all)
        OPT: (int array)   flag2 : second array of flags to check 
        OUT: (float)       array : data of one channel
        """

        self.MessHand.debug("getChanData start...")
        # if the channel is not given take the reference channel as default
        if chan in ['None','none','']:
            chan = self.BolometerArray.RefChannel

        # retrieve the corresponding channel index
        chanIndex = self.BolometerArray.getChanIndex(chan)[0]

        if chanIndex == -1:
            self.MessHand.error("channel not used") 
            return

        # retrieve data
        if dataType in ['FLUX', 'Flux', 'flux']:
            dataArray = self.Data[:,chanIndex]

        elif dataType in ['Weight','weight']:
            dataArray = self.DataWeights[:,chanIndex]
            
        elif dataType in ['LST','lst',\
                          'MJD','mjd',\
                          'UT']:
            dataArray = self.ScanParam.get(dataType=dataType,flag='None')
        elif dataType in ['FOCUS-X','FocusX','FocX','focX','focx',\
                          'FOCUS-Y','FocusY','FocY','focY','focy',\
                          'FOCUS-Z','FocusZ','FocZ','focZ','focz',\
                          'FOCUS-XTILT','Focus-XTILT','FOCUS-YTILT','Focus-YTILT']:
            dataArray = self.ScanParam.get(dataType=dataType,flag='None')
        elif dataType in ['AzimuthOffset', 'azimuthoffset', 'AzOff', 'azoff','azo']:
            dataArray = self.ScanParam.get(dataType=dataType,flag='None')
            # add the channel offset in arcsec
            dataArray = dataArray + self.BolometerArray.getChanSep([chan])[0]
        elif dataType in ['AzimuthOffsetRef']:
            dataArray = self.ScanParam.get(dataType='AzimuthOffset',flag='None')
            # do not add the channel offset!
        elif dataType in ['Azimuth','azimuth','Az','az']:
            dataArray = self.ScanParam.get(dataType=dataType,flag='None')
            # add the channel offset in deg.
            dataArray = dataArray + self.BolometerArray.getChanSep([chan])[0] / 3600.
        elif dataType in ['ElevationOffset', 'elevationoffset', 'ElOff','eloff','elo']:
            dataArray = self.ScanParam.get(dataType=dataType,flag='None')
            # add the channel offset in arcsec
            dataArray = dataArray + self.BolometerArray.getChanSep([chan])[1]
        elif dataType in ['Elevation','elevation','El','el']:
            dataArray = self.ScanParam.get(dataType=dataType,flag='None')
            # add the channel offset in deg.
            dataArray = dataArray + self.BolometerArray.getChanSep([chan])[1] / 3600.
        elif dataType in ['azspeed','Azspeed','AzSpeed','elspeed','Elspeed','ElSpeed']:
            dataArray = self.ScanParam.get(dataType=dataType,flag='None')
        elif dataType in ['acc','Acc','accel','Accel']:
            dataArray = self.ScanParam.get(dataType=dataType,flag='None')
        elif dataType in ['azacc','Azacc','AzAcc','elacc','Elacc','ElAcc']:
            dataArray = self.ScanParam.get(dataType=dataType,flag='None')
        elif dataType in ['speed','Speed']:
            dataArray = self.ScanParam.get(dataType=dataType,flag='None')
        elif dataType in ['ra','Ra','RA','dec','Dec','DEC']:
            dataArray = self.ScanParam.get(dataType=dataType,flag='None')
        elif dataType in ['raoff','RaOff','RAOff','raOffset','RaOffset','RAOffset']:
            dataArray = self.ScanParam.get(dataType=dataType,flag='None')
        elif dataType in ['decoff','DecOff','DECOff','decOffset','DecOffset','DECOffset']:
            dataArray = self.ScanParam.get(dataType=dataType,flag='None')
        elif dataType in ['HeT','HeTemp']:
            dataArray = self.ScanParam.get(dataType=dataType,flag='None')
        elif dataType in ['phase','Phase','wob','wobbler']:
            dataArray = self.ScanParam.get(dataType=dataType,flag='None')
            
        elif dataType in ['Flags', 'flags', 'flag']:
            dataArray = self.FlagHandler.getFlags()[:,chanIndex]
            
        elif dataType in ['SkyNoise', 'skynoise', 'SN', 'sn','CorrelatedNoise','correlatedNoise', 'correlatednoise', 'CN', 'cn']:
            dataArray = self.CorrelatedNoise[:,chanIndex]

        # The following types are computed per scan (and not per integration)
        
        elif dataType in ['Mean','mean']:
            dataArray = self.ChanMean[chanIndex]
            return dataArray

        elif dataType in ['RMS','Rms','rms']:
            dataArray = self.ChanRms[chanIndex]
            return dataArray

        # The following types are computed per subscan (and not per integration)
        # and make sense only on phase diff'ed data
        # WARNING: for these cases, no flag filtering is done, because we only
        # have flags per integration - NO SUBSCAN FLAGGING YET
        # Therefore, return the array before going to second step, which would crash

        elif dataType in ['Mean_s','mean_s']:
            dataArray = self.ChanMean_s[chanIndex]
            return dataArray
        
        elif dataType in ['RMS_s','Rms_s','rms_s']:
            dataArray = self.ChanRms_s[chanIndex]
            return dataArray
        
        elif dataType in ['Subscan','subscan']:
            dataArray = self.ScanParam.SubscanNum
            return dataArray
        
        else:
            self.MessHand.error("Unknown data")
            return

        dataFlags = self.FlagHandler.getFlags()[:,chanIndex].astype(Int32)
        # Check if subscans was asked
        if subscans:
            # Create a mask representing the asked subscans
            mask = zeros(shape(dataArray))
            
            SubscanNum   = self.ScanParam.SubscanNum
            SubscanIndex = self.ScanParam.SubscanIndex
            
            for subscan in subscans:
                if subscan in SubscanNum:
                    isub = SubscanNum.index(subscan)
                    mask[SubscanIndex[0,isub]:SubscanIndex[1,isub]] = 1
                else:
                    self.MessHand.error("subscan "+subscan+" does not exist")
                    return
                
            # Extract that mask from the data
            dataArray = compress(mask, dataArray)
            dataFlags = compress(mask, dataFlags)
            
            
        # .. and only return the desired flag
        dataFlagHandler = BoaFlagHandler.createFlagHandler(dataFlags)
        if flag in ['', 'None']:
            return dataArray
        elif flag in ['Blank','blank']:
            mask = dataFlagHandler.isSetMask()
            return where(mask, float('Nan'), dataArray)
        else:
            if getFlagged:
                mask = dataFlagHandler.isSetMask(flag)
                if flag2:
                    flagHandler2 = BoaFlagHandler.createFlagHandler(flag2)
                    mask2 = flagHandler2.isSetMask(flag)
                    bitwise_or(mask, mask2, mask)
            else:
                mask = dataFlagHandler.isUnsetMask(flag)
                if flag2:
                    flagHandler2 = BoaFlagHandler.createFlagHandler(flag2)
                    mask2 = flagHandler2.isUnsetMask(flag)
                    mask = bitwise_and(mask, mask2)
            return compress(mask, dataArray)
        
        self.MessHand.debug("... end getChanData")

    #----------------------------------------------------------------------------
    def getChanListData(self,type=' ',chanList=[],\
                        channelFlag=[], getFlaggedChannels=0, \
                        dataFlag=[], getFlaggedData=0, dataFlag2=None, \
                        subscans=[]):
        """
        DES: get data for list of channels  
        INP: (string) type = type of data
             (int list) chan = channel list 
             (integer list) channelFlag : retrieve data from channels flagged or unflagged accordingly
             (log)   getFlaggedChannels : channelFlag revers to flagged/unflagged data
             (integer list)    dataFlag : retrieve data flagged or unflagged accordingly
             (log)       getFlaggedData : dataFlag revers to flagged/unflagged data
                                   flag   | getFlagged | Retrieve..
                                   'None' |  0         | all data
                                   []     |  0         | unflagged data (default)
                                   []     |  1         | data with at least one flag set
                                   1      |  0         | data with flag 1 not set
                                   1      |  1         | data with flag 1 set
                                   [1,2]  |  0         | data with neither flag 1 nor flag 2 set
                                   [1,2]  |  1         | data with either flag 1 or flag 2 set
             (int array) dataFlag2 = second array of flags to check (optional)
        OUT: (list of float arrays) = data of the input list of channels
        """

        t0=time.clock()

        if getFlaggedChannels:
            dataFlag = self._removeReservedFlagValues(dataFlag, self.rflags['CHANNEL FLAGGED'])
            if dataFlag==None:
                self.MessHand.error("no valid flags")
                return
        
        chanList = self.BolometerArray.checkChanList(chanList, \
                                                     flag=channelFlag,getFlagged=getFlaggedChannels)
        result=[]

        for i in chanList:
            # Get data chan by chan
            result.append(self.getChanData(type,i,\
                                           flag=dataFlag,getFlagged=getFlaggedData,flag2=dataFlag2,\
                                           subscans=subscans))

        t1=time.clock()
        self.MessHand.debug(" cutting data for a list of channels: "+\
                                   str(t1-t0))
        return result

    #--------------------------------------------------------------------------------

    def _removeReservedFlagValues(self, flag, removeFlags=[]):
        """
        DES: Removes the reserved flag values defined in self.rflag from list of flag values
        INP: (i list) flag: flag values (integers, strings, and [] allowed)
             (i list) removeFlags: flag value to be removed.
        OUT: (i list)     : flag values with reserved values removed.
                            If the resulting is empty, None is returned to avoid confilcts
                            with the notation that [] means all allowed flag values
        """
        if not type(flag)==type(""):
            if flag==[]:
                flag = self.FlagHandler.getValidFlagValues()
            if type(flag)==type(1):
                flag = [flag]

            if removeFlags==[]:
                removeFlags = self.rflags.values()
            if type(removeFlags)==type(1):
                removeFlags = [removeFlags]

            for rflag in removeFlags:
                while rflag in flag:
                    flag.remove(rflag)

            if flag==[]:
                flag = None
        return flag



    #----------------------------------------------------------------------------
    #----- Methods to plot various kinds of data  -------------------------------
    #----------------------------------------------------------------------------
    def signal(self,chanList=[],
               channelFlag=[], plotFlaggedChannels=0, \
               dataFlag=[], plotFlaggedData=0, \
               limitsX=[],limitsY=[], \
               style='l', ci=1, overplot=0, plotMap=0, skynoise=0,
               caption='', subscan=0, noerase=0):
        """
        DES: plot time series of flux density
        INP: (int list) chanList = list of channels
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
             (logical)  skynoise = plot correlated noise (default 0)
             (str)       caption = plot title, default = scan info
             (logical)   subscan = plot vertical lines between subscans 
             (logical)   noerase = do not clear the window
        """

        if plotMap == 1:
            flag='Blank'

        if plotFlaggedChannels:
            dataFlag = self._removeReservedFlagValues(dataFlag, self.rflags['CHANNEL FLAGGED'])
            if dataFlag==None:
                self.MessHand.error("no valid flags")
                return
        
        chanList = self.BolometerArray.checkChanList(chanList, \
                                                     flag=channelFlag,getFlagged=plotFlaggedChannels)
        if len(chanList)<1: 
            self.MessHand.error("no valid channel")
            return
        
        dataX = self.getChanListData('MJD',chanList, \
                                     channelFlag=channelFlag, getFlaggedChannels=plotFlaggedChannels, \
                                     dataFlag=dataFlag, getFlaggedData=plotFlaggedData)
        xLabel = "MJD - MJD(0) [sec]"

        if skynoise:
            dataY = self.getChanListData('skynoise',chanList, \
                                         channelFlag=channelFlag, getFlaggedChannels=plotFlaggedChannels, \
                                         dataFlag=dataFlag, getFlaggedData=plotFlaggedData)
            yLabel = "flux density [arb.u.]"
        else:
            dataY = self.getChanListData('flux',chanList, \
                                         channelFlag=channelFlag, getFlaggedChannels=plotFlaggedChannels, \
                                         dataFlag=dataFlag, getFlaggedData=plotFlaggedData)
            yLabel = "Flux density [Jy]"

        self.MessHand.info("plotting signal")

        if not caption:
            caption = self.ScanParam.caption()
        if not plotMap:
            MultiPlot.plot(chanList,dataX,dataY,\
                           limitsX = limitsX, limitsY = limitsY,\
                           labelX = xLabel, labelY = yLabel, caption=caption,\
                           style=style,ci=ci,overplot=overplot,
                           noerase=noerase)
        else:
            dataX, nNan = compressNan(dataX)
            if nNan != 0:
                nan=1
            else:
                nan=0
            if not caption:
                cap = self.ScanParam.caption()
            Plot.draw(transpose(dataY), \
                      sizeX=[min(dataX),max(dataX)],sizeY=[1,len(chanList)],\
                      limitsX=limitsX,\
                      labelX=xLabel,labelY='Channel #', caption=cap, \
                      wedge=1,nan=nan,noerase=noerase)
        
        if subscan==1:# and not plotMap:
            oldStyle=style
            if limitsY==[]: 
                lower=min(ravel(dataY))
                upper=max(ravel(dataY))
                lower=lower-(upper-lower)*0.2
                upper=upper+(upper-lower)*0.2
            else: 
                upper=limitsY[1]
                lower=limitsY[0]
            mjd=self.getChanData(dataType='MJD',chan=1, flag='None', getFlagged=0,\
                                 flag2=None, subscans=[])
            #mjd= (self.MJD - self.MJD[0])*86400.       
            for i in range(self.ScanParam.SubscanNum[-2]):
                dataY=[]
                dataX=[]
                t=mjd[self.ScanParam.SubscanIndex[1,i]]
                for c in chanList:
                    dataY.append([lower,upper])
                    dataX.append([t,t])   
                MultiPlot.plot(chanList,dataX,dataY,\
                           limitsX = limitsX, limitsY = limitsY,\
                           labelX = xLabel, labelY = yLabel, caption=caption,\
                           style='l',ci=ci+1,overplot=1,\
                           noerase=noerase)
            style=oldStyle         
    #----------------------------------------------------------------------------
    def plotCorrel(self, chanRef=-1, chanList=[],
                   channelFlag=[], plotFlaggedChannels=0, \
                   dataFlag=[], plotFlaggedData=0, \
                   skynoise=0,\
                   limitsX=[], limitsY=[], \
                   style='p', ci=1, overplot=0):
        """
        DES: plot flux density of a list of channels vs. flux density of a
             reference channel
        INP: (int)       chanRef = reference channel number (default : is the first in chanList)
             (int list) chanList = list of channels
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
             (l)        skynoise = plot against the skynoise of chanRef (default :  no)
        """

        chanList = self.BolometerArray.checkChanList(chanList, \
                                                     flag=channelFlag,getFlagged=plotFlaggedChannels)

        if plotFlaggedChannels:
            dataFlag = self._removeReservedFlagValues(dataFlag, self.rflags['CHANNEL FLAGGED'])
            if dataFlag==None:
                self.MessHand.error("no valid flags")
                return
        
        if chanRef == -1:
            chanRef = chanList[0]

	# dataY: remove datapoints where chanRef is flagged
        # retrieve the corresponding channel index
        chanIndex = self.BolometerArray.getChanIndex(chanRef)[0]
	if chanIndex ==  -1:
            self.MessHand.error("chanRef: channel not used") 
            return

        # retrieve flags
        refFlags = self.FlagHandler.getFlags()[:,chanIndex].astype(Int32)
        dataY = self.getChanListData('flux',chanList, \
                                     channelFlag=channelFlag, getFlaggedChannels=plotFlaggedChannels, \
                                     dataFlag=dataFlag, getFlaggedData=plotFlaggedData, \
                                     dataFlag2 = refFlags.astype('i'))

	# dataX: for each chan, remove flagged datapoints
	dataX = []
	for c in chanList:
            chanIndex = self.BolometerArray.getChanIndex(c)[0]
            chanFlags = self.FlagHandler.getFlags()[:,chanIndex].astype(Int32)
            if skynoise:
                dataX.append(self.getChanData('skynoise',chanRef, \
                                              flag=dataFlag, \
                                              getFlagged=plotFlaggedData, \
                                              flag2 = chanFlags.astype('i')))
            else:
                dataX.append(self.getChanData('flux',chanRef, \
                                              flag=dataFlag, \
                                              getFlagged=plotFlaggedData, \
                                              flag2 = chanFlags.astype('i')))


        if skynoise:
            xLabel = "CorrelatedNoise -- channel " + str(chanRef) + " [arb.u.]"
        else:
            xLabel = "flux density -- channel " + str(chanRef) + " [arb.u.]"
            
        yLabel = "flux density [arb.u.]"

        self.MessHand.info("plotting correlation")
        MultiPlot.plot(chanList,dataX,dataY,\
                       limitsX = limitsX, limitsY = limitsY, \
                       labelX = xLabel, labelY = yLabel, caption=self.ScanParam.caption(), \
                       style=style,ci=ci,overplot=overplot)

    #----------------------------------------------------------------------------


    def signalHist(self,chanList=[], \
                   channelFlag=[], plotFlaggedChannels=0, \
                   dataFlag=[], plotFlaggedData=0, \
                   limitsX=[],limitsY=[], \
                   ci=1, overplot=0, caption='', \
                   nbin=60, fitGauss=0, subtractGauss=0, logY=0):
        """
        DES: plot histogram of flux density time series 
        INP: (int list) chanList = list of channels
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
        (int)          nbin = number of bins in histogram
        (l)        fitGauss = fit a gaussian to the data?
        (l)   subtractGauss = subtract gaussian from the data?
        (str)       caption = plot title, default = scan info
        """

        epsilon=0.0001
        ret=[]

        if subtractGauss: fitGauss=1
        
        if plotFlaggedChannels:
            dataFlag = self._removeReservedFlagValues(dataFlag, self.rflags['CHANNEL FLAGGED'])
            if dataFlag==None:
                self.MessHand.error("no valid flags")
                return
        
        chanList = self.BolometerArray.checkChanList(chanList, \
                                                     flag=channelFlag,getFlagged=plotFlaggedChannels)
        chanListIndexes = self.BolometerArray.getChanIndex(chanList)
        
        if len(chanList)<1: 
            self.MessHand.error("no valid channel")
            return
           
        flux = self.getChanListData('flux',chanList, \
                                     channelFlag=channelFlag, getFlaggedChannels=plotFlaggedChannels, \
                                     dataFlag=dataFlag, getFlaggedData=plotFlaggedData)

        nChan=len(chanList)

        extrflux=[]
        for iChan in range(nChan):
            extrflux.append(min(flux[iChan]))
            extrflux.append(max(flux[iChan]))

        minflux=min(extrflux)*0.8
        maxflux=max(extrflux)*0.8
        step=(maxflux-minflux)/nbin

        x=zeros(nbin,'f')
        f=minflux
        for ix in range(len(x)):
            x[ix]=f+step/2.
            f+=step

        dataX=[]
        dataY=[]
        maxval=[]
        for iChan in range(nChan):
            hist=array(fStat.histogram(flux[iChan],minflux,step,nbin))
            dataY.append(hist+epsilon)
            maxval.append(max(hist))
            dataX.append(x)
                
        xLabel = "flux density [arb.u.]"
        yLabel = "distribution of flux"

        if (not subtractGauss):
            
            self.MessHand.info("plotting histogram of signal")

            if logY: limitsY=[0.2,max(maxval)]

            MultiPlot.plot(chanList,dataX,dataY,\
                           limitsX = limitsX, limitsY = limitsY,\
                           labelX = xLabel, labelY = yLabel, \
                           caption=self.ScanParam.caption(),\
                           style='b',ci=ci,overplot=overplot, logY=logY)

            ret=dataY

        if (fitGauss):


            # define a range of x values for a smooth plot of the fit
            npts=200
            x=zeros(npts,'f')
            f=minflux
            step=(maxflux-minflux)/(npts-1)
            for ix in range(len(x)):
                x[ix]=f
                f+=step
                
            dataFit=[]
            plotX=[]
            plotY=[]
            for iChan in range(nChan):
                valid_data_mask = where(dataY[iChan] > 2.*epsilon,1,0)
                xfit=compress(valid_data_mask,dataX[iChan])
                yfit=compress(valid_data_mask,dataY[iChan])
                err=sqrt(yfit)
                result=fitGaussian(xfit,yfit,err)

                dataFit.append(modelgauss(result.params,dataX[iChan]))
                plotX.append(x)
                plotY.append(modelgauss(result.params,x)+epsilon)
                    
            if subtractGauss:

                dataRes=[]
                extrRes=[]
                for iChan in range(nChan):
                    res=(dataY[iChan]-dataFit[iChan])/max(dataY[iChan])
                    dataRes.append(res)
                    extrRes.append(min(res))
                    extrRes.append(max(res))    

                self.MessHand.info("plotting residuals after subtracting gaussian(s)")

                yLabel='residual in units of peak'

                MultiPlot.plot(chanList,dataX,dataRes,\
                               limitsX = limitsX, \
                               limitsY = [min(extrRes)*1.1,max(extrRes)*1.1],\
                               labelX = xLabel, labelY = yLabel, \
                               caption=self.ScanParam.caption(),\
                               style='b',ci=ci,overplot=0)

                ret=dataRes
                    
            else:

                self.MessHand.info("plotting fitted gaussian(s)")
            
                MultiPlot.plot(chanList,plotX,plotY,\
                               limitsX = limitsX, limitsY = limitsY,\
                               labelX = xLabel, labelY = yLabel,\
                               caption=self.ScanParam.caption(),\
                               style='l',ci=2,overplot=1)

                ret=plotY

        return ret
