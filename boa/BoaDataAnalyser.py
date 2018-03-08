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
NAM: BoaDataAnalyser.py (file)
DES: contains BoA channel analyser class
"""
__version__=  '$Revision: 2811 $'
__date__=     '$Date: 2015-03-04 15:32:01 +0100 (Wed, 04 Mar 2015) $'


# TODO Gain should be removed here

#----------------------------------------------------------------------------------
#----- Import ---------------------------------------------------------------------
#----------------------------------------------------------------------------------
import time, os
from Numeric       import *
from boa.Bogli     import Plot, MultiPlot, Forms, DeviceHandler,BogliConfig
from boa.fortran   import fStat, fSNF, fUtilities, fFlag, fBaseline, fMap, fWavelets
from boa.Utilities import Timing, compress2d, Cr2p, Cp2r, principalComponentAnalysis
from boa.Utilities import adjustDataPCA, compressNan, tolist_boa, inPolygon, outPolygon
import arrayfns, FFT

from boa import BoaDataEntity, BoaFlagHandler, BoaConfig


#----------------------------------------------------------------------------------
#----- FFT related function  ------------------------------------------------------
#----------------------------------------------------------------------------------
def applyWindow(inData,window=4,undo=0):
    """
    DES: Utility to apply a window function to a chunk of data
    INP: (f array) inData: the data to apply the function to
              (i) window : the type window function to apply
              1 - Barlett; 2 - Welch ; 3 - Hanning ; 4 - Hamming (default);
              5 - Blackman; 6 - flat-top
              (log) undo : if set, un-apply the window function
    OUT: (f array) the resulting array (same size as input)
    """
    meanY = fStat.f_mean(inData)
    N = len(inData)
    
    #w = arange(N, typecode=Float32)/(N-1)*N
    w = arange(N, typecode=Float32)

    if window == 1:    # Barlett (0 ended)
        w = 1.0 - abs((w-N*1/2)/(N*1/2))
    elif window == 2:  # Welch (0 ended)
        w = 1.0 - ((w-N*1/2)/(N*1/2))**2
    elif window == 3:  # Hanning (0 ended)
        alpha = 0.5; w = alpha - (1-alpha)*cos(2*pi*w/N)
    elif window == 4:  # Hamming
        alpha = 0.53836 ; w = alpha - (1-alpha)*cos(2*pi*w/N)
    elif window == 5:  # Blackman
        w = 0.42-0.5*cos(2*pi*w/N)-0.008*cos(4*pi*w/N)
    elif window == 6:  # flat-top
        w = 1.0 - 1.933*cos(2*pi*w/N) + 1.286*cos(4*pi*w/N) - \
            0.388*cos(6*pi*w/N) + 0.0322*cos(8*pi*w/N)
    else:
        return

    # Compute the amplitude correction factor
    #corr = fStat.f_mean(w)
    corr = sqrt(sum(w**2)/N)
    
    if undo:
        if w[0] == 0 or w[-1] == 0 :
            result = 0.*inData
            # exclude 1st and last elements to avoid division by zero
            #result[1:-1] = (inData[1:-1]-meanY)/w[1:-1]*corr + meanY
            result[1:-1] = inData[1:-1] /w[1:-1] * corr
        else:
            #result = (inData-meanY)/w*corr+meanY
            result = inData / w * corr
    else:
        #result = (inData-meanY)*w/corr+meanY
        result = inData * w / corr

    return result

#----------------------------------------------------------------------------------
#----- FilterFFT class ------------------------------------------------------------
#----------------------------------------------------------------------------------
class FilterFFT:
    """
    DES: To easily do FFT filtering
    INF: make the assumption that the input signal is real, so do not
         care about negative frequencies...
    """

    # from boa.Bogli import Plot

    def __init__(self,X,Y):
        
        # Starting values
        self.X         = array(X)
        self.Y         = array(Y)

        # Flag on the timing precision
        self.Timing    = 0
        self.SamplFreq = 0
        self.__checkTiming()
        self.Timing    = 1

        # ALl the computation will be done on the following data
        self.OutX = array(X)
        self.OutY = array(Y)

        # the number of point used for the FFT
        self.N         = 0
        
        # Result of the FFT
        self.Freq      = 0.0
        self.Amplitude = 0.0
        self.Phase     = 0.0
        self.Power     = 0.0

        # This will contain the datagram when computed
        self.DataGram  = None

    # -----------------------------------------------------------------------------
    def __checkTiming(self):
        """
        DES: check if time sampling is precise enough to allow FFTs
        """
        
        X = self.X
        Y = self.Y

        dX = X[1::]-X[:-1]
        nu = 1./dX
        # use median value (should be most frequent one)
        med = fStat.f_median(nu)
        # compute rms relative to this median
        diff = nu - med
        mean = fStat.f_mean(diff)
        rms  = fStat.f_rms(nu,mean)
        
        # Sampling frequency and error
        f  = med
        df = rms

        # store the sampling Frequency
        self.SamplFreq = f

        # relative error "should" be below 0.1 per cent
        self.Timing = df/f*100 < 0.1
            
        
    # -----------------------------------------------------------------------------
    def __interpolate(self, interval=0):
        """
        DES: linear interpolation of the data to get a regulary gridded dataset
        INP: (f) interval : force the interval (default : median)
        """
        
        X = self.X
        Y = self.Y
        nData = len(X)
        
        if not interval:
            interval = 1./self.SamplFreq

        # let OutX be a little longer than X so that there is
        # no more a 'last point problem' when doing an invFFT
        
        outX = arange(X[0],X[-1]+interval,interval)
        outY = arrayfns.interp(Y,X,outX)

        self.OutX  = outX
        self.OutY  = outY

        self.Timing = 1

    # -----------------------------------------------------------------------------
    def __uninterpolate(self):
        """
        DES: Rebin the data to the initial grid
        """
        
        X         = self.X
        outX      = self.OutX
        outY      = self.OutY
        #self.OutY = (arrayfns.interp(outY,outX,X)).astype(Float32)
        self.Y    = (arrayfns.interp(outY,outX,X)).astype(Float32)

    # -----------------------------------------------------------------------------
    def doFFT(self, interpolate=0, windowing=4, windowSize=0, Xstart=0, Xend=0):
        """
        DES : perform all the necessary steps to do a forward FFT
        INP : (l) interpolate : force an interpolation to be done
                               (default no : check of timing quality better than 0.1 %)
              (i) windowing   : windowing function used to compute FFTs (default: Hamming)
              (i) windowSize  : length of chunks to compute FFT on and to average
                                (default: 0 = compute on the entire data serie)
              (i) Xstart, Xend: optional indices for using only part of the data
        """

        if not self.Timing or interpolate:
            self.__interpolate()
        
        if windowSize:
            fullSize = len(self.OutY)
            # Check that we have enough data points
            if fullSize >= windowSize:
                startPoint = 0
                endPoint   = windowSize
                nbChunk = 0
                while(endPoint <= fullSize):
                    self.doFFT(interpolate=interpolate, windowing=windowing, windowSize=0,\
                               Xstart = startPoint, Xend = endPoint)
                    if startPoint == 0:
                        chunkAmp = self.Amplitude
                        chunkPha = self.Phase
                    else:
                        chunkAmp = chunkAmp + self.Amplitude
                        chunkPha = chunkPha + self.Phase
                    startPoint += windowSize/2
                    endPoint    = startPoint + windowSize
                    nbChunk    += 1
                
                # Average the results by chunk
                self.Amplitude = chunkAmp / float(nbChunk)
                self.Phase     = chunkPha / float(nbChunk)
                # self.Freq already contains the correct Freq array
            else:
                self.doFFT(interpolate=interpolate, windowing=windowing, windowSize=0,\
                           Xstart = Xstart, Xend = Xend)
                
        else:
            self.__forwardFFT(windowing=windowing, Xstart = Xstart, Xend = Xend)

    # -----------------------------------------------------------------------------
    def invFFT(self,windowing=4):
        """
        DES : perform all the necessary steps to do a backward FFT
        """
        self.__backwardFFT()
        if windowing:
            self.OutY = applyWindow(self.OutY,window=windowing,undo=1)
        #self.__uninterpolate()
        self.Y = self.OutY[:len(self.X)]

    # -----------------------------------------------------------------------------
    def doDataGram(self, interpolate=0, n=1024, window=4):
        """
        DES : Compute the Datagram of the data
        INP : (l) interpolate : force an interpolation to be done
                               (default no : check of timing quality better then  0.1 %)
              (i) n           : Number of points for the ffts
        """
        
        if not self.Timing or interpolate:
            self.__interpolate()
        self.__computeDataGram(n=n,window=window)

    # -----------------------------------------------------------------------------
    def __computeDataGram(self, n=1024, window=4):
        """
        DES : Compute the Datagram of the data
        INP  (i) n Number of points for the ffts
        """
        
        # Define the highest power of 2 containing the data
        nData = len(self.OutX)
        # higher power of 2 closer to n to perform the ffts
        N = 2**(int(ceil(log(n)/log(2))))

        # Number of frequency values
        midPoint = int(N/2+1)
        DataGram = zeros((nData-n,midPoint),Float32)

        for index in range(0,nData-n):
            self.__forwardFFT(optimize=1,windowing=window,Xstart=index,Xend=index+n)
            DataGram[index,::] = (sqrt(self.Power)).astype(Float32)

        self.DataGram = DataGram
        
    # -----------------------------------------------------------------------------
    def __forwardFFT(self,optimize=1,windowing=4,Xstart=0,Xend=0):
        """
        DES: perform the FFT
        INP: (i) optimize :  0, will use the full data set
                             1, will zero-pad the data till next power of 2 (default)
             (i) Xstart,Xend: optional indices for using only part of the data
        """

        if Xend:
            X = self.OutX[Xstart:Xend]
            Y = self.OutY[Xstart:Xend]
        else:
            X = self.OutX[Xstart::]
            Y = self.OutY[Xstart::]
            
        if windowing:
            Y = applyWindow(Y,window=windowing)

        # Define the highest power of 2 containing the data
        nData = len(X)
        
        if optimize:
            N = 2**(int(ceil(log(nData)/log(2))))
        else:
            N=nData

        # Save the number of point used for the FFT
        self.N = N

        # Sampling frequency
        midPoint = int(N/2+1)
        samplFreq = self.SamplFreq
        freq = arange(midPoint,typecode='float')*samplFreq/N

        # Do the FFT
        ff   = FFT.real_fft(Y,n=N)

        # # In case we need to switch back to complex ffts...
        # # Shift to have negative frequencies first
        # ff = concatenate([ff[midPoint::],ff[:midPoint]])
        # u = (arange(N,typecode='float')-(midPoint-2))*f/N

        amp, phase = Cr2p(ff)

        # Normalization of the amplitude
        self.Amplitude = amp * sqrt(2.) / float(nData)
        self.Phase     = phase
        self.Freq      = freq
        # PSD: sqrt(self.Power) will be rms/SQRT(Hz)
        self.Power     = self.Amplitude**2 / (samplFreq/N)  # Power density

    # -----------------------------------------------------------------------------
    # compute the inverse FFT 
    def __backwardFFT(self): 
        """
        DES: perform the inverse FFT
        """

        amp       = self.Amplitude
        phase     = self.Phase
        samplFreq = self.SamplFreq
        outX      = self.OutX
        N         = self.N
        
        # Remove the normalization of amplitude
        #amp     = amp / sqrt(2/samplFreq/N)
        amp     = amp / sqrt(2.) * float(len(self.X))
        
        ff = Cp2r(amp,phase)
        outY = FFT.inverse_real_fft(ff,n=N)
        self.OutY = outY[:len(outX)]

    # -----------------------------------------------------------------------------
    def plotFFT(self,plotPhase=0, \
               labelX='Frequency [Hz]', labelY='Amplitude (a.b.u/sqrt(Hz)', \
               limitsX=[],limitsY=[],logX=1, logY=1, overplot=0, ci=1, \
               returnSpectrum=0):
       """
       DES: Plot the fft
       INP: (str)  labelX/Y  : the X/Y label
            (2d f) limitsX/Y : the plot limits for X/Y 
            (bol)  plotPhase : plot phase instead of amplitude (default no)
        (logical)  returnSpectrum : return the values of freq. and amplitude?
                                    (default no)
       """

       X = self.Freq
       Y = self.Amplitude
       if plotPhase:
           Y = self.Phase
           labelY='Phase'
           logY=0

       if logX:
           X = X[1::]
           Y = Y[1::]

       Plot.plot(X,Y, \
                 labelX=labelX,labelY=labelY, \
                 limitsX=limitsX,limitsY=limitsY, \
                 logX = logX, logY = logY, \
                 overplot=overplot, style='l',ci=ci)

       if returnSpectrum:
           return(X,Y)

    # -----------------------------------------------------------------------------
    def plotDataGram(self, interpolate=0, n=1024, window=4,limitsZ=[]):
        """
        DES: Plot the Datagram of the Data
        INP  (i) n Number of points for the ffts
        """

        if not self.DataGram:
            self.doDataGram(interpolate=interpolate,n=n,window=window)
        DataGram  = self.DataGram
        freq      = self.Freq
        X         = self.X
            
        Plot.draw(DataGram,style='idl4',sizeX=[X[0],X[-n]],sizeY=[min(freq),max(freq)], \
                  labelX="Time [s]", labelY="Frequency [Hz]", \
                  limitsZ=limitsZ, wedge=1, caption='sqrt(PSD) [rms/sqrt(Hz)]')
        
        
    # -----------------------------------------------------------------------------
    def blankAmplitude(self, below='?', above='?'):
        """
        DES: blank the amplitude below and/or after a certain frequency
        """

        frequency = self.Freq
        amplitude = self.Amplitude
        
        if above == '?':
            above = min(frequency)
        if below == '?':
            below = max(frequency)
            
        mask = nonzero(where(bitwise_and(frequency >= above,frequency <= below),1,0))

        if mask[0] > 0 and mask[-1] < len(amplitude)-1:
            value = (amplitude[mask[0]-1] + amplitude[mask[-1]+1]) / 2.
        else:
            value = min(amplitude)
            
        if len(mask) > 0:
            put(amplitude,mask,(value))
            
        self.Amplitude = amplitude

    # -----------------------------------------------------------------------------
    def taperAmplitude(self, above='?', N=2):
        """
        DES: Butterworth taper the amplitude above a certain frequency
        INP: (f)  above:  frequency above which to taper
             (f)  N:      steepness parameter
        """

        frequency = self.Freq
        amplitude = self.Amplitude
        
        if above == '?':
            above = min(frequency)

        amplitude = amplitude / sqrt((1.+(frequency/above)**(2*N)))

        self.Amplitude = amplitude

    # -----------------------------------------------------------------------------
    def reduceAmplitude(self, center=50., width=1., factor=10., dB=0):
        """
        DES: multiply the Fourrier spectrum with a filter function
        INP: (f) center: central frequency, in Hz
             (f)  width: window FWHM
             (f) factor: attenuation factor
             (l)    dB : is factor expressed in dB? (default: no)
        """
        frequency = self.Freq
        amplitude = self.Amplitude

        # generate the filter function, sampled at 0.1 Hz
        c    = width / (2.*sqrt(log(2.)))  # convert FWHM
        gaus = exp(-1.*(frequency - center)**2 / c**2)
        filt = 1./(1.+(factor-1.)*gaus)
        #Plot.plot(frequency,self.Amplitude,logY=1)
        self.Amplitude = self.Amplitude * filt
        #Plot.plot(frequency,self.Amplitude,overplot=1,ci=2,style='l')
        #raw_input()
        
#----------------------------------------------------------------------------------
#----- BoA Channel Analyser -------------------------------------------------------
#----------------------------------------------------------------------------------

class DataAna(BoaDataEntity.DataEntity):
  """
  DES: An object of this class is responsible for the flagging of
       individual channels, i.e. it sets the values in the
       Channel_Flag array of the corresponding DataEntity object. It
       provides methods to derive the rms of each channel and to
       automatically search for bad or noisy channels.  Channels might
       be flagged according to a given input file. This object
       provides methods to derive the correlation matrix.

  """

  def __init__(self):
    """
    DES: initialise an instance
    """

    BoaDataEntity.DataEntity.__init__(self)

    # float and integer parameter attributes: 
    self.__statisticsDone = 0     # (logical) are statistics up-to-date?
    self.__corMatrixDone  = 0     # (logical) is correlation Matrix up-to-date ?
    self.__pcaDone        = 0     # (logical) is PCA up-to-date ?

      
  def _coadd(self,other):
      BoaDataEntity.DataEntity._coadd(self,other)
      self.__resetStatistics()
    
  #--------------------------------------------------------------------------------
  #----- Reading: DataEntity.read + auto-flagging  --------------------------------
  #--------------------------------------------------------------------------------
  def read(self,inFile='',febe='',baseband=0,subscans=[],update=0,phase=0, \
           channelFlag=1, integrationFlag=9, \
           blanking=1,readHe=0,readAzEl0=0,readT=0,readWind=0,readBias=0,readPWV=0):
      """
      DES: fill a BoA data object from an MB-FITS file
      INP:     (string) inFile : path to the dataset to be read
                 (string) febe : FE-BE name to select
                (int) baseband : baseband to select
           (int list) subscans : subscan numbers to read (default: all)
              (logical) update : if true, do not reset previous entity object
                   (int) phase : phase to be stored (default: phase diff)
                (log) blanking : automatic flagging of blanked data (default: yes)
          channelFlag (i list) : flag for not connected feeds (default: 1 'NOT CONNECTED')
      integrationFlag (i list) : flag for blanked integrations (default: 9 'BLANK DATA')
                (log) readHe   : do we need the He3 temparatures? (default: no)
               (log) readAzEl0 : do we read monitor Az,El at start? (default: no)
             (logical)   readT : do we read T_amb from monitor? (def: no)
            (logical) readWind : do we read wind speed, dir... (def: no)
            (logical)  readPWV : do we read pwv? (def: no)
            (logical) readBias : do we need ASZCa bias settings? (def: no)
      OUT:        (int) status : 0 if reading ok, <> 0 if an error occured
           Possible error codes are:
                -1 = file could not be openned
                -2 = something wrong with FEBE
                -3 = something wrong with basebands
                -4 = something wrong with subscans
       """
      status = BoaDataEntity.DataEntity.read(self,inFile,febe,baseband,subscans, \
                                             update,phase,channelFlag,integrationFlag, \
                                             readHe,readAzEl0,readT,readWind,readBias,readPWV)

      if status:  # if reading failed
        return status
      
      nUsedChannels = self.BolometerArray.NUsedChannels
      nSub          = self.ScanParam.NObs
      nInt          = self.ScanParam.NInt

      dataShape = shape(self.Data)
      nUsedChannels = dataShape[1]

      self.CorrelatedNoise = fUtilities.as_column_major_storage(zeros((nInt,nUsedChannels), Float32))

      self.FFCF_CN         = ones((nUsedChannels,nUsedChannels), Float32)
      self.CorMatrix       = ones((nUsedChannels,nUsedChannels), Float32)

      self.FFCF_Gain       = ones(nUsedChannels, Float32)
      self.FF_Median       = ones(nUsedChannels, Float32)

      self.ChanRms         = ones(nUsedChannels,Float32)
      self.ChanMean        = ones(nUsedChannels,Float32)
      self.ChanMed         = ones(nUsedChannels,Float32)

      self.ChanRms_s       = ones((nUsedChannels,nSub),Float32)
      self.ChanMean_s      = ones((nUsedChannels,nSub),Float32)
      self.ChanMed_s       = ones((nUsedChannels,nSub),Float32)

      return status

  #--------------------------------------------------------------------------------
  #----- flagging methods ---------------------------------------------------------
  #--------------------------------------------------------------------------------
  def flagRms(self,chanList=[],below=0,above=1e10,flag=2): 
    """
    DES: Flag channels with rms below 'below' or above 'above'
    INP: (i list) chanList : list of channel to flag (default: current list)
         (f)      below    : flag channels with rms < 'below'
         (f)      above    : flag channels with rms > 'above'
         (i)       flag    : flag value to set (default: 2 'BAD SENSITIVITY')
    """

    self.MessHand.debug("flagRms start...")

    # recompute stat if needed
    if not self.__statisticsDone:
      self.__statistics()

    # check channel list
    chanList = self.BolometerArray.checkChanList(chanList)
    # this should have returned a sublist of valid channels or all valid channels
    # if input chanList was empty
    if len(chanList)<1: 
      self.MessHand.error("no valid channel")
      return

    chanRms = array(self.getChanListData('rms',chanList))

    badChan = []
    badChan.extend(compress(where(chanRms < below,1,0),chanList))
    badChan.extend(compress(where(chanRms > above,1,0),chanList))

    if badChan:
      self.flagChannels(chanList=badChan,flag=flag)
      
    self.MessHand.debug("...flagRMS ends")

  #--------------------------------------------------------------------------------
  def flagFractionRms(self,chanList=[],ratio=10.,flag=2,plot=0,above=1,below=1):
      """
      DES: flag according to rms, with limits depending on median rms
      INP: (i list) chanList : list of channel to flag (default: current list)
           (f)         ratio : channels with rms below median/ratio and above
                               median*ratio will be flagged
           (i)          flag : value of flag to set (default: 2 'BAD SENSITIVITY')
           (b)          plot : plot the results
           (b)         above : should we flag above median * ratio? (default yes)
           (b)         below : should we flag below median / ratio? (default yes)
      """
      self.MessHand.debug("flagFractionRms start...")

      # recompute stat if needed
      if not self.__statisticsDone:
          self.__statistics()

      # check channel list
      chanList = self.BolometerArray.checkChanList(chanList)
      if len(chanList)<1: 
          self.MessHand.error("no valid channel")
          return

      chanRms = array(self.getChanListData('rms',chanList))
      med = fStat.f_median(chanRms)

      if plot:
          self.plotRmsChan()
          Forms.shadeY(med*ratio,max(chanRms),ci=2)
          Forms.shadeY(min(chanRms),med/ratio,ci=2)

      if below:
          self.flagRms(chanList,below=med/ratio,flag=flag)
      if above:
          self.flagRms(chanList,above=med*ratio,flag=flag)

      self.MessHand.debug("...flagFractionRMS ends")
      
  #--------------------------------------------------------------------------------

  
  def flagAutoRms(self,chanList=[],threshold=3.,flag=2): 
    """
    DES: Automatic flagging of channels, based on their rms
    INP: (i list) chanList : list of channel to flag (default: current list)
         (f)     threshold : flag outliers channels w.r.t. threshold
         (i)          flag : flag value to set  (default: 2 'BAD SENSITIVITY')
    
    """

    self.MessHand.debug("flagRms start...")

    # recompute stat if needed
    if not self.__statisticsDone:
      self.__statistics()

    # check channel list
    chanList = self.BolometerArray.checkChanList(chanList)
    # this should have returned a sublist of valid channels or all valid channels
    # if input chanList was empty
    if len(chanList)<1: 
      self.MessHand.error("no valid channel")
      return

    chanRms = array(self.getChanListData('rms',chanList))
    meanRms, sigmaRms, medianSDev, medianMDev = fStat.f_stat(chanRms)
    
    UsedChannels = self.BolometerArray.UsedChannels

    badChan = []
    badChan.extend(compress(where(chanRms < medianSDev - threshold * sigmaRms,1,0),chanList))
    badChan.extend(compress(where(chanRms > medianSDev + threshold * sigmaRms,1,0),chanList))

    if badChan:
      self.flagChannels(chanList=badChan,flag=flag)
      
    self.MessHand.debug("...flagRMS ends")

#--------------------------------------------------------------------------------
  
  def flagChannels(self,chanList=[],flag=8): 
    """
    DES: assign flags to a list of channels
    INP: (i list) chanList : list of channels to be flagged (default current list)
         (i list)     flag : flag values (default: 8 'TEMPORARY')
    """

    self.MessHand.debug("flagChannels start...")

    chanList = self.BolometerArray.checkChanList(chanList)

    if len(chanList)<1: 
        self.MessHand.error("no valid channel")
        return

    # Flag in the BolometerArray class and...
    nbFlag = self.BolometerArray.flag(chanList,flag=flag)

    # ... report that on self.DataFlag
    chanIndexes = self.BolometerArray.getChanIndex(chanList)
    #for index in chanIndexes:
    #    if index != -1:
    #        if self.BolometerArray.FlagHandler.isSetOnIndex(index):
    #            self.FlagHandler.setAll(self.rflags['CHANNEL FLAGGED'], dim=1, index=index)
    # Buggy ! corrected 2017-07-04
    for i in range(len(chanIndexes)):
        if chanIndexes[i] != -1:
            if self.BolometerArray.FlagHandler.isSetOnIndex(chanList[i]-1):
                self.FlagHandler.setAll(self.rflags['CHANNEL FLAGGED'], dim=1, index=chanIndexes[i])

    self.MessHand.debug("...flagChannels ends")

  #--------------------------------------------------------------------------------
  def unflagChannels(self,chanList=[],flag=[]): 
    """
    DES: unflags a list of channels
    INP: (i list) chanList : list of channels to be unflagged (default current list)
         (i list)     flag : flag values (default []: unset all flags)
    """

    self.MessHand.debug("unflagChannels start...")
    
    if type(chanList) == int:
        chanList = [chanList]

    if chanList == []:
        chanList = self.BolometerArray.UsedChannels
    if len(chanList)<1: 
        self.MessHand.error("no valid channel")
        return

    # Unflag in the BolometerArray class and...
    nbFlag = self.BolometerArray.unflag(chanList,flag=flag)

    # ... report that on self.DataFlag
    chanIndexes = self.BolometerArray.getChanIndex(chanList)
    #for index in chanIndexes:
    #    if index != -1:
    #        if self.BolometerArray.FlagHandler.isUnsetOnIndex(index):
    #            self.FlagHandler.unsetAll(self.rflags['CHANNEL FLAGGED'],
    #                                      dim=1, index=index)
    # Buggy ! corrected 2017-07-04
    for i in range(len(chanIndexes)):
        if chanIndexes[i] != -1:
            if self.BolometerArray.FlagHandler.isUnsetOnIndex(chanList[i]-1):
                self.FlagHandler.unsetAll(self.rflags['CHANNEL FLAGGED'],
                                          dim=1, index=chanIndexes[i])

    self.MessHand.debug("...unflagChannels ends")

  #--------------------------------------------------------------------------------
  def getFlaggedChannels(self):
      """
      Function which returns the list of channels currently flagged.
      """
      nonFlagged = self.BolometerArray.checkChanList([])  # List of non-flagged
      nbTotal = self.BolometerArray.NChannels
      result = range(1,nbTotal+1)
      for chan in nonFlagged:
          result.remove(chan)

      return result

  # ---------------------------------------------------------------------
  def flagRCP(self,rcpFile,flag=1):
      """
      NAM: flagRCP
      DES: flag channels not present in the given RCP file
      INP: (string) rcpFile: name of input RCP file
           (int)       flag: value used to flag channels (def.: 1)
      """
      try:
          f = file(os.path.join(BoaConfig.rcpPath,rcpFile))
      except IOError:
          self.MessHand.error("could not open file %s"%(rcpFile))
          return
      
      param = f.readlines()
      f.close()

      channels = []
      for i in range(len(param)):
          if param[i][0] not in ['!','#','n']:    # skip comments
              tmp = string.split(param[i])
              num = string.atoi(tmp[0])
              channels.append(num)

      missing = []
      for i in range(self.BolometerArray.NChannels):
          if i+1 not in channels:
              missing.append(i+1)

      if missing:
          self.flagChannels(chanList=missing,flag=flag)

  #--------------------------------------------------------------------------------

  def deglitch_old(self,chanList=[],window=10,above=5,flag=1,maxIter=10,minTimeSampInSubscan=100):
      """
      DES: Flag yet unflagged data where glitches occur
           IT IS HIGHLY RECOMMENDED TO REMOVE SKYNOISE BEFORE DEGLITCHING.
      INP: (i list) chanList : list of channel to flag (default: current list)
                (int) window : compute sliding rms in this window
                   (f) above : flag data where the sliding rms > 'above'*rms
               (i list) flag : flag values (default: 1 'SPIKE')
      """

  
      chanList = self.BolometerArray.checkChanList(chanList)
      chanListIndices = self.BolometerArray.getChanIndex(chanList)
      setFlags=take(self.FlagHandler.getFlags(),chanListIndices,axis=1)
      setFlags=take(setFlags,range(shape(setFlags)[0]-window),axis=0)
      
      slideRms0=self.slidingRms(nbInteg=window,channel=[],getFlagged=0,flag='None')
      
      #set flags on slideRms0
      mask=where(setFlags > 0,float('NaN'),1)
      slideRms0=slideRms0*mask
      
      meanRms=slideRms0[:,0]*0.0

      for ts in range(len(slideRms0[:,0])):
          timeSlice=slideRms0[ts,::]
          timeSlice,nnan=compressNan([timeSlice])
          meanRms[ts]=fStat.f_mean(timeSlice)

      addBefore=int(window/2.)
      addAfter=window-addBefore
      slideRms=(array(range(addBefore))*0.0+meanRms[0]).tolist()
      slideRms.extend(meanRms[::].tolist())
      slideRms.extend((array(range(addAfter))*0.0+meanRms[shape(meanRms)[0]-1]).tolist())
      
      slideRms=array(slideRms)

      done=0
      counter=0
      totalflags=array(range(shape(slideRms)[0]))*0

      while(done == 0):
          if (counter >= maxIter):
              break
          median = fStat.f_median(slideRms)
          offsets = array(slideRms)-median
          pos_mask = where(offsets > 0,1,0)
          offsets = compress(pos_mask,offsets)
          std_dev = fStat.f_mean(offsets)
          hiVal = median+above*std_dev
          mask = where(array(slideRms) >= hiVal,1,0)
          if (not(1 in mask)):
              done=1
          else:
            counter+=1
          putmask(totalflags,mask,1)
          putmask(slideRms,mask,float('NaN'))

      n_flagged=len(nonzero(totalflags))

      if n_flagged > 0:
        
          totalflags = convolve(totalflags,(array(range(window*2))*0+1).tolist(),mode=1)
          mask = where(totalflags > 0,1,0)

          subscan_start=[0]
          subscan_end=[]
          oldFlagSet=mask[0]
          for i in range(len(mask)):
              if (oldFlagSet != mask[i]):
                  oldFlagSet = mask[i]
                  subscan_end.extend([i-1])
                  subscan_start.extend([i])
                  
          subscan_end.extend([len(mask)-1])                  
          subscan_len=array(subscan_end)-array(subscan_start)
        
          for i in range(len(subscan_len)):
              if ((subscan_len[i] < minTimeSampInSubscan)):
                  mask[subscan_start[i]:subscan_end[i]+1]=1
          
          # put the flags on the ScanParam flag array
          self.ScanParam.FlagHandler.setOnMask(mask,iFlags=flag)
          # put the flags on the main flag array
          for chan in chanListIndices:
              self.FlagHandler.setOnMask(mask, self.rflags['INTEGRATION FLAGGED'], \
                                         dim=1, index=chan)
              
          nflagged=shape(compress(mask,mask))[0]
          self.MessHand.info("%5i timestamps flagged with flag %s" % (nflagged, str(flag)))

          self.__resetStatistics()

      return n_flagged


  #--------------------------------------------------------------------------------

  def deglitch2(self,chanList='all',above=10,flag=1,maxIter=10,minTimeSampInSubscan=200):
      # EXPERIMENTAL BUT FAST

      # check chan list
      chanList = self.BolometerArray.checkChanList([])
      chanListIndices = self.BolometerArray.getChanIndex(chanList)

      # get flux
      flux = array(self.getChanListData('flux',chanList=chanList,\
                                        dataFlag='None'))
      
      # get flags, select those corresponding to chanList
      timeFlags = self.ScanParam.FlagHandler.getFlags()
      dataFlags = self.FlagHandler.getFlags()
      dataFlags = take(dataFlags,chanListIndices,1)
      
      absflux=abs(flux)
      meanabs=sum(absflux,0)/(float(len(chanListIndices)))
      med=fStat.f_median(meanabs)
      meanabs=meanabs-med

      done=0
      counter=0
      totalflags=array(range(shape(flux)[1]))*0
      while(done == 0):
          if (counter >= maxIter):
              break
          m=fStat.f_median(meanabs)
          #compute the MEDIAN deviation!
          s=fStat.f_median(abs(meanabs))
      
          mask=where((meanabs > m+above*s),1,0)

          if (not(1 in mask)):
              done=1
          else:
              counter+=1
              putmask(totalflags,mask,1)
              putmask(meanabs,mask,m)
              
      n_flagged=len(nonzero(totalflags))

      nflagged=0
      if n_flagged > 0:
        
          totalflags = convolve(totalflags,(array(range(20))*0+1).tolist(),mode=1)
          mask = where(totalflags > 0,1,0)

          subscan_start=[0]
          subscan_end=[]
          oldFlagSet=mask[0]
          for i in range(len(mask)):
              if (oldFlagSet != mask[i]):
                  oldFlagSet = mask[i]
                  subscan_end.extend([i-1])
                  subscan_start.extend([i])
                  
          subscan_end.extend([len(mask)-1])                  
          subscan_len=array(subscan_end)-array(subscan_start)
        
          for i in range(len(subscan_len)):
              if ((subscan_len[i] < minTimeSampInSubscan)):
                  mask[subscan_start[i]:subscan_end[i]+1]=1
          
          # put the flags on the ScanParam flag array
          self.ScanParam.FlagHandler.setOnMask(mask,iFlags=flag)
          # put the flags on the main flag array
          for chan in chanListIndices:
              self.FlagHandler.setOnMask(mask, self.rflags['INTEGRATION FLAGGED'], \
                                         dim=1, index=chan)
              
          nflagged=shape(compress(mask,mask))[0]
          self.MessHand.info("%5i timestamps flagged with flag %s" % (nflagged, str(flag)))
          self.MessHand.info("deglitch: %4i iterations" % (counter))

          self.__resetStatistics()

      return nflagged


  #--------------------------------------------------------------------------------
 

  def deglitch(self,chanList='all',above=5,flag=1,maxIter=10,window=20,minTimeSampInSubscan=200,plot=0):
      """
      DES: Flag yet unflagged data where glitches occur. Iterative method.
      INP: (i list) chanList : list of channels to flag (default: all valid channels)
                   (f) above : flag data where the sliding rms > 'above'*rms (defualt 5)
               (i list) flag : flag values (default: 1 'SPIKE')
                 (i) maxIter : maximum number of iterations (default 10)
                  (i) window : flag only in windows of this many time stamps (default 20)
    (i) minTimeSampInSubscan : minimum allowed number of time samples in continuous unflagged regions
      """

      # get ALL data
      chanList = self.BolometerArray.checkChanList(chanList)
      chanListIndices = self.BolometerArray.getChanIndex(chanList)
      
      flux = array(self.getChanListData('flux',chanList=chanList,\
                                        dataFlag='None'))

      # normalise by channel rms
      rms=self.getChanListData('rms',chanList=chanList)

      for i in range(len(chanList)):
          flux[i,::]=(flux[i,::]/rms[i]).astype(flux.typecode())
      
      # get time series flags
      timeFlags = self.ScanParam.FlagHandler.getFlags()
      dataFlags = self.FlagHandler.getFlags()
      dataFlags = take(dataFlags,chanListIndices,1)
      # initialize array of standard deviations
      sdev = array(range(shape(flux)[1])).astype('f')
      sdev[::] = 0.#float('NaN')
      for i in range(len(sdev)):
          if (timeFlags[i] == 0):
              fluxmask=where(dataFlags[i,::] == 0, 1, 0)
              fl=flux[:,i]
              fl=compress(fluxmask,fl)
              if (shape(fl)[0] > 3):
                  mean=fStat.f_mean(fl)
                  sdev[i]=fStat.f_rms(fl,mean)

      done=0
      counter=0
      totalflags=array(range(shape(sdev)[0]))*0

      while(done == 0):
          if (counter >= maxIter):
              break
          mask=where(sdev > 0.,1,0)
          sdev_good=compress(mask,sdev)
          sdev_mean=fStat.f_mean(sdev_good)
          sdev_rms=fStat.f_rms(sdev_good,sdev_mean)
          sdev_median=fStat.f_median(sdev_good)
          hiVal = sdev_mean+above*sdev_rms
          mask = where(sdev >= hiVal,1,0)
          if (not(1 in mask)):
              done=1
          else:
              counter+=1
              putmask(totalflags,mask,1)
              putmask(sdev,mask,0.)

      n_flagged=len(nonzero(totalflags))
      nflagged=0

      if n_flagged > 0:
        
          totalflags = convolve(totalflags,(array(range(window))*0+1).tolist(),mode=1)
          mask = where(totalflags > 0,1,0)

          subscan_start=[0]
          subscan_end=[]
          oldFlagSet=mask[0]
          for i in range(len(mask)):
              if (oldFlagSet != mask[i]):
                  oldFlagSet = mask[i]
                  subscan_end.extend([i-1])
                  subscan_start.extend([i])
                  
          subscan_end.extend([len(mask)-1])                  
          subscan_len=array(subscan_end)-array(subscan_start)
        
          for i in range(len(subscan_len)):
              if ((subscan_len[i] < minTimeSampInSubscan)):
                  mask[subscan_start[i]:subscan_end[i]+1]=1
          
          # put the flags on the ScanParam flag array
          self.ScanParam.FlagHandler.setOnMask(mask,iFlags=flag)
          # put the flags on the main flag array
          for chan in chanListIndices:
              self.FlagHandler.setOnMask(mask, self.rflags['INTEGRATION FLAGGED'], \
                                         dim=1, index=chan)
              
          nflagged=shape(compress(mask,mask))[0]
          self.MessHand.info("%5i timestamps flagged with flag %s" % (nflagged, str(flag)))
          self.MessHand.info("deglitch: %4i iterations" % (counter))

          self.__resetStatistics()

      return nflagged

  #--------------------------------------------------------------------------------
 

  def glwDetect(self,chanList=[],scale=5,nsigma=5,window=25,plotCh='?',collapse=1,updateFlags=0):
      """
      NAM: glwDetect (function)
      DES: detects glitchy time intervals using wavelets
      INP: (f array) data structure
           (i list) chanList = list of channels to consider [def. all]
	   (i)         scale = wavelet scale considered [def. 5]
	   (i)        nsigma = [def. 5]
	   (i)        window = window size for flag smoothing [def. 25]
	   (i)        plotCh = if set plot the result for channel plotCh
	   (b)      collapse = whether to collapse channel flags together [def. 1]
	   (b)   updateFlags = whether to update the data flags accordingly [def. 0]
      OUT: Mask of channels to be flagged
      """    

      self.MessHand.debug('glwDetect start...')

      # check channel list
      chanList = self.BolometerArray.checkChanList(chanList)
      if len(chanList)<1: 
        self.MessHand.error("no valid channel")
      	return 0

      nchan=len(chanList)
   
      # Create time step flags
      chFlags=zeros([self.ScanParam.NInt,self.BolometerArray.NChannels],typecode='int4')    
      for ch in chanList:
	nFl=0
        Flags=zeros(size(self.ScanParam.MJD))
	ts = array(self.Data[:,ch-1])
	wtrans = fWavelets.wavetransform1d(ts,scale+1)
	wmax = wtrans[:,scale-1]
	Mean, Med, SDev, MDev = fStat.arraystat(wmax,Flags)
        hiVal = Med+nsigma*MDev
        loVal = Med-nsigma*MDev
	Flags = where(bitwise_or(wmax >= hiVal,wmax <= loVal),1,0)
	while (nFl < sum(chFlags)):
		nFl=sum(chFlags)
		Mean, Med, SDev, MDev = fStat.arraystat(wmax,Flags)
        	hiVal = Med+nsigma*MDev
        	loVal = Med-nsigma*MDev
		Flags = where(bitwise_or(wmax >= hiVal,wmax <= loVal),1,0)
	chFlags[:,ch-1] = Flags.astype(type(chFlags))
    
      # Smooth (and eventually collapse) flags
      kernel  = zeros(window) + 1
      if collapse:
    	chFlags = where(sum(transpose(chFlags)) < 0.2*float(nchan),0,1)
    	chFlags = convolve(chFlags,kernel,mode=1)
      else:
    	for ch in chanList:
		chFlags[:,ch-1] = convolve(chFlags[:,ch-1],kernel,mode=1)
	
      # Plot the channel chPlot with flagged time intervals
      if (plotCh != '?'):
        MJD = self.ScanParam.get(dataType='mjd',flag='None')
	if collapse:
		Flags = chFlags
	else:
		Flags = chFlags[:,plotCh-1]
	d = array(self.getChanData('flux',plotCh,flag='None'))
	Mean, Med, SDev, MDev = fStat.arraystat(d,Flags)
	pl = where(Flags < 0.5,d,0)
        Plot.plot(MJD,d-Med ,style='l',ls=1,labelX='MJD - MJD(0) [sec]',labelY='flux density [arb.u.]',limitsY=[min(pl-Med)*2.,max(pl-Med)*2.])
        Plot.plot(MJD,pl-Med,style='l',overplot=1,ci=2)
   
	
      # Update flags
      if updateFlags:
    	if collapse:
		chans = zeros(self.BolometerArray.NChannels) + 1
		Flags = multiply.outer(chFlags,chans)
		self.FlagHandler.setOnMask(Flags,1)
	else:
		self.FlagHandler.setOnMask(chFlags,1)
	    
      self.MessHand.debug('... glwDetect end')

      return chFlags

#--------------------------------------------------------------------------------

  def despike(self,chanList=[],below=-5,above=5,flag=1):
      """
      DES: Flag yet unflagged data below 'below'*rms and above 'above'*rms.
      INP: (i list) chanList : list of channel to flag (default: current list)
           (f)      below    : flag data with value < 'below'*rms
           (f)      above    : flag data with value > 'above'*rms
           (i list)     flag : flag values (default: 1 'SPIKE')
      """
    
      self.MessHand.debug('despike start...')

      # check channel list
      chanList = self.BolometerArray.checkChanList(chanList)
      if len(chanList)<1: 
          self.MessHand.error("no valid channel")
          return 0

      # recompute stat if needed
      if not self.__statisticsDone:
          self.__statistics()

      flag = self._removeReservedFlagValues(flag)
      if flag==None:
          self.MessHand.error("no valid flags")
          return

      # process the DataFlags array
      totalFlag = 0
    
      for chan in chanList:
          # remember that physical channel numbers do not correspond to
          # index of Data*, so...
          chanIndex = self.BolometerArray.getChanIndex(chan)[0]

          # dataType (unflagged)
          dataTest = self.getChanData(dataType='flux',chan=chan,flag='Blank')

          # We want to flag all the data which do not have this
          # particular flag, so check what is left
          dataTest_noFlag = ravel(self.getChanData(dataType='flux', chan=chan, flag=flag))

          if len(dataTest_noFlag)>0:
              rms  = self.getChanData('rms',chan=chan)
              mean = self.getChanData('mean',chan=chan)

              # ...per subscan
              # to work on subscan, one need to wait for subscanIndex fixes
              # a loop on subscans will the been needed to flag properly

              # rms = self.getChanData('rms_s',chan=chan,flag=flag)
              # mean = self.getChanData('mean_s',chan=chan,flag=flag)

              hiVal = mean+above*rms
              loVal = mean+below*rms

              mask = where(bitwise_or(dataTest >= hiVal,dataTest <= loVal),1,0)

              if len(nonzero(mask)) > 0:
                  n0 = self.FlagHandler.nSet(flag, dim=1, index=chanIndex)
                  self.FlagHandler.setOnMask(mask, flag, dim=1, index=chanIndex)
                  n1 = self.FlagHandler.nSet(flag, dim=1, index=chanIndex)
                  totalFlag += (n1-n0)
                  if (n1-n0):
                      self.MessHand.debug("Channel %i"%chan +\
                                          " %5i timestamps flagged" % (n1-n0) + \
                                          " with flag %s" % str(flag))
                  else:
                      self.MessHand.debug("Channel %i "%chan+\
                                          " nothing to flag")

              
      if totalFlag > 0:
          self.MessHand.info("%5i samples flagged with flag %s" % (totalFlag, str(flag)))
          self.__resetStatistics()
      else:
          self.MessHand.warning("Nothing flagged")
      
      self.MessHand.debug("... despike end")
    
      return totalFlag

    

  #--------------------------------------------------------------------------------

  def iterativeDespike(self,chanList=[],below=-5,above=5,flag=1,maxIter=100):
    """
    DES: Iteratively flag yet unflagged data below 'below'*rms and above 'above'*rms.
    INP: (i list) chanList : list of channel to flag (default: current list)
         (f)      below    : flag data with value < 'below'*rms
         (f)      above    : flag data with value > 'above'*rms
         (i)      maxIter  : maximum number of iteration (default 100)
         (i list)     flag : flag values (default: 1 'SPIKE')
    """
    
    # Initialize loop variables
    despiked=1
    i=0
    
    # Run loop
    while ((despiked > 0) and (i < maxIter)):
        i += 1
	despiked = self.despike(chanList=chanList,below=below,above=above,flag=flag)
    
  #--------------------------------------------------------------------------------

  def unflag(self,channel=[],flag=[]):
      """
      NAM: unflag (method)  
      DES: Unflag data, i.e. reset to 0.
      INP: (i list) channel : list of channel to flag (default: current list)
           (i)         flag : unflag only this value (default []: all non-reserved flag values)
 
      """

      self.MessHand.debug('unflag start...')

      # check channel list, special one, since even flagged channel has
      # to be taken into account
      
      CurrChanList = self.BolometerArray.CurrChanList
      UsedChannels = self.BolometerArray.UsedChannels

      if channel in ['all','al','a']:
          channel = UsedChannels
      elif channel == []:
          channel = CurrChanList
      elif type(channel) == type(1):
          channel = [channel]
      
      chanList = []
      for num in channel:
          if num in UsedChannels:
              chanList.append(num)

      if len(chanList)<1: 
          self.MessHand.error("no valid channel")
          return

      flag = self._removeReservedFlagValues(flag)
      if flag==None:
          self.MessHand.error("no valid flags")
          return

      totalUnFlag=0
      for chan in chanList:

          chanIndex = self.BolometerArray.getChanIndex(chan)[0]

          n0 = self.FlagHandler.nUnset(flag, dim=1, index=chanIndex)
          self.FlagHandler.unsetAll(flag, dim=1, index=chanIndex)
          n1 = self.FlagHandler.nUnset(flag, dim=1, index=chanIndex)
          totalUnFlag += (n1-n0)
    
      if totalUnFlag > 0:
          self.MessHand.info("%5i samples unflagged with flag %s" % (totalUnFlag,str(flag)))
          self.__resetStatistics()
      else:
          self.MessHand.warning("Nothing unflagged")

      self.MessHand.debug('... unflag end')

  #----------------------------------------------------------------------------
  def flag(self,dataType='', channel='all', below='?', above='?', flag=8):
      """
      DES: flag data based on dataType, general flagging routine, may be slow
      INP: (s)     dataType : flag based on this dataType
           (i list) channel : list of channel to flag (default: all)
           (f)        below : flag dataType < below (default max; or 5*RMS)
           (f)        above : flag dataType > above (default min; or -5*RMS)
           (i)         flag : flag value (default 8 'TEMPORARY')

          below and above should be in unit of the flagged data,
          except for 'Lon' and 'Lat' where they should be in arcsec
      """
    
      self.MessHand.debug('flag start...')
    
      # check channel list
      chanList = self.BolometerArray.checkChanList(channel)
      if len(chanList)<1: 
          self.MessHand.error("no valid channel")
          return

      if channel in ['all', 'al', 'a'] and \
             ( dataType in ['LST','lst'] or \
               dataType in ['MJD','mjd'] or \
               dataType in ['UT'] \
               ):
      
          # this is a channel independant and ScanParam based flagging ...
          self.MessHand.warning("You should use the flagInTime method")
          self.flagInTime(dataType=dataType,below=below,above=above,flag=flag)
          return

      flag = self._removeReservedFlagValues(flag)
      if flag==None:
          self.MessHand.error("no valid flags")
          return

      # process the DataFlags array
      totalFlag = 0
      for chan in chanList:
          self.MessHand.debug("flaging for chan: "+str(chan))

          # remember that physical channel numbers do not correspond to
          # index of Data*, so...
          chanIndex = self.BolometerArray.getChanIndex(chan)[0]

          # dataType (unflagged)
          dataTest = self.getChanData(dataType=dataType,chan=chan,flag='None')

          # We want to flag all the data which do not have this
          # particular flag, so check what is left
          dataTest_noFlag = ravel(self.getChanData(dataType=dataType, chan=chan, flag=flag))

          if len(dataTest_noFlag)>0:
              self.MessHand.debug(" found something to flag")

              hiVal = max(dataTest_noFlag)
              loVal = min(dataTest_noFlag)

              # default inputs
              if above != '?':
                  loVal = above
              if below != '?':
                  hiVal = below

              mask = where(bitwise_and(dataTest >= loVal,dataTest <= hiVal),1,0)

              if len(nonzero(mask)) > 0:
                  n0 = self.FlagHandler.nSet(flag, dim=1, index=chanIndex)
                  self.FlagHandler.setOnMask(mask, flag, dim=1, index=chanIndex)
                  n1 = self.FlagHandler.nSet(flag, dim=1, index=chanIndex)
                  totalFlag += (n1-n0)
                  if (n1-n0):
                      self.MessHand.debug("Channel %i"%chan +\
                                          " %5i timestamps flagged" % (n1-n0) +\
                                          " with flag %s" % str(flag))
                  else:
                      self.MessHand.debug("Channel %i "%chan+\
                                          " nothing to flag")

      if totalFlag > 0:
          self.MessHand.info("%5i samples flagged with flag %s" % (totalFlag, str(flag)))
          self.__resetStatistics()
      else:
          self.MessHand.warning("Nothing flagged")

      self.MessHand.debug("... flag end")

  #--------------------------------------------------------------------------------
  def flagInTime(self,dataType='LST', below='?', above='?', flag=8):
      """
      DES: Flag data in time interval
      INP: (float)    below    = flag data below this value (default end of the scan) 
           (float)    above    = flag data above this value (default start of the scan)
           (int)      flag     = flag to be set (default: 8 'TEMPORARY')
      """
      self.MessHand.debug('flagInTime start...')
    
      # check channel list
      chanList = self.BolometerArray.checkChanList("all", flag='None')
      if len(chanList)<1: 
          self.MessHand.error("no valid channel")
          return

      # return Data index
      chanListIndexes = self.BolometerArray.getChanIndex(chanList) 

      # check dataType
      if ( dataType not in ['LST','lst'] and
           dataType not in ['MJD','mjd'] and
           dataType not in ['UT'] and
           dataType not in ['speed','Speed'] and
           dataType not in ['accel','Accel']):
          self.MessHand.error("use only this method to flag in time")
          return

      # get the default values
      dataTest = self.ScanParam.get(dataType=dataType,flag='None')
      if above == '?':
          above = min(dataTest)
      if below == '?':
          below = max(dataTest)

      self.ScanParam.flag(dataType=dataType,below=below,above=above,flag=flag)
      # Report the flag on the ScanParam attribute
      mask = self.ScanParam.FlagHandler.isSetMask()
      if len(nonzero(mask)) > 0:
          for chan in chanListIndexes:
              self.FlagHandler.setOnMask(mask, \
                                         self.rflags['INTEGRATION FLAGGED'], \
                                         dim=1, index=chan)
          self.__resetStatistics()
      else:
          self.MessHand.warning("Nothing flagged")
      
      self.MessHand.debug("... flagInTime end")
    
  #--------------------------------------------------------------------------------
  def unflagInTime(self,dataType='LST', below='?', above='?', flag=[]):
      """
      DES: Unflag data in time interval
      INP: (float)    below    = unflag data below this value (default end of the scan) 
           (float)    above    = unflag data above this value (default start of the scan)
           (int)      flag     = flag to be unset (default []: all flag values)
      """
      self.MessHand.debug('unflagInTime start...')
    
      # check channel list
      chanList = self.BolometerArray.checkChanList("all", flag='None')
      if len(chanList)<1: 
          self.MessHand.error("no valid channel")
          return

      # return Data index
      chanListIndexes = self.BolometerArray.getChanIndex(chanList) 

      # check dataType
      if ( dataType not in ['LST','lst'] and
           dataType not in ['MJD','mjd'] and
           dataType not in ['UT'] and
           dataType not in ['speed','Speed'] and
           dataType not in ['accel','Accel']):
          self.MessHand.error("use only this method to unflag in time")
          return

      # get the default values
      dataTest = self.ScanParam.get(dataType=dataType,flag='None')
      if above == '?':
          above = min(dataTest)
      if below == '?':
          below = max(dataTest)

      self.ScanParam.unflag(dataType=dataType,below=below,above=above,flag=flag)
      # Report the flag on the ScanParam attribute
      mask = self.ScanParam.FlagHandler.isUnsetMask()
      if len(nonzero(mask)) > 0:
          for chan in chanListIndexes:
              self.FlagHandler.unsetOnMask(mask, \
                                           self.rflags['INTEGRATION FLAGGED'], \
                                           dim=1, index=chan)
      
      self.MessHand.debug("... unflagInTime end")

  #--------------------------------------------------------------------------------
  def flagMJD(self, below='?', above='?', flag=8):
    """
    DES: Flag data in time interval
    INP: (float)    below    = flag data below this value (default end of the scan) 
         (float)    above    = flag data above this value (default start of the scan)
         (int)      flag     = flag to be set (default: 8 'TEMPORARY')
    """
    self.flagInTime(dataType='MJD', below=below, above=above, flag=flag)

  #--------------------------------------------------------------------------------
  def unflagMJD(self, below='?', above='?', flag=[]):
    """
    DES: Unflag data in time interval
    INP: (float)    below    = flag data below this value (default end of the scan) 
         (float)    above    = flag data above this value (default start of the scan)
         (int)      flag     = flag to be unset (default []: all flag values)
    """
    self.unflagInTime(dataType='MJD', below=below, above=above, flag=flag)

  #--------------------------------------------------------------------------------
  def flagPolygon(self, channel='all', system='EQ', poly=zeros((1,2)), inout='IN', flag=8):
      """
      DES: flag a position in the sky within or outside a given polygon
      INP: (int list) channel : list of channels to flag (default: 'all')
           (str)       system : coord. system, one of 'HO' (Az,El *OFFSETS*) or 
                                'EQ' (RA, Dec absolute coord.), default='EQ' 
           (float array) poly : vertices of polygon
           (str)        inout : inside/outside the polygon, one of 'IN' or 'OUT' 
           (int)         flag : flag to be set (default 8 'TEMPORARY')
      """
    
      self.MessHand.debug('flagPolygon start...')
    
      # check channel list
      chanList = self.BolometerArray.checkChanList(channel)
      if len(chanList)<1: 
          self.MessHand.error("no valid channel")
          return
      chanListIndexes = self.BolometerArray.getChanIndex(chanList)

      # check system
      system = string.upper(system)
      if system not in ['EQ','HO']: 
          self.MessHand.error("no valid coordinate system: "+system)
          return

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

      # check flags
      flag=self._removeReservedFlagValues(flag)
      if flag==None:
          self.MessHand.error("no valid flags")
          return       

      OffsetsUsed=array(self.BolometerArray.getChanSep(self.BolometerArray.UsedChannels))

      totalFlag = 0

      if system=='EQ':
          RADec = array([self.ScanParam.get('RA',flag='None'), \
                         self.ScanParam.get('Dec',flag='None')])
          OffsetsUsed = OffsetsUsed/3600.

          for chan in chanListIndexes:
              dRA = RADec[0]+(-1.*cos((3.1416/180.)*self.ScanParam.ParAngle)*OffsetsUsed[0,chan]\
                              +   sin((3.1416/180.)*self.ScanParam.ParAngle)*OffsetsUsed[1,chan])\
                             /array(cos(RADec[1]*3.1416/180.)) 
              dDec= RADec[1] + sin((3.1416/180.)*self.ScanParam.ParAngle)*OffsetsUsed[0,chan]\
                             + cos((3.1416/180.)*self.ScanParam.ParAngle)*OffsetsUsed[1,chan]
              mask = self.maskPolygon(dRA,dDec,poly,inout)
              if len(nonzero(mask)):
                  n0 = self.FlagHandler.nSet(flag, dim=1, index=chan)
                  self.FlagHandler.setOnMask(mask, flag, dim=1, index=chan)
                  n1 = self.FlagHandler.nSet(flag, dim=1, index=chan)
                  totalFlag += (n1-n0)

      elif system=='HO':
          AzEl = array([self.ScanParam.get('AzimuthOffset',flag='None'), \
                        self.ScanParam.get('ElevationOffset',flag='None')])

          for chan in chanListIndexes:
              dAz = AzEl[0] + OffsetsUsed[0,chan] 
              dEl = AzEl[1] + OffsetsUsed[1,chan]
              mask = self.maskPolygon(dAz,dEl,poly,inout)
              if len(nonzero(mask)):
                  n0 = self.FlagHandler.nSet(flag, dim=1, index=chan)
                  self.FlagHandler.setOnMask(mask, flag, dim=1, index=chan)
                  n1 = self.FlagHandler.nSet(flag, dim=1, index=chan)
                  totalFlag += (n1-n0)                        

      if totalFlag > 0:
          self.MessHand.info("%5i samples flagged with flag %s" % (totalFlag,str(flag)))
          self.__resetStatistics()
      else:
          self.MessHand.warning("Nothing flagged")
      
      self.MessHand.debug("... flagPolygon end")

  #--------------------------------------------------------------------------------
  def unflagPolygon(self, channel='all', system='EQ', poly=zeros((1,2)), inout='IN', flag=[]):
      """
      DES: unflag a position in the sky within or outside a given polygon
      INP: (int list) channel : list of channels to flag (default: 'all')
           (str)       system : coord. system, one of 'HO' (Az,El *OFFSETS*) or 
                                'EQ' (RA, Dec absolute coord.), default='EQ' 
           (float array) poly : vertices of polygon
           (str)        inout : inside/outside the polygon, one of 'IN' or 'OUT' 
           (int)         flag : flag to be set (default 8 'TEMPORARY')
      """
    
      self.MessHand.debug('unflagPolygon start...')
    
      # check channel list
      chanList = self.BolometerArray.checkChanList(channel)
      if len(chanList)<1: 
          self.MessHand.error("no valid channel")
          return
      chanListIndexes = self.BolometerArray.getChanIndex(chanList)

      # check system
      system = string.upper(system)
      if system not in ['EQ','HO']: 
          self.MessHand.error("no valid coordinate system: "+system)
          return

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

      # check flags
      flag=self._removeReservedFlagValues(flag)
      if flag==None:
          self.MessHand.error("no valid flags")
          return       

      OffsetsUsed=array(self.BolometerArray.getChanSep(self.BolometerArray.UsedChannels))

      totalFlag = 0

      if system=='EQ':
          RADec = array([self.ScanParam.get('RA',flag='None'), \
                         self.ScanParam.get('Dec',flag='None')])
          OffsetsUsed = OffsetsUsed/3600.

          for chan in chanListIndexes:
              dRA = RADec[0]+(-1.*cos((3.1416/180.)*self.ScanParam.ParAngle)*OffsetsUsed[0,chan]\
                              +   sin((3.1416/180.)*self.ScanParam.ParAngle)*OffsetsUsed[1,chan])\
                             /array(cos(RADec[1]*3.1416/180.)) 
              dDec= RADec[1] + sin((3.1416/180.)*self.ScanParam.ParAngle)*OffsetsUsed[0,chan]\
                             + cos((3.1416/180.)*self.ScanParam.ParAngle)*OffsetsUsed[1,chan]
              mask = self.maskPolygon(dRA,dDec,poly,inout)
              if len(nonzero(mask)):
                  n0 = self.FlagHandler.nUnset(flag, dim=1, index=chan)
                  self.FlagHandler.unsetOnMask(mask, flag, dim=1, index=chan)
                  n1 = self.FlagHandler.nUnset(flag, dim=1, index=chan)
                  totalFlag += (n1-n0)

      elif system=='HO':
          AzEl = array([self.ScanParam.get('AzimuthOffset',flag='None'), \
                        self.ScanParam.get('ElevationOffset',flag='None')])

          for chan in chanListIndexes:
              dAz = AzEl[0] + OffsetsUsed[0,chan] 
              dEl = AzEl[1] + OffsetsUsed[1,chan]
              mask = self.maskPolygon(dAz,dEl,poly,inout)
              if len(nonzero(mask)):
                  n0 = self.FlagHandler.nUnset(flag, dim=1, index=chan)
                  self.FlagHandler.unsetOnMask(mask, flag, dim=1, index=chan)
                  n1 = self.FlagHandler.nUnset(flag, dim=1, index=chan)
                  totalFlag += (n1-n0)                        

      if totalFlag > 0:
          self.MessHand.info("%5i samples unflagged with flag %s" % (totalFlag,str(flag)))
          self.__resetStatistics()
      else:
          self.MessHand.warning("Nothing unflagged")
      
      self.MessHand.debug("... unflagPolygon end")

  #--------------------------------------------------------------------------------
  def maskPolygon(self,x,y,poly,inout='IN'):
      """
      DES: create an array of zeros and ones for a list of points being inside/outside 
             a polygon
      INP: (float array)  x/y : coordinates of points 
           (float array) poly : vertices of polygon
           (str)        inout : inside/outside the polygon
      OUT: (float array)      : array to be used for masking data points
      """
      mask=[]
      if inout=='IN':            # look for points inside the polygon 
        for i in range(len(x)):
          mask.append(inPolygon(x[i],y[i],poly))
      elif inout=='OUT':         # look for points outside the polygon 
        for i in range(len(x)):
          mask.append(outPolygon(x[i],y[i],poly))
      else:                      # wrong or missing information 
        self.MessHand.warning("Inside or outside the polygon?")
      mask=array(mask)
      return mask

  #--------------------------------------------------------------------------------
  def flagPosition(self, channel='all', Az = 0, El = 0, radius = 0, flag = 8,
                   offset=1, outer=0, relative=1):
      """
      DES: flag a position in the sky within a given radius
      INP: (int list) channel : list of channel to flag (default: 'all')
           (float)      Az/El : the horizontal reference position (arcsec for offsets, deg for absolute) 
           (float)     radius : aperture to flag in unit of the reference position
           (int)         flag : flag to be set (default 8 'TEMPORARY')
           (logical)   offset : flag on the offsets (default yes,)
           (logical)   outer  : flag OUTSIDE the given radius? default: no
           (logical) relative : use bolometer offsets w.r.t. to reference channel
                                (relative=1, default) or use absolute offsets (relative=0)
      """
    
      self.MessHand.debug('flagPosition start...')
    
      # check channel list
      chanList = self.BolometerArray.checkChanList(channel)
      if len(chanList)<1: 
          self.MessHand.error("no valid channel")
          return

      chanListIndexes = self.BolometerArray.getChanIndex(chanList)

      if offset:
          AzEl            = array([self.ScanParam.get('AzimuthOffset',flag='None'), \
                                   self.ScanParam.get('ElevationOffset',flag='None')])
          OffsetsUsed     = array(self.BolometerArray.getChanSep(self.BolometerArray.UsedChannels))
      else:
          AzEl            = array([self.ScanParam.get('Azimuth',flag='None'), \
                                   self.ScanParam.get('ElevationOffset',flag='None')])
          OffsetsUsed     = array(self.BolometerArray.getChanSep(self.BolometerArray.UsedChannels))/3600.

      if not relative:
          OffsetsUsed = self.BolometerArray.Offsets

      flag = self._removeReservedFlagValues(flag)
      if flag==None:
          self.MessHand.error("no valid flags")
          return

      totalFlag = 0
      for chan in chanListIndexes:
          dAz = AzEl[0] + OffsetsUsed[0,chan] - Az 
          dEl = AzEl[1] + OffsetsUsed[1,chan] - El
          if outer:
              mask = where((dAz**2 + dEl**2) > radius**2, 1, 0)
          else:
              mask = where((dAz**2 + dEl**2) <= radius**2, 1, 0)
          if len(nonzero(mask)):
              n0 = self.FlagHandler.nSet(flag, dim=1, index=chan)
              self.FlagHandler.setOnMask(mask, flag, dim=1, index=chan)
              n1 = self.FlagHandler.nSet(flag, dim=1, index=chan)
              totalFlag += (n1-n0)

      if totalFlag > 0:
          self.MessHand.info("%5i samples flagged with flag %s" % (totalFlag,str(flag)))
          self.__resetStatistics()
      else:
          self.MessHand.warning("Nothing flagged")
      
      self.MessHand.debug("... flagPosition end")


  #--------------------------------------------------------------------------------
  def flagRadius(self, channel='all', radius = 0, flag = 8, outer=0):
      """
      DES: flag time series (all channels) by reference offset in Az/El
      INP: (int list) channel : list of channel to flag (default: 'all')
           (float)     radius : aperture to flag in ARCSECONDS
           (int)         flag : flag to be set (default 8 'TEMPORARY')
           (logical)   outer  : flag OUTSIDE the given radius? default: no
      """
    
      self.MessHand.debug('flagRadius start...')

      chanList = self.BolometerArray.checkChanList(channel)
      if len(chanList)<1: 
          self.MessHand.error("no valid channel")
          return

      chanListIndices = self.BolometerArray.getChanIndex(chanList)
    
      AzEl            = array([self.ScanParam.get('AzimuthOffset',flag='None'), \
                               self.ScanParam.get('ElevationOffset',flag='None')])
      
      flag = self._removeReservedFlagValues(flag)
      if flag==None:
          self.MessHand.error("no valid flags")
          return

      
      if outer:
          mask = where((AzEl[0,::]**2 + AzEl[1,::]**2) > radius**2, 1, 0)
      else:
          mask = where((AzEl[0,::]**2 + AzEl[1,::]**2) <= radius**2, 1, 0)
      if len(nonzero(mask)):
                    
          self.ScanParam.FlagHandler.setOnMask(mask,iFlags=flag)
      
          for chan in chanListIndices:
              self.FlagHandler.setOnMask(mask, self.rflags['INTEGRATION FLAGGED'], \
                                         dim=1, index=chan)
              
          nflagged=shape(compress(mask,mask))[0]
          self.MessHand.info("%5i timestamps flagged with flag %s" % (nflagged, str(flag)))
          self.__resetStatistics()

      else:
          self.MessHand.warning("Nothing flagged")
      
      self.MessHand.debug("... flagPosition end")

  #--------------------------------------------------------------------------------
  def unflagPosition(self, channel='all', Az = 0, El = 0, radius = 0, flag = [], offset=1):
      """
      DES: unflag a position in the sky within a given radius
      INP: (int list) channel : list of channel to unflag (default: 'all')
           (float)      Az/El : the horizontal reference position (arcsec for offsets, deg for absolute) 
           (float)     radius : aperture to unflag in unit of the reference position
           (int)         flag : unflag to be set (default []: unflag all non-reserved flag values)
           (logical)   offset : unflag on the offsets (default yes,) 
      """
    
      self.MessHand.debug('unflagPosition start...')
    
      # check channel list
      chanList = self.BolometerArray.checkChanList(channel)
      if len(chanList)<1: 
          self.MessHand.error("no valid channel")
          return

      chanListIndexes = self.BolometerArray.getChanIndex(chanList)

      if offset:
          AzEl            = array([self.ScanParam.get('AzimuthOffset',flag='None'), \
                                   self.ScanParam.get('ElevationOffset',flag='None')])
          OffsetsUsed     = array(self.BolometerArray.getChanSep(self.BolometerArray.UsedChannels))
      else:
          AzEl            = array([self.ScanParam.get('Azimuth',flag='None'), \
                                   self.ScanParam.get('ElevationOffset',flag='None')])
          OffsetsUsed     = array(self.BolometerArray.getChanSep(self.BolometerArray.UsedChannels))/3600.
 
      flag = self._removeReservedFlagValues(flag)
      if flag==None:
          self.MessHand.error("no valid flags")
          return

      totalFlag = 0
      for chan in chanListIndexes:
          dAz = AzEl[0] + OffsetsUsed[0,chan] - Az 
          dEl = AzEl[1] + OffsetsUsed[1,chan] - El
          mask = where((dAz**2 + dEl**2) <= radius**2, 1, 0)
          if len(nonzero(mask)):
              n0 = self.FlagHandler.nUnset(flag, dim=1, index=chan)
              self.FlagHandler.unsetOnMask(mask, dim=1, index=chan)
              n1 = self.FlagHandler.nUnset(flag, dim=1, index=chan)
              totalFlag += (n1-n0)

      if totalFlag > 0:
          self.MessHand.info("%5i samples unflagged with flag %s" % (totalFlag,str(flag)))
          self.__resetStatistics()
      else:
          self.MessHand.warning("Nothing unflagged")
      
      self.MessHand.debug("... unflagPosition end")


  #--------------------------------------------------------------------------------
  def flagSubscan(self,subList,flag=7):
    """
    DES: flag subscans
    INP: (int list) subList = list of subscan numbers (or single number)
                              to be flagged
              (int) flag    = value of flags to set (default: 7 'SUBSCAN FLAGGED')
    """
    # If a single number given, convert to list
    if type(subList) == type(1):
        subList = [subList]
    # check if subscan numbers exist
    for subNum in subList:
        if subNum not in self.ScanParam.SubscanNum:
            subList.remove(subNum)
    if len(subList) == 0:
        self.MessHand.warning("No valid subscan number given")
        return

    numberFlags = 0
    # Use MJD rather than LST: it's never -999!
    mjd = self.getChanData('mjd',flag='None')  # get all mjd values, including
                        # flagged datapoints (same convention in SubscanIndex)
    for subNum in subList:
        numIndex = self.ScanParam.SubscanNum.index(subNum)
        mjd1 = mjd[self.ScanParam.SubscanIndex[0,numIndex]]
        mjd2 = mjd[self.ScanParam.SubscanIndex[1,numIndex]-1]
        self.flagInTime('mjd',above=mjd1,below=mjd2,flag=flag)
        numberFlags += 1

    self.MessHand.info(str("%i subscans flagged with flag %s"%(numberFlags,str(flag))))

  #--------------------------------------------------------------------------------
  def unflagSubscan(self,subList,flag=[]):
    """
    DES: unflag subscans
    INP: (int list) subList = list of subscan numbers (or single number)
                              to be unflagged
              (int) flag    = value of flags to unset (default []: all flag values)
    """
    # If a single number given, convert to list
    if type(subList) == type(1):
        subList = [subList]
    # check if subscan numbers exist
    for subNum in subList:
        if subNum not in self.ScanParam.SubscanNum:
            subList.remove(subNum)
    if len(subList) == 0:
        self.MessHand.warning("No valid subscan number given")
        return

    numberFlags = 0
    # Use MJD rather than LST: it's never -999!
    mjd = self.getChanData('mjd',flag='None')  # get all mjd values, including
                        # flagged datapoints (same convention in SubscanIndex)
    for subNum in subList:
        numIndex = self.ScanParam.SubscanNum.index(subNum)
        mjd1 = mjd[self.ScanParam.SubscanIndex[0,numIndex]]
        mjd2 = mjd[self.ScanParam.SubscanIndex[1,numIndex]-1]
        self.unflagInTime('mjd',above=mjd1,below=mjd2,flag=flag)
        numberFlags += 1

    self.MessHand.info(str("%i subscan unflagged with flag %s"%(numberFlags,str(flag))))


  #--------------------------------------------------------------------------------
  def flagLon(self, channel='all', below='?', above='?', flag=8):
    """
    NAM: flagLon (method)  
    DES: Flag data in Longitude interval
    INP: (int list) channel = list of channel to flag (default: all)
         (float)    below   = flag data below this value
         (float)    above   = flag data above this value
         (int)      flag    = flag to be set (default 8 'TEMPORARY')
    """
    self.flag(dataType='azoff', channel=channel, below=below, above=above, flag=flag)

  #--------------------------------------------------------------------------------
  def unflagLon(self, channel='all', below='?', above='?', flag=[]):
    """
    NAM: unflagLon (method)  
    DES: Unflag data in Longitude interval
    INP: (int list) channel = list of channel to flag (default: all)
         (float)    below   = flag data below this value
         (float)    above   = flag data above this value
         (int)      flag    = flag to be unset (default []: all non-reserved flag values)
    """
    self.flag(dataType='azoff', channel=channel, below=below, above=above, flag=flag)

  #--------------------------------------------------------------------------------

  def flagTurnaround(self, flag=1):
    """
    NAM: flagTurnaround (method)  
    DES: flag subscans where azimuth offset changes sign
    INP: (int)      flag    = flag to be set (default 1 'TURNAROUND')
    """
    try:
        sublist=[]
        for i in range(shape(self.ScanParam.SubscanPos)[0]):
            if (self.ScanParam.SubscanPos[i] != 0):
                sublist.append(self.ScanParam.SubscanNum[i])
                
        self.flagSubscan(subList=sublist,flag=flag)
    except:
        self.MessHand.warning('No subscan information, nothing will be done')

  #--------------------------------------------------------------------------------

  def unflagTurnaround(self, flag=[]):
    """
    NAM: unflagTurnaround (method)  
    DES: unflag subscans where azimuth offset changes sign
    INP: (int)      flag    = flag to be unset (default []: all flag values)
    """
    try:
        sublist=[]
        for i in range(shape(self.ScanParam.SubscanPos)[0]):
            if (self.ScanParam.SubscanPos[i] != 0):
                sublist.append(self.ScanParam.SubscanNum[i])
                
        self.unflagSubscan(subList=sublist,flag=flag)
    except:
        self.MessHand.warning('No subscan information, nothing will be done')

  #--------------------------------------------------------------------------------
  def flagSparseSubscans(self,minLiveFrac=0.3):
      """
      DES: flag whole subscans with few live time stamps
      INP: (f)   minLiveFrac : minimum fraction of live time stamps
      """

      self.__statistics()
      subnum=self.ScanParam.SubscanNum
      subin=self.ScanParam.SubscanIndex
      sflags=[]
      tflags=self.ScanParam.FlagHandler.getFlags()
      for s in range(len(subnum)):
          n_in_sub=subin[1,s]-subin[0,s]
          flags_in_sub=tflags[subin[0,s]:subin[1,s]]
          mask=where(flags_in_sub == 0,1,0)
          ok_in_sub=shape(compress(mask,mask))[0]
          if (float(ok_in_sub) > 0):
              if (float(ok_in_sub)/float(n_in_sub) < minLiveFrac):
                  sflags.append(subnum[s])

      if sflags:
          self.flagSubscan(sflags)

          self.__resetStatistics()

      return shape(sflags)[0]

  

  def flagSubscanByRms(self,above=2.,maxIter=20):
      """
      DES: iteratively flag subscans with high rms. Subscan rms is
      determined as the mean of all channels.
      INP: (f)   above : flag data with value > 'above'*rms
      (i) maxIter : maximum number of iterations
      """
      
      if not self.__statisticsDone:
          self.__statistics()
          
      subnum=self.ScanParam.SubscanNum
      subflags=subnum*0
      subrms=array(self.getChanListData('rms_s'))
      mean_subrms=subrms[0,::]*0.0
      for i in range(len(mean_subrms)):
          mean_subrms[i]=fStat.f_mean(subrms[:,i])       
      computation_mask=where(mean_subrms > 0,1,0)
      done=0
      counter=0
      while (done == 0):
          if (counter > maxIter):
              done=1
          rms_considered=compress(computation_mask,mean_subrms)
          meanRms=fStat.f_median(rms_considered)
          rmsRms=fStat.f_rms(rms_considered,meanRms)     
          subflagmask=where(mean_subrms > meanRms+above*rmsRms, 1,0)
          putmask(mean_subrms,subflagmask,-1)
          new_comp_mask=where(mean_subrms > 0,1,0)
          if (new_comp_mask.tolist() == computation_mask.tolist()):
              done=1
          else:
              computation_mask = new_comp_mask
              counter+=1      
      totmask=where(mean_subrms < 0,1,0)
      subflaglist=compress(totmask,subnum)
      totmask2=where(mean_subrms > 0,1,0)

      if subflaglist:
          self.flagSubscan(subflaglist)
          self.__resetStatistics()

      nflagged=shape(compress(totmask,totmask))[0]
      

      return nflagged

    #--------------------------------------------------------------------------------

  
  #--------------------------------------------------------------------------------
      
  
  def flagSpeed(self, below='?', above='?', flag=3):
    """
    DES: Flag data according to telescope speed
    INP: (float)    below   = flag data below this value
         (float)    above   = flag data above this value
         (int)      flag    = flag to be set (default 3 'ELEVATION VELOCITY THRESHOLD')
    """
    self.flagInTime(dataType='speed', below=below, above=above, flag=flag)

  #--------------------------------------------------------------------------------
  
  def unflagSpeed(self, below='?', above='?', flag=[]):
    """
    DES: Unflag data according to telescope speed
    INP: (float)    below   = unflag data below this value
         (float)    above   = unflag data above this value
         (int)      flag    = flag to be unset (default []: all flag values)
    """
    self.unflagInTime(dataType='speed', below=below, above=above, flag=flag)

  #--------------------------------------------------------------------------------
  def flagAccel(self, channel='all', below='?', above='?', flag=2):
    """
    DES: Flag data according to telescope acceleration
    INP: (float)    below   = flag data below this value
         (float)    above   = flag data above this value
         (int)      flag    = flag to be set (default 2 'ACCELERATION THRESHOLD')
    """
    self.flagInTime(dataType='accel', below=below, above=above, flag=flag)

  #--------------------------------------------------------------------------------
  def unflagAccel(self, channel='all', below='?', above='?', flag=[]):
    """
    DES: Unflag data according to telescope acceleration
    INP: (float)    below   = unflag data below this value
         (float)    above   = unflag data above this value
         (int)      flag    = flag to be unset (default []: all flag values)
    """
    self.unflagInTime(dataType='accel', below=below, above=above, flag=flag)


  #--------------------------------------------------------------------------------
  def flagNan(self, channel='all', flag=1):

      # check channel list
      chanList = self.BolometerArray.checkChanList(channel)
      if len(chanList)<1: 
          self.MessHand.error("no valid channel")
          return

      chanListIndexes = self.BolometerArray.getChanIndex(chanList)
      nInt = self.ScanParam.NInt
      for c in chanListIndexes:
          mask = fUtilities.masknan(self.Data[:,c])
          if len(nonzero(mask)):
              self.FlagHandler.setOnMask(mask, flag, dim=1, index=c)

  #--------------------------------------------------------------------------------
  # FFT filtering methods
  #--------------------------------------------------------------------------------
  def blankFreq(self, channel='all', below='?', above='?'):
      """
      DES: Permanently remove some frequency interval in the Fourrier spectrum
           of the signal. This is computed subscan by subscan.
      INP: (int list) channel = list of channel to flagprocess (default: all)
           (float)    below   = filter data below this value
           (float)    above   = filter data above this value
      """
      self.MessHand.debug('blankFreq start...')
    
      # check channel list
      chanList = self.BolometerArray.checkChanList(channel)
    
      # We use only MJD values for fft. They are the correct time stamps.
      mjd = self.getChanData('mjd',chanList,flag='None')
      nbSub = len(self.ScanParam.SubscanNum)   # number of subscans

      for chan in chanList:
          flux = self.getChanData('flux',chan,flag='None')
          num = self.BolometerArray.getChanIndex(chan)[0]
          for i in range(nbSub):
              ind1 = self.ScanParam.SubscanIndex[0,i]
              ind2 = self.ScanParam.SubscanIndex[1,i]
              theTime = mjd[ind1:ind2]
              theFlux = flux[ind1:ind2]
              if not theFlux.iscontiguous():
                  theFlux = theFlux.copy()
              # select only unflagged data - zeros for flagged ones
              maskOk = nonzero(self.FlagHandler.isUnsetMask(dim=1, index=num)[ind1:ind2])
              theFluxOk = zeros(shape(theFlux),'f')
              put(theFluxOk,maskOk,take(theFlux,maskOk))

              oneFFT = FilterFFT(theTime,theFluxOk)
              oneFFT.doFFT()
              # Do the filtering
              oneFFT.blankAmplitude(above=above,below=below)
              oneFFT.invFFT()  # this also "unbins" the data

              put(theFlux,maskOk,take(oneFFT.Y,maskOk))  # flagged values unchanged
              # Now update flagged values with inverse FFT + previous value
              maskBad = nonzero(self.FlagHandler.isSetMask(dim=1,index=num)[ind1:ind2])
              oldVal = take(theFlux,maskBad)
              put(theFlux,maskBad,take(oneFFT.Y,maskBad)+oldVal)
              self.Data[ind1:ind2,num] = theFlux

      self.__resetStatistics()

  #--------------------------------------------------------------------------------
  def reduceFreq(self, channel='all', center=50., width=1., factor=10.,optimize=1, window=4):
      """
      DES: Permanently reduce some frequency interval in the Fourrier spectrum
           of the signal. This is computed subscan by subscan.
      INP: (int list) channel = list of channel to process (default: all)
                   (f) center = central frequency, in Hz
                   (f)  width = line FWHM
                   (f) factor = attenuation factor
      """
      self.MessHand.debug('reduceFreq start...')
    
      # check channel list
      chanList = self.BolometerArray.checkChanList(channel)
    
      # We use only MJD values for fft. They are the correct time stamps.
      mjd = self.getChanData('mjd',chanList,flag='None')
      nbSub = len(self.ScanParam.SubscanNum)   # number of subscans

      for chan in chanList:
          flux = self.getChanData('flux',chan,flag='None')
          num = self.BolometerArray.getChanIndex(chan)[0]
          for i in range(nbSub):
              ind1 = self.ScanParam.SubscanIndex[0,i]
              ind2 = self.ScanParam.SubscanIndex[1,i]
              theTime = mjd[ind1:ind2]
              theFlux = flux[ind1:ind2]
              if not theFlux.iscontiguous():
                  theFlux = theFlux.copy()
              # select only unflagged data - zeros for flagged ones
              maskOk = nonzero(self.FlagHandler.isUnsetMask(dim=1, index=num)[ind1:ind2])
              theFluxOk = zeros(shape(theFlux),'f')
              put(theFluxOk,maskOk,take(theFlux,maskOk))

              oneFFT = FilterFFT(theTime,theFluxOk)
              oneFFT.doFFT(windowing=window)
              # Do the filtering
              oneFFT.reduceAmplitude(center=center,width=width,factor=factor)
              oneFFT.invFFT(windowing=window)
              put(theFlux,maskOk,take(oneFFT.Y,maskOk))  # flagged values unchanged
              # Now update flagged values with inverse FFT + previous value
              maskBad = nonzero(self.FlagHandler.isSetMask(dim=1,index=num)[ind1:ind2])
              oldVal = take(theFlux,maskBad)
              put(theFlux,maskBad,take(oneFFT.Y,maskBad)+oldVal)
              self.Data[ind1:ind2,num] = theFlux

      self.__resetStatistics()

  #--------------------------------------------------------------------------------
  def taperFreq(self, channel='all', above='?', N=2, window=4):
      """
      DES: Permanently taper off Fourier spectrum above given value 
           of the signal
      INP: (int list) channel = list of channel to flagprocess (default: all)
           (float)    above   = filter data above this value
           (int)      N       = Butterworth steepenes order
      """
      self.MessHand.debug('taperFreq start...')
    
      # check channel list
      chanList = self.BolometerArray.checkChanList(channel)
    
      # We use only MJD values for fft. They are the correct time stamps.
      mjd   = self.getChanData('mjd',chanList,flag='None')
      nbSub = len(self.ScanParam.SubscanNum)   # number of subscans
    
      for chan in chanList:
          flux = self.getChanData('flux',chan,flag='None')
          num = self.BolometerArray.getChanIndex(chan)[0]

          for i in range(nbSub):
              ind1 = self.ScanParam.SubscanIndex[0,i]
              ind2 = self.ScanParam.SubscanIndex[1,i]
              theTime = mjd[ind1:ind2]
              theFlux = flux[ind1:ind2]
              if not theFlux.iscontiguous():
                  theFlux = theFlux.copy()

              # select only unflagged data - zeros for flagged ones
              maskOk = nonzero(self.FlagHandler.isUnsetMask(dim=1, index=num)[ind1:ind2])
              theFluxOk = zeros(shape(theFlux),'f')
              put(theFluxOk,maskOk,take(theFlux,maskOk))

              oneFFT = FilterFFT(theTime,theFluxOk)
              oneFFT.doFFT(windowing=window)
              oneFFT.taperAmplitude(above=above,N=N)
              
              oneFFT.invFFT(windowing=window)
              put(theFlux,maskOk,take(oneFFT.Y,maskOk))  # flagged values unchanged
              # Now update flagged values with inverse FFT + previous value
              maskBad = nonzero(self.FlagHandler.isSetMask(dim=1,index=num)[ind1:ind2])
              oldVal = take(theFlux,maskBad)
              put(theFlux,maskBad,take(oneFFT.Y,maskBad)+oldVal)

              self.Data[ind1:ind2,num] = theFlux

      self.__resetStatistics()

  #--------------------------------------------------------------------------------
  def flattenFreq(self, channel='all', below=0.1, hiref=1.,optimize=1, window=4):
      """
      DES: flatten the 1/F part of the FFT using constant amplitude
      INP: (int list) channel = list of channels to process (default: all)
           (float)    below   = filter data below this value
           (float)    hiref   = amplitudes at f < below will be replaced with
                                the average value between below and hiref
      """
      self.MessHand.debug('flattenFreq start...')
    
      # check channel list
      chanList = self.BolometerArray.checkChanList(channel)
    
      # We use only MJD values for fft. They are the correct time stamps.
      mjd = self.getChanData('mjd',chanList,flag='None')
      nbSub = len(self.ScanParam.SubscanNum)   # number of subscans

      for chan in chanList:
          flux = self.getChanData('flux',chan,flag='None')
          num = self.BolometerArray.getChanIndex(chan)[0]

          for i in range(nbSub):
              ind1 = self.ScanParam.SubscanIndex[0,i]
              ind2 = self.ScanParam.SubscanIndex[1,i]
              theTime = mjd[ind1:ind2]
              theFlux = flux[ind1:ind2]
              if not theFlux.iscontiguous():
                  theFlux = theFlux.copy()

              # select only unflagged data - zeros for flagged ones
              maskOk = nonzero(self.FlagHandler.isUnsetMask(dim=1, index=num)[ind1:ind2])
              theFluxOk = zeros(shape(theFlux),'f')
              put(theFluxOk,maskOk,take(theFlux,maskOk))
                  
              oneFFT = FilterFFT(theTime,theFluxOk)
              oneFFT.doFFT(windowing=window)
              # compute median value between below and hiref
              mask = nonzero(greater(oneFFT.Freq,below) and less(oneFFT.Freq,hiref))
              inside = take(oneFFT.Amplitude,mask)
              meanAmp = fStat.f_median(inside)
              # replace amplitudes at lower freq with this value
              mask = nonzero(less(oneFFT.Freq,below))
              for k in mask:
                  oneFFT.Amplitude[k] = meanAmp
              # inverse FFT
              oneFFT.invFFT(windowing=window)
              put(theFlux,maskOk,take(oneFFT.Y,maskOk))  # flagged values unchanged
              # Now update flagged values with inverse FFT + previous value
              maskTmp = self.FlagHandler.isSetMask(dim=1,index=num)[ind1:ind2]
              # ... but only if not NaN
              maskNan = fUtilities.masknan(self.Data[ind1:ind2,num])
              maskBad = nonzero(logical_and(maskTmp,logical_not(maskNan)))
              oldVal = take(theFlux,maskBad)
              put(theFlux,maskBad,take(oneFFT.Y,maskBad)+oldVal)

              self.Data[ind1:ind2,num] = theFlux   # store results in data

      self.__resetStatistics()

  #--------------------------------------------------------------------------------
  #----- statistics methods -------------------------------------------------------
  #--------------------------------------------------------------------------------

  def __resetStatistics(self):
      """
      DES: to be called every time the data are altered:
      statistics and correlation matrix should be recomputed
      """
      self.__statisticsDone = 0
      self.__corMatrixDone  = 0
      self.__pcaDone = 0

  def __statistics(self):
    """
    NAM: statistics (method)  
    DES: computes mean, median, rms for all scans and subscans for all used channels
    """

    self.MessHand.debug('statistics start...')
    
    myTiming = Timing()

    if self._existData(): 

      myTiming.setTime()

      Mean, Med, SDev, MDev = fStat.arraystat(self.Data,self.FlagHandler.getFlags())
      self.ChanMean = Mean
      self.ChanMed  = Med
      self.ChanRms  = SDev

      Mean_s, Med_s, SDev_s, MDev_s = fStat.arraystat_s(self.Data,self.FlagHandler.getFlags(),\
                                                        self.ScanParam.SubscanIndex)
      self.ChanMean_s = Mean_s
      self.ChanMed_s  = Med_s
      self.ChanRms_s  = SDev_s

      # Mark Statistics as done
      self.__statisticsDone = 1

      self.MessHand.debug(" statistics by scan in "+ str(myTiming))

    self.MessHand.debug('... statistics end')


  def slidingRms(self,nbInteg=10, channel=[], flag=[], getFlagged=0):
    """
    NAM: slidingRms (method)
    DES: compute rms in a sliding window
    INP: (int)    nbInteg : number of elements on which one rms is computed (=window size)
         (i list) channel : list of channel to flag (default: all; [] : current list)
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
    OUT: (array) the rms are returned
    """
    # check channel list
    chanList = self.BolometerArray.checkChanList(channel)
    if len(chanList)<1:
      self.MessHand.error("no valid channel")
      return

    result = []
    for chan in chanList:
      # get data for this chan - that's a 1D array
      chanData = self.getChanData('flux',chan=chan,flag=flag,getFlagged=getFlagged)
      lenData = shape(chanData)[0]
      # build a 2D array corresponding to the sliding window
      # there will be nb_data - window_size possible windows
      slidingData = zeros((nbInteg,lenData-nbInteg),Float32)
      for i in range(lenData-nbInteg):
        slidingData[:,i] = chanData[i:i+nbInteg]
      # Now call arraystat on this 2D array, with all flags set to zero
      # (filtering already done in call to getChanData)
      Mean, Med, SDev, MDev = fStat.arraystat(slidingData,0*slidingData)
      # store Rms in result
      result.append(SDev)

    # Return result as an array, with channels along 2nd dimension
    return transpose(array(result))


      
  #--------------------------------------------------------------------------------
  # Correlated Noise methods
  #-----------------------------------------------------------------------
  def medianCorrelations(self,chanList=[],numCorr=0):

    """
    DES: returns the  median correlation of each channel with all other channels
    INP: (i list) chanList : the list of channels to consider
             (int) numCorr : if set to non-zero, takes the median correlation of the 
                             numCorr most correlated channels   
    """

    chanList = self.BolometerArray.checkChanList(chanList)
    chanListIndexes = self.BolometerArray.getChanIndex(chanList)
      
    if not chanList:
      self.MessHand.error('no valid channel')

    if not self.__corMatrixDone:
        self.computeWeights(chanList)

    matrix = compress2d(self.CorMatrix,chanListIndexes)

    corrs=matrix[0,::]

    for rownum in range(len(matrix[:,0])):
        row=matrix[rownum,::]
        mask=where(row < 1.0,1,0)
        row=compress(mask,row)
        if numCorr:
            maxCorrs=tolist_boa(take(row,argsort(row)[-(abs(numCorr))::]))
            corrs[rownum]=fStat.f_median(maxCorrs)
        else:
            corrs[rownum]=fStat.f_median(row)

    return corrs
    



  def plotCorMatrix(self,chanList=[], check = 1, distance=0, weights=0, xLabel='Channels', style='idl4', limitsZ=[]):
    """
    DES: plot the correlation matrix
    INP: (i list) chanList : the list of channel to plot
         (l)         check : check the chanList first ( default : yes)
         (l)      distance : sort the second dimension by distance (default : no)
         (l)       weights : plot weights instead of correlation matrix (default: no)
    """


    if not self.__corMatrixDone:
        self.computeWeights()

    if check:
      chanList = self.BolometerArray.checkChanList(chanList)
      
    if not chanList:
      self.MessHand.error('no channel to plot')

    chanListIndexes = self.BolometerArray.getChanIndex(chanList)
    
    if not weights:
        matrix = compress2d(self.CorMatrix,chanListIndexes)
        subCaption = ' correlation matrix'
    else:
        matrix = compress2d(self.Weight,chanListIndexes)
        subCaption = ' weight matrix'
        
    yLabel = xLabel

    if distance:
      yLabel ='Distance to channel'

      nBolo = len(chanList)
      subChanSep = compress2d(self.BolometerArray.ChannelSep,\
                              chanList-1)

      for i in range(nBolo):
        indexArray = arange(nBolo)
        sortedIndex = take(indexArray,(argsort(subChanSep[i,::])))
        matrix[i,::] = take(matrix[i,::],(sortedIndex))

    
    Plot.draw(matrix, \
              labelX=xLabel,labelY=yLabel, \
              caption=self.ScanParam.caption()+subCaption, \
              style=style, wedge=1,\
              limitsZ=limitsZ, nan=1)

    
    #-----------------------------------------------------------------------
  def plotCorDist(self,chanList=[],average=1,upperlim=-1.,check=1,style='p',ci=1,overplot=0,limitsX=[],limitsY=[],pointsize=3.,plot=1):
    """
    DES: plot correlations (correlation matrix) as a function
         of channel separation
    INP: (i list) chanList : the list of channels to plot
         (i)       average : number of data to average over (for easier viewing)
         (f)      upperlim : return only distances in arcsec below this value
                             (negative value means no limit, which is the default)
         (l)         check : check the chanList first (default: yes)
         (l)          plot : actually produce a plot? (defult: yes)
    """

    if check:
      chanList = self.BolometerArray.checkChanList(chanList)
      
    if not chanList:
      self.MessHand.error('no channel to plot')

    chanListIndexes = self.BolometerArray.getChanIndex(chanList)

    self.__corMatrixDone=0
    corMatrix = fSNF.cmatrix(self.Data,self.FlagHandler.getFlags(),chanListIndexes)
    self.CorMatrix        = corMatrix
    self.__corMatrixDone  = 1
    
    subChanSep = compress2d(self.BolometerArray.ChannelSep,chanListIndexes)
    nBolo = len(chanList)
    matrix = compress2d(self.CorMatrix,chanListIndexes)
    for i in range(nBolo):
        indexArray = arange(nBolo)
        sortedIndex = take(indexArray,(argsort(subChanSep[i,::])))
        matrix[i,::] = take(matrix[i,::],(sortedIndex))
        subChanSep[i,::] = take(subChanSep[i,::],(sortedIndex))

    dataX=take(ravel(subChanSep),(argsort(ravel(subChanSep))))
    dataY=take(ravel(matrix),(argsort(ravel(subChanSep))))

    if (upperlim > 0):
        mask=where((dataX < upperlim),1,0)
        dataX=compress(mask,dataX)
        dataY=compress(mask,dataY)

    if (average > 1):
        average=int(average)
        newshape=len(dataX)/average
        dataX=sum(resize(dataX,(newshape,average)),1)/float(average)
        dataY=sum(resize(dataY,(newshape,average)),1)/float(average)
        
    if plot:
        BogliConfig.point['size']=pointsize
        Plot.plot(dataX,dataY,overplot=overplot,ci=ci,style=style,\
                  limitsX=limitsX,limitsY=limitsY,\
                  labelX='Channel separation (arcsec)',labelY='correlation',\
                  caption=self.ScanParam.caption())
        BogliConfig.point['size']=0.01
    else:
        return dataX, dataY
    

  #-----------------------------------------------------------------------
  def computeWeights(self, chanList=[], minCorr=0., a=0.95, b=2.0, core=10., beta=2.):
    """
    DES: compute correlation and weight matrix of the used channels
         Weight is a non-linear rescaling of the correlation coefficient
         
                       weight_nm = ( CM_nm - a * min_m( CM_nm ) )**b
                       
         an additionnal weighting factor is applied with channel separation

                       weight_nm  = weight_nm * 1.0 / ( 1 + ( dist_nm / core )**beta )
         
    INP: (i list)    : chanList restrict the computation to certain channel (default : all used channel)
         (f) minCorr : minimum correlation coefficient (defaut:0, should be positiv)
         (f) a       : parameter for weights, usually = 0.90-0.98
         (f) b       : parameter for weights, usually = 1
         (f) core    : core radius in arcmin for radial weighting (weight = 0.5)
         (f) beta    : beta for beta profile for radial weighting
    """


    myTiming      = Timing()
    corMatrix     = self.CorMatrix
    nUsedChannels = self.BolometerArray.NUsedChannels

    if not chanList :
        chanList  = self.BolometerArray.UsedChannels

    chanList        = self.BolometerArray.checkChanList(chanList)
    chanListIndexes = self.BolometerArray.getChanIndex(chanList)
    #channelFlags    = self.BolometerArray.FlagHandler.getFlags()

    # Compute the Correlation Matrix first
    if not self.__corMatrixDone:
        corMatrix = fSNF.cmatrix(self.Data,self.FlagHandler.getFlags(),chanListIndexes)
        # Set the corMatrixDone
        self.CorMatrix        = corMatrix
        self.__corMatrixDone  = 1
        self.MessHand.debug("corMatrix computed in "+str(myTiming))

    myTiming.setTime()

    # Compute the Weights now
    boloWeight = zeros((nUsedChannels,nUsedChannels),Float32)
    chanSep    = compress2d(self.BolometerArray.ChannelSep,self.BolometerArray.UsedChannels-1)
    
    boloWeight = fSNF.wmatrix(corMatrix,chanSep,chanListIndexes ,minCorr,a,b,core,beta)

    # no need to normalize since in SNF this is done
    # Yes, master Yoda
    self.Weight = boloWeight

    self.MessHand.debug("weights computed in "+str(myTiming))

  # ---------------------------------------------------------------------
  # ---------------------------------------------------------------------
  def writeFFCF(self,outFile='ffcf.txt'):
        """
        NAM: writeFFCF (method)
        DES: store current correlated noise flat field to a file
        INP: (string) file: complete name of output file
        """

        # filename = "ffcf_%s.txt"%(scanNum)
        
        try:
            # f = file(BoaConfig.rcpPath+outFile,'w') # how to addrsss rcppath?
            f = file(outFile,'w')
        except IOError:
            self.MessHand.error("could not open file %s in write mode"%(outFile))
            return

        # Write header
        f.write("! FFCF_CN \n")
        # Write parameters for all channels
	for i in range(len(self.BolometerArray.FlagHandler.getFlags())):
            f.write("%i %f %i \n"% \
                    (i+1, self.FFCF_CN[i,i],int(self.BolometerArray.FlagHandler.getFlags()[i])))
        f.close()
  # ---------------------------------------------------------------------
  def readFFCF(self,inFile='ffcf.txt'):
        """
        NAM: readFFCF (method)
        DES: 
        INP: (string) inFile: complete name of file to read in
        """

        try:
            f = file(inFile)
        except IOError:
            self.MessHand.error("could not open file %s"%(inFile))
            return

        # read and process file
        param = f.readlines()
        f.close()
        ff, flag, chan = [], [], []   # local lists to store FFCF and flag
        
        for i in range(len(param)-1):	        # -1: skip last line
            if param[i][0] != '!':              # skip comments
                tmp = string.split(param[i])
                chan.append(string.atof(tmp[0]))
                ff.append(string.atof(tmp[1]))
                flag.append(string.atof(tmp[2]))

        for i in range(len(ff)):
            ic = int(chan[i])
            self.FFCF_CN[ic,ic] = ff[i]
            if (flag[i] != 0.):
                self.flagChannels(ic)
  #-----------------------------------------------------------------------
  def __correlatedNoiseFFCF(self, chanList=[], skynoise=0, minSlope=0.1, maxSlope=10.0, plot=0, chanRef=-1):
    """
    DES: compute correlation factor relative to given reference channel or skynois
    INP: (l i) chanList : list of channel to correlate (default current list)
         (l)   skynoise : correlate skynoise[channel] not signal[channel] to signal
         (f)   minSlope : limit slope of least squares fit (default 0.1)
         (f)   maxSlope : limit slope of least squares fit (default 10)
         (l)       plot : plot the correlation and the fit (default no)
         (i)    chanRef : reference channel in case of plotting (default : first in chanList, must be in chanList!)
    """

    # This method is very robust and does not depend on the chanRef
    # TODO: it is possible to make it faster by needing a chanRef, i.e
    # just correlate to this one, however, the choice of the chanRef
    # can be problematic
    #
    # FB20070401 fixed fit plotting
    
    if not self._existData(): 
      self.MessHand.error(" (correlate) no data available! ") 
      return

    chanList        = self.BolometerArray.checkChanList(chanList)
    chanListIndexes = self.BolometerArray.getChanIndex(chanList)

    
    if chanRef == -1:
        chanRef       = chanList[0]
        chanRefIndex  = chanListIndexes[0]
    else:
        chanRef = self.BolometerArray.checkChanList(chanRef) [0]
        if not chanRef: 
            self.MessHand.error("not a valid reference channel")
            return
        else:
            #chanRefIndex = list(chanList).index(chanRef)
            chanRefIndex=self.BolometerArray.getChanIndex(chanRef)[0]

    if not chanRef in chanList:
        self.MessHand.error("reference channel must also be in chanList")
        return

    chanRefIndexList = nonzero( equal(chanListIndexes,chanRefIndex) )[0]
    
    dataFlags   = self.FlagHandler.getFlags()
    CorrelateTo = self.Data
                       
    if skynoise:    # if correlation against skynoise
        CorrelateFrom = self.CorrelatedNoise
    else:     # if correlation against a channel
        CorrelateFrom = self.Data

    slopes, intercep, FFCF = fSNF.correlationfit(CorrelateFrom,CorrelateTo,dataFlags,chanListIndexes, \
                                                 array([minSlope,maxSlope]))

    self.FFCF_CN = FFCF   # make public              dimension: arraysize x arraysize
    self.Slopes = slopes  # safe for diagnostic only dimension: chanList  x chanList
    
    if plot:
      self.plotCorrel(chanRef=chanRef,chanList=chanList,skynoise=skynoise)
      dataX = [Plot.xAxis['limits']]*len(chanList)
      dataY = []

      # note that slopes and intercep have the size of chanList in this case
      for i in range(len(chanList)):
          # loop i=0,1,2,.... through nonflg chanList
          dataY.append( array(Plot.xAxis['limits']) * slopes[chanRefIndexList,i] + intercep[chanRefIndexList,i] )
          
      MultiPlot.plot(chanList,dataX,dataY,overplot=1,ci=2,style='l')

  #-----------------------------------------------------------------------
  def __computeCorrelatedNoise(self, chanList=[], clip = 4.0, fastnoise=0):
    """
    DES: compute correlated noise must be run after computeWeight() and correlatedFFCF() with the same chanList
    INP: (i list) chanList : list of channel to use (default: all; [] : current list)
         (f)          clip : the limit where to use the data (+- clip*rms - default 5)
    """
    
    chanList        = self.BolometerArray.checkChanList(chanList)
    chanListIndexes = self.BolometerArray.getChanIndex(chanList)

    data            = self.Data
    dataFlags       = self.FlagHandler.getFlags()
    boloWeight      = self.Weight
    FFCF            = self.FFCF_CN
    nChan           = self.BolometerArray.NUsedChannels
    nInt            = self.ScanParam.NInt    # number of integration points
    
    correlatedNoise = zeros((nInt,nChan),Float32)
    
    correlatedNoise      = fSNF.correlatednoise(data, dataFlags, chanListIndexes, \
                                                boloWeight, FFCF, clip)
    
    # Take only the 7 largest correlations for each channel (columns)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #print shape(correlatedNoise)
    if fastnoise:
        for i in range(size(correlatedNoise,1)):
            corrNoiseColumn = correlatedNoise[:,i]
            corrNoiseSort = sort(corrNoiseColumn)
            #    #corrNoiseSort = corrNoiseSort.reverse()
            
            maskC = where(corrNoiseColumn >= corrNoiseSort[-7], 1, 0).astype(correlatedNoise.typecode())
            correlatedNoise[:,i] = maskC * correlatedNoise[:,i]
    
    self.CorrelatedNoise = correlatedNoise

  #-----------------------------------------------------------------------
  def correlatedNoiseRemoval(self, chanList=[], threshold=1.e-3, iterMax=4, plot=0, \
                             coreRadius=30, beta=2., chanRef=17, fastnoise=0):
    """
    DES: remove the correlated noise from the data. NOTE: THIS METHOD IS EXPERIMENTAL AND MAY NOT WORK PROPERLY
         ON ALL INSTALLATIONS! If you are unsure, use medianNoiseRemoval or corrPCA for the removal of
         correlated noise. 
    INP: (i list) chanList : list of channel to flag (default: all; [] : current list)
         (f)     threshold : threshold value for the Flat Field Correction Factor (in %, default 1.e-3)
         (i)       iterMax : maximum number of iteration
         (i)         plot  : plot or not to plot (def 0)
         (i)      coreRadius: core radius for weight taper beta profile
         (i)      chanRef :  reference channel to start with
    """

    if not self.__corMatrixDone:
        self.computeWeights()

    chanList        = self.BolometerArray.checkChanList(chanList)
    chanListIndexes = self.BolometerArray.getChanIndex(chanList)
    data            = self.Data
    chanRefIndex    = self.BolometerArray.getChanIndex(chanRef)[0]
    chanRefIndexList = nonzero( equal(chanListIndexes,chanRefIndex) )[0]
    #print chanRefIndexList
    
    # First iteration to compute the FFCF
    if plot:
        DeviceHandler.selectDev(1)
    self._DataAna__correlatedNoiseFFCF(chanList=chanList, skynoise=0, plot=plot, chanRef=chanRef)
    if plot:
        if size(DeviceHandler.DevList)==1:
            DeviceHandler.openDev()
        DeviceHandler.selectDev(2)
        BogliConfig.point['size']=3
        Plot.plot(self.Slopes[chanRefIndexList,::],labelY='slope',limitsY=[0,2],ci=1,style='l')
        BogliConfig.point['size']=0.01
    ref = self.FFCF_CN
    self.computeWeights(minCorr=0., a=0.95, b=2.0, core=coreRadius, beta=beta)
    self._DataAna__computeCorrelatedNoise(chanList=chanList,clip=5.,fastnoise=fastnoise)

    # -------------------------------------------- Main loop to estimate the FFCF
    iterNum = 1
    while iterNum <= iterMax:
        self.MessHand.longinfo(" - iteration : %i"%(iterNum))

        if not plot:
            self._DataAna__correlatedNoiseFFCF(chanList=chanList, skynoise=1, plot=0)

        if plot:
            DeviceHandler.selectDev(1)
            self._DataAna__correlatedNoiseFFCF(chanList=chanList, skynoise=1, plot=1)
            #matrix = compress2d(self.FFCF_CN,chanListIndexes)
            #Plot.draw(matrix,wedge=1,limitsZ=[0,2],nan=0,style='idl4', \
            #          caption='Flat Field Correction Factor',labelY = 'good Channel index', \
            #          labelX = 'good Channel index')
            DeviceHandler.selectDev(2)
            BogliConfig.point['size']=3
            # ci = iterNum+1
            # plot slopes for first channel
            Plot.plot(self.Slopes[chanRefIndexList,::],labelY='slope',limitsY=[0,2],overplot=0,style='l',ci=1)
            Plot.plot(diagonal(self.Slopes),overplot=1,style='l',ci=2)
            BogliConfig.point['size']=0.01

        change = fStat.f_mean(abs(ravel((ref-self.FFCF_CN)/ref)))
        self.MessHand.info(" FFCF relative change= %7.3f"%(change))
        if change < threshold:
            self.MessHand.info(" FFCF relative change limit reached: break")
            break
        else:
            ref     = self.FFCF_CN
            iterNum = iterNum+1

        self.computeWeights(minCorr=0., a=0.95, b=2.0, core=coreRadius, beta=2.)
        self._DataAna__computeCorrelatedNoise(chanList=chanList,clip=5.)
    # ------------------------------------------------------------------------------
    if iterNum == iterMax:
        self.MessHand.info("maximum number of iteration reached")
        
    correlatedNoise = self.CorrelatedNoise
    FFCF            = self.FFCF_CN

    for iChan in chanListIndexes:    # diagonal FFCF are *Gains* thus multiply !!
        correlatedNoise[:,iChan] = (correlatedNoise[:,iChan]/FFCF[iChan,iChan]).astype(Float32)

    if plot: # plotting correlated noise over signals
        print "... NEXT showing signal and skynoise ...PRESS <Enter>"; raw_input()
        self.signal(skynoise=0)
        self.signal(skynoise=1,overplot=1,ci=2)

    self.MessHand.longinfo("subtracting CN from Data")

    for iChan in chanListIndexes:
        data[:,iChan] = (data[:,iChan]-correlatedNoise[:,iChan]).astype(Float32)

    if plot: 
        print "... NEXT showing residual signal ...PRESS <Enter>"; raw_input()
        self.signal(skynoise=0)

    self.Data       = data
    self.__resetStatistics()



  # -------------------------------------------------------------------
  # ---- correlated noise reduction by PCA ----------------------------
  # -------------------------------------------------------------------


  def corrPCA(self,chanList=[],order=1,subscan=0,minChanNum=0):
      """
      DES: remove the correlated noise from the data
           by principal component analysis, subscan by suscan
      INP: (i list) chanList : list of channel to flag
           (i) order : number of principal components to remove
           (l) subscan : do the PCA subscan by subscan? default no
           (i) minChanNum : minimum number of valid channels to do PCA (default order+2)
      """
      
      chanList        = self.BolometerArray.checkChanList(chanList)
      chanListIndexes = self.BolometerArray.getChanIndex(chanList)

      if not minChanNum:
          minChanNum=order+2

      if (shape(chanList)[0] <= minChanNum):
          self.MessHand.warning('PCA: not enough valid channels, nothing will be done')
          return

      if subscan:      
          for i in range(shape(self.ScanParam.SubscanNum)[0]):
              self.MessHand.info('PCA: Processing subscan '+str(self.ScanParam.SubscanNum[i]))
              # concatenate data of good channels
              first=self.ScanParam.SubscanIndex[0,i]
              last=self.ScanParam.SubscanIndex[1,i]
              pcadata=transpose(take(self.Data[first:last,::],chanListIndexes,axis=1))
              pcadata,eigenvals,eigenvect=principalComponentAnalysis(array(pcadata),order)
              
              j=0
              for iChan in chanListIndexes:
                  #self.CorrelatedNoise[:,iChan] = self.Data[:,iChan] - \
                  #     pcadata[:,j].astype(self.Data.typecode()).copy()
                  self.Data[first:last,iChan] = pcadata[:,j].astype(self.Data.typecode()).copy()
                  j+=1
      else:
          # concatenate data of good channels
          self.MessHand.info('Doing PCA of order '+str(order))
          pcadata=transpose(take(self.Data,chanListIndexes,axis=1))
          pcadata,eigenvals,eigenvect=principalComponentAnalysis(array(pcadata),order)
          self.pca_eigenvalues=eigenvals
          self.pca_eigenvectors=eigenvect
          
          j=0
          for iChan in chanListIndexes:
              self.CorrelatedNoise[:,iChan] = self.Data[:,iChan] - \
                   pcadata[:,j].astype(self.Data.typecode()).copy()
              self.Data[:,iChan] = pcadata[:,j].astype(self.Data.typecode()).copy()
              j+=1

      self.__resetStatistics()


  def corrPCA_old(self,chanList=[],order=1):
      """
      DES: remove the correlated noise from the data
           by principal component analysis
      INP: (i list) chanList : list of channel to flag
           (i) order : number of principal components to remove
                       ---negative value means choose the optimal number---not yet!!!!
      """
      
      chanList        = self.BolometerArray.checkChanList(chanList)
      chanListIndexes = self.BolometerArray.getChanIndex(chanList)

      # concatenate data of good channels
      pcadata=transpose(take(self.Data,chanListIndexes,axis=1))
      
      if not self.__pcaDone:
          # do the PCA
          pcadata,eigenvals,eigenvect=principalComponentAnalysis(array(pcadata),order)
          self.pca_eigenvalues=eigenvals
          self.pca_eigenvectors=eigenvect
          self.pca_components_removed=order
          self.__pcaDone=1

      else:
          m,n=shape(pcadata)
          #m - number of channels
          #n - number of time samples

          componentsRemoved=self.pca_components_removed

          if (componentsRemoved >= order):
              self.MessHand.warning('These components have already been removed;')
              self.MessHand.warning('  nothing will be done')
              return

          eigenvals=self.pca_eigenvalues
          eigenvect=self.pca_eigenvectors

          rawDataAdjust, originalMean = adjustDataPCA(pcadata)          

          eig_index=argsort(eigenvals)
          featureVector=transpose(take(eigenvect,eig_index[(m-order):(m-componentsRemoved)]))
          undesiredFeature=fUtilities.matrixmultiply(transpose(featureVector),rawDataAdjust)

          pcadata=transpose(pcadata-(matrixmultiply(featureVector,\
                undesiredFeature)+originalMean))

          self.pca_components_removed=order
        
      j=0
      for iChan in chanListIndexes:
          self.CorrelatedNoise[:,iChan] = self.Data[:,iChan] - \
                     pcadata[:,j].astype(self.Data.typecode()).copy()
          self.Data[:,iChan] = pcadata[:,j].astype(self.Data.typecode()).copy()
          j+=1

      self.__resetStatistics()
    
  # -------------------------------------------------------------------
  # ---- baseline methods ---------------------------------------------
  # -------------------------------------------------------------------
  def polynomialBaseline(self,chanList=[],order=0, subscan=1, plot=0, subtract=1, secant=0, returnCoeff=0):
      """
      DES: polynomial baseline removal on the Data.
      INP: (i list) channel     : list of channel to flag (default: all; [] : current list)
               (i) order        : polynomial order,>0 
               (l) subscan      : compute baseline per subscan (default: yes)
               (l) plot         : plot the signal and the fitted polynomials (default: no)
               (l) subtract     : subtract the polynomial from the data (default: yes)
               (l) secant       : also fit a csc(elevation) normalization?
                                  (useful for removing scan-synchronous signals in circular scan patterns)
                                  (default: no)
               (l) returnCoeff  : return polynomial coefficients (default: no)
 """

      self.MessHand.debug('basePoly start...')
      
      # check polynomial order
      if order<0:
          self.MessHand.error("polynomial order must be positive! ") 
          return
      if order>15:
          self.MessHand.warning("fitting is known to fail with high-order polynomial")

      # check channel list
      chanList = self.BolometerArray.checkChanList(chanList)
      if len(chanList)<1: 
          self.MessHand.error("no valid channel")
          return

      chanListIndexes = self.BolometerArray.getChanIndex(chanList)
  
      Data            = self.Data
      DataFlags       = self.FlagHandler.getFlags().copy()
      timeRaw         = self.ScanParam.get('MJD',flag='None')
      timeSubRaw      = self.ScanParam.get('MJD',flag='None')

      mean,stddev=fStat.f_mean(timeRaw),fStat.f_rms(timeRaw,fStat.f_mean(timeRaw))
      time=(timeRaw-mean)/stddev
      timeSub=(timeSubRaw-mean)/stddev
      #time=timeRaw
      #timeSub=timeSubRaw

      if secant:
          elev=self.ScanParam.get('El',flag='None')
          #sec_el=1./(1.-sin(elev*(pi/180.))*sin(elev*(pi/180.)))
          # note that in spite of the name, sec_el is the COSECANT of the elevation
          sec_el=1./(sin(elev*(pi/180.)))
          
      # retrieve subscan index
      if subscan:
          subscanIndex = self.ScanParam.SubscanIndex
      else:
          subscanIndex = array([[0],[self.ScanParam.NInt]])

      if subscan and len(subscanIndex[0]) > 0:
          for i in range(1,len(subscanIndex[0])):
              timeSub[subscanIndex[0,i]:subscanIndex[1,i]]-=time[subscanIndex[0,i]]

              
      if not secant:
          poly = fBaseline.arrayfitpoly_s(chanListIndexes, Data, DataFlags,
                                          timeSub, subscanIndex, order)
      else:
          poly = fBaseline.arrayfitpolysecant_s(chanListIndexes, Data, DataFlags,
                                                timeSub, subscanIndex, sec_el, order)


      # TODO FIX :
      # In the case of scan or single subscan fortran receive a 1d
      # array for iPoly instead of a [1,x] array.
      
      if plot:
          self.signal(chanList)
          dataX = []
          dataY = []
          
          for i in range(len(chanList)):

              iPoly = poly[i,:,::]
   
              iData = Data[:,chanListIndexes[i]]
              X, n = fUtilities.compress(time,DataFlags[:,chanListIndexes[i]],0)
              X = (X * stddev) + mean
              dataX.append(X[:n])
              if secant:
                  Y = fBaseline.evalchunkedpolysecant(timeSub, sec_el, subscanIndex, iPoly)
              else:
                  Y = fBaseline.evalchunkedpoly(timeSub, subscanIndex, iPoly)

              Y, n = fUtilities.compress(Y,DataFlags[:,chanListIndexes[i]],0)
              dataY.append(Y[:n])
              
          MultiPlot.plot(chanList,dataX,dataY,overplot=1,ci=2,style='l')

      if subtract:
          if len(chanListIndexes) > 0:
              if secant:
                  Data = fBaseline.subtractpolysecant(chanListIndexes, Data,
                                                      timeSub, sec_el, subscanIndex, poly)
              else:
                  Data = fBaseline.subtractpoly(chanListIndexes, Data,
                                                timeSub, subscanIndex, poly)
              self.Data = Data
              self.__resetStatistics()
          else:
              self.MessHand.warning("Due to some bug (f2py?) subtraction is not possible on single channel")

      if returnCoeff:
          return(poly)

      self.MessHand.debug('... basePoly end')

  # -----------------------------------------------------------------
  def medianBaseline(self,chanList=[],subscan=1,order=0):
      """
      DES: baseline: Remove median value per channel and per subscan
      INP: (i list) channel : list of channels to process (default: [] = current list)
                (l) subscan : compute baseline per subscan (default: yes)
                (i)   order : polynomial order (default: 0)
      """

      self.MessHand.debug('medianBaseline start...')

      # check channel list
      chanList = self.BolometerArray.checkChanList(chanList)
      if len(chanList)<1: 
          self.MessHand.error("no valid channel")
          return
  
      if order > 1:
          self.MessHand.warning("order > 1 not implemented yet - subtracting order 1")
          order = 1
          
      if not self.__statisticsDone:
          self.__statistics()

      if subscan:
          nbSub    = len(self.ScanParam.SubscanNum)
          for chan in chanList:
              chanNum = self.BolometerArray.getChanIndex(chan)[0]
              for sub in range(nbSub):
                  lo = self.ScanParam.SubscanIndex[0,sub]
                  hi = self.ScanParam.SubscanIndex[1,sub]
                  if order:
		      # compute median of delta_signal / delta_time
                      t = self.getChanData('mjd',chan,
                                           subscans=[self.ScanParam.SubscanNum[sub]])
                      tall = self.getChanData('mjd',chan,
                                              subscans=[self.ScanParam.SubscanNum[sub]],
                                              flag='None')
                      f = self.getChanData('flux',chan,
                                           subscans=[self.ScanParam.SubscanNum[sub]])

		      dt = t-t[0]
		      df = f-f[0]
		      med = fStat.f_median(df[1::]/dt[1::])
		      # compute product median * time
		      med_t = array(med,'d')*tall
                      self.Data[lo:hi,chanNum] -= med_t.astype('f')
	          else:
		      # order 0
		      self.Data[lo:hi,chanNum] -= self.ChanMed_s[chanNum,sub]
      else:
          for chan in chanList:
              chanNum = self.BolometerArray.getChanIndex(chan)[0]
              self.Data[:,chanNum] -= array(self.ChanMed[chanNum]).astype(Float32)

      self.__resetStatistics()
      if order:
          # a first-order has been subtracted, need to subtract 0-order afterwards
	  self.medianBaseline(order=0,subscan=subscan,chanList=chanList)
	  
      self.MessHand.debug('... medianBaseline end')


  # -----------------------------------------------------------------
  def zeroStart(self,chanList=[],subscan = 0):
      """
      DES: make signal start at zero
      INP: (i list) channel : list of channels to process (default: [] = current list)
               (l) subscan  : compute zero per subscan? (default: no)
      """

      self.MessHand.debug('zeroStart start...')

      # check channel list
      chanList = self.BolometerArray.checkChanList(chanList)
      if len(chanList)<1: 
          self.MessHand.error("no valid channel")
          return
      
      if subscan:
          nbSub    = len(self.ScanParam.SubscanNum)
          for chan in chanList:
              chanNum = self.BolometerArray.getChanIndex(chan)[0]
              for sub in range(nbSub):
                  # get only non-flagged data
                  subnum = self.ScanParam.SubscanNum[sub]
                  flux = self.getChanData('flux',chan,subscans=[subnum])
                  lo = self.ScanParam.SubscanIndex[0,sub]
                  hi = self.ScanParam.SubscanIndex[1,sub]
                  self.Data[lo:hi,chanNum] = self.Data[lo:hi,chanNum] - \
                                             array(flux[0],'f')
      else:
          for chan in chanList:
              chanNum = self.BolometerArray.getChanIndex(chan)[0]
              flux = self.getChanData('flux',chan)
              self.Data[:,chanNum] = self.Data[:,chanNum] - array(flux[0],'f')

      self.__resetStatistics()
      self.MessHand.debug('... zeroStart end')

  # -----------------------------------------------------------------
  def zeroEnds(self,chanList=[],subscan=0):
      """
      DES: make signal start AND end at zero, by subtracting an order-1 baseline
      INP: (i list) channel : list of channels to process (default: [] = current list)
               (l) subscan  : compute baseline per subscan? (default: no)
      """

      self.MessHand.debug('zeroEnds start...')

      # check channel list
      chanList = self.BolometerArray.checkChanList(chanList)
      if len(chanList)<1: 
          self.MessHand.error("no valid channel")
          return
      
      # get full time array (also on flagged timestamps)
      mjd = self.ScanParam.get('mjd',flag='None')

      if subscan:
          nbSub    = len(self.ScanParam.SubscanNum)
          for chan in chanList:
              chanNum = self.BolometerArray.getChanIndex(chan)[0]
              for sub in range(nbSub):
                  lo = self.ScanParam.SubscanIndex[0,sub]
                  hi = self.ScanParam.SubscanIndex[1,sub]
                  tt = self.getChanData('mjd',chan,
                                        subscans=[self.ScanParam.SubscanNum[sub]])
                  ss = self.getChanData('flux',chan,
                                        subscans=[self.ScanParam.SubscanNum[sub]])
                  if len(tt) > 0.:
                      slope = (ss[-1]-ss[0]) / (tt[-1]-tt[0])
                      zero  = ss[0] - slope*tt[0]
                      base1 = slope * mjd[lo:hi] + zero
                      base1 = base1.astype('f')
                      self.Data[lo:hi,chanNum] = self.Data[lo:hi,chanNum] - base1
      else:
          for chan in chanList:
              chanNum = self.BolometerArray.getChanIndex(chan)[0]
              tt = self.getChanData('mjd',chan)
              ss = self.getChanData('flux',chan)
              slope = (ss[-1]-ss[0]) / (tt[-1]-tt[0])
              zero  = ss[0] - slope*tt[0]
              base1 = slope * mjd + zero
              self.Data[:,chanNum] = self.Data[:,chanNum] - base1.astype('f')

      self.__resetStatistics()
      self.MessHand.debug('... zeroEnds end')

  # -----------------------------------------------------------------
  def medianFilter(self,chanList=[],window=20,subtract=1,plot=0,limitsX=[], limitsY=[]):
    """
    DES: median filtering: remove median values computed over sliding window
    INP: (i list)       chanList : list of channels to process (default: [] = current list)
    OPT: (i)              window : number of samples to compute median
         (l)            subtract : subtract from data? (default: yes)
         (l)                plot : plot the result? (default: no)
         (2elts array) limitsX/Y : limits to use in X/Y for the plot 
    """

    self.MessHand.debug('medianFilter start...')

    # check channel list
    chanList = self.BolometerArray.checkChanList(chanList)
    if len(chanList)<1: 
      self.MessHand.error("no valid channel")
      return

    if plot:
	self.signal(chanList,limitsX=limitsX, limitsY=limitsY)

    median = []
    for chan in chanList:
	chanNum = self.BolometerArray.getChanIndex(chan)[0]
        data = self.getChanData('flux',chan)
        nbSamp = len(data)
        median1Chan = []
        for i in range(nbSamp):
            if i < window/2:
                median1Chan.append(fStat.f_median(data[:i+1]))
            elif i > nbSamp-window/2:
                median1Chan.append(fStat.f_median(data[i::]))
            else:
                median1Chan.append(fStat.f_median(data[i-window/2:i+window/2]))

        if subtract:
            #TODO: take care of flagged data!
	    self.Data[:,chanNum] = self.Data[:,chanNum] - array(median1Chan).astype(Float32)
	
	median.append(array(median1Chan).astype(Float32))
    
    if plot:
	dataX = self.getChanListData('MJD',chanList)        
	MultiPlot.plot(chanList,dataX,median,style='l',ci=2,overplot=1, \
		labelX="MJD - MJD(0) [sec]",labelY="Flux density [arb.u.]")
	
    if subtract:
        self.__resetStatistics()
    self.MessHand.debug('... medianFilter end')
    
    return median
      
  # -----------------------------------------------------------------
  def __computeMeanSignal(self,chanList=[]):
      """
      DES: compute mean value of the signals at each timestamp
      INP: (i list) channel: list of channels to process (default: [] = current list)
      OUT: this function returns a 1D array containing the mean value of all non-flagged
           data at each timestamp
      """

      self.MessHand.debug('computeMeanSignal start...')
      chanList = self.BolometerArray.checkChanList(chanList)
      if len(chanList)<1:
          self.MessHand.error("no valid channel")
          return

      # get all data, also flagged ones - flags will be handled in fortran
      fluxes = self.getChanListData('flux',chanList,channelFlag='None',dataFlag='None')
      flags  = self.getChanListData('flag',chanList,channelFlag='None',dataFlag='None')

      # Here the 1st dimension correspond to channels - that's what we want
      Mean = fStat.arraymean(fluxes,flags)

      # Now keep only values where timestamps are not flagged
      good = nonzero(self.ScanParam.FlagHandler.isUnsetMask())
      result = take(Mean,good)
      self.MessHand.debug('...computeMeanSignal end')
      return result
    
  # -----------------------------------------------------------------
  def __computeMedianSignal(self,chanList=[]):
      """
      DES: compute median value of the signals at each timestamp
      INP: (i list) channel: list of channels to process (default: [] = current list)
      OUT: this function returns a 1D array containing the median value of all non-flagged
           data at each timestamp
      """

      self.MessHand.debug('computeMedianSignal start...')
      chanList = self.BolometerArray.checkChanList(chanList)
      if len(chanList)<1:
          self.MessHand.error("no valid channel")
          return

      # get all data, also flagged ones - flags will be handled in fortran
      fluxes = self.getChanListData('flux',chanList,channelFlag='None',dataFlag='None')
      flags  = self.getChanListData('flag',chanList,channelFlag='None',dataFlag='None')

      # Here the 1st dimension correspond to channels - that's what we want
      Med = fStat.arraymedian(fluxes,flags)

      # Now keep only values where timestamps are not flagged
      good = nonzero(self.ScanParam.FlagHandler.isUnsetMask())
      result = take(Med,good)
      self.MessHand.debug('...computeMedianSignal end')
      return result
    
  # -----------------------------------------------------------------
  def __computeMedianAbsSignal(self,chanList=[]):
      """
      DES: compute median value of absolute values of signals at each timestamp
      INP: (i list) channel: list of channels to process (default: [] = current list)
      OUT: this function returns a 1D array containing the median value of all non-flagged
           data at each timestamp
      """

      self.MessHand.debug('computeMedianAbsSignal start...')
      chanList = self.BolometerArray.checkChanList(chanList)
      if len(chanList)<1:
          self.MessHand.error("no valid channel")
          return

      # get all data, also flagged ones - flags will be handled in fortran
      fluxes = self.getChanListData('flux',chanList,channelFlag='None',dataFlag='None')
      fluxes = absolute(fluxes)
      flags  = self.getChanListData('flag',chanList,channelFlag='None',dataFlag='None')

      # Here the 1st dimension correspond to channels - that's what we want
      Med = fStat.arraymedian(fluxes,flags)

      # Now keep only values where timestamps are not flagged
      good = nonzero(self.ScanParam.FlagHandler.isUnsetMask())
      result = take(Med,good)
      self.MessHand.debug('...computeMedianAbsSignal end')
      return result
    
  # -----------------------------------------------------------------
  def __computeMedianFlatField(self,chanList=[],chanRef=0):
    """
    DES: compute flat field as median of signal ratio between channels
    INP: (i list) channel: list of channels to process (default: [] = current list)
         (i)      chanRef: reference channel number (def: 0, i.e. use Reference chan
                           can also be -1, i.e. use mean signal as reference
                                    or -2, i.e. use median signal as reference
    """

    self.MessHand.debug('computeMedianFlatField start...')

    # check channel list
    chanList = self.BolometerArray.checkChanList(chanList)
    if len(chanList)<1: 
        self.MessHand.error("no valid channel")
        return

    if chanRef == 0:
        chanRef = self.BolometerArray.RefChannel
    if chanRef > -1:
        chanRefIndex = self.BolometerArray.getChanIndex(chanRef)[0]
        if chanRefIndex ==  -1:
            self.MessHand.error("chanRef: channel not used") 
            return
        refFlags = self.FlagHandler.getFlags()[:,chanRefIndex].astype(Int32)
    elif chanRef == -1:
        refFlags = self.ScanParam.FlagHandler.getFlags()
        refSignalFull = self._DataAna__computeMeanSignal(chanList)
                        # full non-flagged time stream
        goodFlagRef   = nonzero(self.ScanParam.FlagHandler.isUnsetMask())
                        # corresponding indices
    elif chanRef == -2:
        refFlags = self.ScanParam.FlagHandler.getFlags()
        refSignalFull = self._DataAna__computeMedianSignal(chanList)
        goodFlagRef   = nonzero(self.ScanParam.FlagHandler.isUnsetMask())

    for chan in chanList:
        num = self.BolometerArray.getChanIndex(chan)[0]
        if chan == chanRef:
            self.FF_Median[num] = 1.
        else:
            # we consider only datapoints where both ref chan. and current
            # chan are not flagged
            chanFlags  = self.FlagHandler.getFlags()[:,num].astype(Int32)
            chanSignal = self.getChanData('flux',chan,flag2 = refFlags)
            if chanRef > -1:
                refSignal  = self.getChanData('flux',chanRef,flag2 = chanFlags)
            else:
                chanFlags = take(chanFlags,goodFlagRef)
                flagHandler = BoaFlagHandler.createFlagHandler(chanFlags)
                good = nonzero(flagHandler.isUnsetMask())
                refSignal = take(refSignalFull,good)

            if len(refSignal) == 0 or len(chanSignal) == 0:
                # this happens if too much is flagged
                self.MessHand.warning("Too much flags, could not determine FF for chan. %s"%(chan))
                ratio = 1.
            else:
                ratio = fStat.f_median(chanSignal / refSignal)
            self.FF_Median[num] = ratio

    self.MessHand.debug('... computeMedianFlatField end')

  # -----------------------------------------------------------------
  def __applyFlatField(self,chanList=[]):
    """
    DES: divide signals by bolo gains to normalise them
    INP: (i list) channel: list of channels to process (default: [] = current list)
    """

    self.MessHand.debug('applyFlatField start...')

    # check channel list
    chanList = self.BolometerArray.checkChanList(chanList)
    if len(chanList)<1: 
        self.MessHand.error("no valid channel")
        return

    for chan in chanList:
        num = self.BolometerArray.getChanIndex(chan)[0]
        # self.Data[:,num] = self.Data[:,num] / array((self.FFCF_Gain[num]),Float32)
        # corrected 2017-07-04 (FSc)
        self.Data[:,num] = self.Data[:,num] / array((self.FFCF_Gain[chan-1]),Float32)
        
    self.__resetStatistics()
    self.MessHand.debug('... applyFlatField end')

  def flatfield(self,chanList=[],method='point'):
      """
      DES: divide signals by bolo gains to normalise them
      INP: (i list) channel: list of channels to process (default: [] = current list)
           (str)    method : choose which flat field to apply:
                             - point [default] = use point source relative gains
                             - median = use correlate noise relative gains
                             - extend = use relative gains to extended emission
      """
      self.MessHand.debug('flatfield start...')
      
      # check channel list
      chanList = self.BolometerArray.checkChanList(chanList)
      if len(chanList)<1: 
          self.MessHand.error("no valid channel")
          return

      if method == 'point':
        self.FFCF_Gain = self.BolometerArray.Gain
      elif method == 'median':
        self.FFCF_Gain = self.FF_Median
      elif method == 'extend':
        self.FFCF_Gain = self.BolometerArray.ExtGain
      self.__applyFlatField(chanList=chanList)
      
  # -----------------------------------------------------------------
  def __computeMedianNoise(self,chanList=[],mean=0):
    """
    DES: compute median noise, i.e. median of all bolos (normalised!) at
         each individual timestamp
    INP: (i list) channel: list of channels to process (default: [] = current list)      
         (log)  mean : compute mean noise instead of median noise
    """

    self.MessHand.debug('computeMedianNoise start...')

    # check channel list
    chanList = self.BolometerArray.checkChanList(chanList)
    if len(chanList)<1: 
        self.MessHand.error("no valid channel")
        return
    alldata = self.getChanListData('flux',chanList,channelFlag='None',dataFlag='None')
    alldata = array(alldata,Float32)
    allFlags = array(self.getChanListData('flag',chanList,channelFlag='None',dataFlag='None'))

    # Correct for median flat field
    for i in range(len(chanList)):
        num = self.BolometerArray.getChanIndex(chanList[i])[0]
        alldata[i,::] = alldata[i,::] / array((self.FF_Median[num]),Float32)
    
    # Number of integrations
    nInt = self.ScanParam.NInt
    #skynoise = fStat.arraymedian(alldata,allFlags)
    if mean:
        skynoise = fStat.arraymean(alldata,allFlags)
    else:
        skynoise = fStat.arraymedian(alldata,allFlags)
    skynoise = repeat(reshape(skynoise,(1,nInt)),self.BolometerArray.NUsedChannels)
    self.Skynoise = fUtilities.as_column_major_storage(transpose(skynoise))
    
    self.MessHand.debug('... computeMedianNoise end')
        
  #-----------------------------------------------------------------------

  def medianNoiseLocal(self, chanList=[], chanRef=-2, computeFF=1,factor=1.,
                          numCorr=7, minDist=0., selByDist=0,outputChanList=0,mean=0):

      """
      DES: remove median noise from the data by using only the n most correlated channels w.r.t. each bolometer
      INP: (i list)    chanList : list of channels (default: [] = current list)
           (int)        chanRef : -1 = compute relative gains w.r.t. mean signal
                                  -2 = compute relative gains w.r.t. median signal (default)
           (log)      computeFF : compute skynoise FF (def.) or use existing FF_Median?
           (float)       factor : fraction of skynoise to be subtracted (default: 100%)
           (int)        numCorr : number of (most correlated) channels to use to compute the sky noise for each channel
           (float)      minDist : minimum distance on sky, in ARCSEC, between channels to be considered (useful for extended emission) 
           (int)      selByDist : set this to select the n closest channels (outside minDist) instead of the most correlated ones
           (int) outputChanList : set this to obtain the list of most correlated channels in output 
           (log)          mean  : use mean instead of median (useful for linear filtering)
      """

      
      # check channel list
      chanList = self.BolometerArray.checkChanList(chanList)
      if len(chanList)<1: 
          self.MessHand.error("no valid channel")
          return
      chanListIndices = self.BolometerArray.getChanIndex(chanList)

      # compute the correlation matrix (stored in self.CorMatrix)
      self.computeWeights()

      # copy the data to have something to work on 
      # (while data array itself is kept intact until the end)
      wdata=copy.deepcopy(self.Data)

      #loop over channles, remove noise from one channel at a time
      chanIndexList=[]
      noiseChanList=[]
      for chanIn in chanList:

          # get chan index
          chanIndex=self.BolometerArray.getChanIndex(chanIn)[0]
          chanIndexList.append(chanIndex)
	  
          # get channel separations
          chanSep=take(self.BolometerArray.ChannelSep[:,chanIndex], chanListIndices)

          #determine correlations
          corr=take(self.CorMatrix[:,chanIndex], chanListIndices)

          if selByDist:
              mask=where((corr > -1.) and (corr <= 1.) and (chanSep > minDist),1,0)
              newList=compress(mask,chanList)
              newSep=compress(mask,chanSep)
              maxCorList=tolist_boa(take(newList,argsort(newSep)[:(abs(numCorr))]))
          else:
              mask=where((corr > -1.) and (corr <= 1.) and (chanSep > minDist),1,0)
              corr=compress(mask,corr)
              newList=compress(mask,chanList)
              # select
              maxCorList=tolist_boa(take(newList,argsort(corr)[-(abs(numCorr))::]))

          # make sure current channel is in list
          if not chanIn in maxCorList:
              maxCorList.append(chanIn)

          print str(chanIn)+' '+str(maxCorList)
	  noiseChanList.append(maxCorList)
          
	  # compute flatfield and median noise for new list
          if computeFF:
              self.__computeMedianFlatField(maxCorList,chanRef=chanRef)

          self.__computeMedianNoise(maxCorList,mean=mean)
          
          wdata[:,chanIndex] = self.Data[:,chanIndex] - self.Skynoise[:,chanIndex] * \
                               array(self.FF_Median[chanIndex]*factor,Float32)

      self.Data=copy.deepcopy(wdata)
          
      self.__resetStatistics()
      
      if outputChanList:
          noiseLocalOut = []
	  noiseLocalOut.append(chanIndexList)
	  noiseLocalOut.append(noiseChanList)
	  return noiseLocalOut
   
         
  #-----------------------------------------------------------------------

  def medianNoiseFromList(self, cList, chanRef=-2, computeFF=1,factor=1.):

      """
      DES: remove median noise from the data by using only the channels provided in input
      INP: (i list)       cList : list of channels as returned by MedianNoiseLocal
           (int)        chanRef : reference channel number (default: RefChannel;
                                  -1 = compute relative gains w.r.t. mean signal
                                  -2 = compute relative gains w.r.t. median signal
           (log)      computeFF : compute skynoise FF (def.) or use existing FF_Median?
           (float)       factor : fraction of skynoise to be subtracted (default: 100%)
      """

      
      # check input channel list
      # Must consist of two sub-lists
      if (len(cList)!=2): 
          self.MessHand.error("invalid channel list")
          return
      # First list contains the channel indices
      chanListIndices = cList[0]
      nChan = len(chanListIndices)
      # Second list contains the channels from which to estimate the noise
      noiseChanList = cList[1]
      if (len(noiseChanList)!=nChan):
          self.MessHand.error("invalid channel list")
          return

      # copy the data to have something to work on 
      # (while data array itself is kept intact until the end)
      wdata=copy.deepcopy(self.Data)

      #loop over channles, remove noise from one channel at a time
      for ind,maxCorList in enumerate(noiseChanList):
	  
	  # compute flatfield and median noise for new list
          if computeFF:
              self.__computeMedianFlatField(maxCorList,chanRef=chanRef)
          self.__computeMedianNoise(maxCorList)
          
          # Remove median noise
	  wdata[:,chanListIndices[ind]] = self.Data[:,chanListIndices[ind]] - \
	  		self.Skynoise[:,chanListIndices[ind]] * array(self.FF_Median[chanListIndices[ind]]*factor,Float32)

      self.Data=copy.deepcopy(wdata)
          
      self.__resetStatistics()
       
  #-----------------------------------------------------------------------

  
  def medianNoiseRemoval(self, chanList=[], chanRef=0, computeFF=1,
                         factor=1.,nbloop=1,mean=0):
    """
    DES: remove median noise from the data
    INP: (i list) chanList : list of channels (default: [] = current list)
         (int)     chanRef : reference channel number (default: RefChannel;
                             -1 = compute relative gains w.r.t. mean signal
                             -2 = compute relative gains w.r.t. median signal
         (log)   computeFF : compute skynoise FF (def.) or use existing FF_Median?
         (float)    factor : fraction of skynoise to be subtracted (default: 100%)
         (int)      nbloop : number of iterations (default: 1)
         (log)        mean : use mean instead of median (useful for linear filtering)
    """

    self.MessHand.debug('medianNoiseRemoval start...')

    # check channel list
    chanList = self.BolometerArray.checkChanList(chanList)
    if len(chanList)<1: 
        self.MessHand.error("no valid channel")
        return

    for i in range(nbloop):
        if computeFF:
            self.__computeMedianFlatField(chanList,chanRef=chanRef)

        self.__computeMedianNoise(chanList,mean=mean)
        for chan in chanList:
            num = self.BolometerArray.getChanIndex(chan)[0]
            self.Data[:,num] = self.Data[:,num] - self.Skynoise[:,num] * \
                               array(self.FF_Median[num]*factor,Float32)

    self.__resetStatistics()
    self.MessHand.debug('... medianNoiseRemoval end')
    
  #-----------------------------------------------------------------------
  def averageNoiseRemoval(self, chanList=[], chanRef=0):
    """
    DES: remove correlated noise computed as average value of all but
         the reference channel
    INP: (i list) chanList : list of channels (default: [] = all valid channels)
         (int)     chanRef : reference channel number, not used to compute
                             the average noise (default: none
    """
    self.MessHand.debug('averageNoiseRemoval start...')

    # check channel list
    chanList = self.BolometerArray.checkChanList(chanList)
    chanList = chanList.tolist()
    if len(chanList)<1: 
        self.MessHand.error("no valid channel")
        return

    tmpList = copy.copy(chanList)
    if chanRef:
        if not chanList.count(chanRef):
            self.MessHand.warning("Ref. channel not in Chanlist, average")
            self.MessHand.warning("noise not subtracted from ref. channel")
        else:
            tmpList.remove(chanRef)

    if tmpList:
        tmpData = self.getChanListData('flux',chanList=tmpList,dataFlag='None')
        tmpFlag = self.getChanListData('flag',chanList=tmpList,dataFlag='None')
        noise = fStat.arraymean(tmpData,tmpFlag)

        for chan in chanList:
            chanNum = self.BolometerArray.getChanIndex(chan)[0]
            self.Data[:,chanNum] = self.Data[:,chanNum] - noise
    else:
        self.MessHand.warning("No channel available to compute noise")
        
  #-----------------------------------------------------------------------
  #-----------------------------------------------------------------------
  def correctOpacity(self,tau=0.):
      """
      DES: correct for atmospheric opacity
      """
      self.MessHand.debug('correctOpacity start...')
      if not tau:
          self.MessHand.warning("No tau value provided - exiting")
          return

      # use first non-flagged channel to extract elevation
      chan0 = self.BolometerArray.checkChanList ([])[0]
      el = self.getChanData('el',chan0)
      # use median value
      med_el = fStat.f_median(el)
      tau_los = tau/sin(med_el * pi / 180.)
      self.Data *= array(exp(tau_los),'f')

      self.__resetStatistics()
      self.MessHand.debug('... correctOpacity end')

  #--------------------------------------------------------------------------------
  def addSourceModel(self,model,chanList='all',factor=1.):
      """
      DES: add data according to a model map
      INP: (Imgae object) model : the input model map (with WCS)
           (i list)    chanList : the list of channels to work with
           (f)           factor : add model data multiplied with this factor

      """
      # check channel list
      self.MessHand.info('adding source model ...')
      chanList = self.BolometerArray.checkChanList(chanList)
      if len(chanList)<1: 
          self.MessHand.error("no valid channel")
          return
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

      # get the new data + factor x model array
      tmp = fMap.addsource(chanListIndexes, self.Data, model.Data, \
                           XYOffsets, OffsetsAzEl, rotAngles, refChOffsets, \
                           AXIS1, AXIS2, factor)
      # replace self.Data with updated one
      #self.Data = copy.copy(tmp)
      self.Data=tmp
      tmp = 0  # free memory       
      self.MessHand.info('...done')


   # excise / select data --------------------------------------

  def select(self, begin, end):
    """
    DES: retain a portion of the data, defined between two indices
         of the time stream; excise the rest permanetly (no flagging)
    INP: (int) begin : start index
         (int)   end : end index 
    """

    #NInt = self.ScanParam.NInt
    #NCh  = shape(self.Data)[1]

    newData = self.Data[begin:end,:]
    newWeig = self.DataWeights[begin:end,:]
    newFlag = self.FlagHandler.getFlags()[begin:end,:]

    newBasLat = self.ScanParam.BasLat[begin:end]
    newBasLon = self.ScanParam.BasLon[begin:end]
    newAz     = self.ScanParam.Az[begin:end]
    newEl     = self.ScanParam.El[begin:end]
    newAzOff  = self.ScanParam.AzOff[begin:end]
    newElOff  = self.ScanParam.ElOff[begin:end]
    newLatOff = self.ScanParam.LatOff[begin:end]
    newLonOff = self.ScanParam.LonOff[begin:end]
    newRot    = self.ScanParam.Rot[begin:end]
    newPar    = self.ScanParam.ParAngle[begin:end]
    newLST    = self.ScanParam.LST[begin:end]
    newMJD    = self.ScanParam.MJD[begin:end]
    newRA     = self.ScanParam.RA[begin:end]
    newDec    = self.ScanParam.Dec[begin:end]
    newRAOff  = self.ScanParam.RAOff[begin:end]
    newDecOff = self.ScanParam.DecOff[begin:end]
    newFlags  = self.ScanParam.FlagHandler.getFlags()[begin:end]
    newFocX   = self.ScanParam.FocX[begin:end]
    newFocY   = self.ScanParam.FocY[begin:end]
    newFocZ   = self.ScanParam.FocZ[begin:end]

    self.Data        = fUtilities.as_column_major_storage(newData)
    self.DataWeights = fUtilities.as_column_major_storage(newWeig)
    self.FlagHandler = BoaFlagHandler.createFlagHandler(newFlag)

    self.ScanParam.BasLat = fUtilities.as_column_major_storage(newBasLat)
    self.ScanParam.BasLon = fUtilities.as_column_major_storage(newBasLon)
    self.ScanParam.Az     = fUtilities.as_column_major_storage(newAz)
    self.ScanParam.El     = fUtilities.as_column_major_storage(newEl)
    self.ScanParam.AzOff  = fUtilities.as_column_major_storage(newAzOff)
    self.ScanParam.ElOff  = fUtilities.as_column_major_storage(newElOff)
    self.ScanParam.LatOff = fUtilities.as_column_major_storage(newLatOff)
    self.ScanParam.LonOff = fUtilities.as_column_major_storage(newLonOff)
    self.ScanParam.Rot    = fUtilities.as_column_major_storage(newRot)
    self.ScanParam.ParAngle= fUtilities.as_column_major_storage(newPar)
    self.ScanParam.LST    = fUtilities.as_column_major_storage(newLST)
    self.ScanParam.MJD    = fUtilities.as_column_major_storage(newMJD)
    self.ScanParam.RA     = fUtilities.as_column_major_storage(newRA)
    self.ScanParam.Dec    = fUtilities.as_column_major_storage(newDec)
    self.ScanParam.RAOff  = fUtilities.as_column_major_storage(newRAOff)
    self.ScanParam.DecOff = fUtilities.as_column_major_storage(newDecOff)
    self.ScanParam.FocX   = fUtilities.as_column_major_storage(newFocX)
    self.ScanParam.FocY   = fUtilities.as_column_major_storage(newFocY)
    self.ScanParam.FocZ   = fUtilities.as_column_major_storage(newFocZ)
    self.ScanParam.FlagHandler = BoaFlagHandler.createFlagHandler(newFlags)

    self.ScanParam.NInt = end-begin   
    self.ScanParam.findSubscan(combine=10000)
    #for i in range(len(self.ScanParam.SubscanNum)):
    #    self.ScanParam.SubscanIndex[0][i] = self.ScanParam.SubscanIndex[0][i]/2
    #    self.ScanParam.SubscanIndex[1][i] = self.ScanParam.SubscanIndex[1][i]/2

    self._DataAna__resetStatistics()


      
  #-----------------------------------------------------------------------
  # rebin: down sample the data by a factor 2
  #-----------------------------------------------------------------------
  def rebin(self):
    """
    DES: average integrations 2 by 2
    """
    NInt = self.ScanParam.NInt
    NCh  = shape(self.Data)[1]

    newData = zeros((NInt/2,NCh),Float32)
    newWeig = ones((NInt/2,NCh),Float32)
    newFlag = zeros((NInt/2,NCh),Int8)

    newBasLat = zeros((NInt/2),Float32)
    newBasLon = zeros((NInt/2),Float32)
    newAz     = zeros((NInt/2),Float32)
    newEl     = zeros((NInt/2),Float32)
    newAzOff  = zeros((NInt/2),Float32)
    newElOff  = zeros((NInt/2),Float32)
    newLatOff = zeros((NInt/2),Float32)
    newLonOff = zeros((NInt/2),Float32)
    newRot    = zeros((NInt/2),Float32)
    newPar    = zeros((NInt/2),Float32)
    newLST    = zeros((NInt/2),Float32)
    newMJD    = zeros((NInt/2),Float64)
    newRA     = zeros((NInt/2),Float32)
    newDec    = zeros((NInt/2),Float32)
    newRAOff  = zeros((NInt/2),Float32)
    newDecOff = zeros((NInt/2),Float32)
    newFlags  = zeros((NInt/2),Int32)
    newFocX   = zeros((NInt/2),Float32)
    newFocY   = zeros((NInt/2),Float32)
    newFocZ   = zeros((NInt/2),Float32)

    two = array((2.),Float32)
    
    for i in range(NInt/2):
        newData[i,::] = (self.Data[2*i,::] + self.Data[2*i+1,::])/two
        newWeig[i,::] = (self.DataWeights[2*i,::] + self.DataWeights[2*i+1,::])/two
        # Flags: flag new datapoint if one of the two is flagged
        tmpFlag      = bitwise_or(self.FlagHandler.getFlags()[2*i,::], \
                                  self.FlagHandler.getFlags()[2*i+1,::])
        newFlag[i,::] = tmpFlag#.astype(Int8)
        
        # ScanParam attributes
        newFlags[i]  = bitwise_or(self.ScanParam.FlagHandler.getFlags()[2*i], \
                                  self.ScanParam.FlagHandler.getFlags()[2*i+1])
        newBasLat[i] = (self.ScanParam.BasLat[2*i] + self.ScanParam.BasLat[2*i+1])/2.
        newBasLon[i] = (self.ScanParam.BasLon[2*i] + self.ScanParam.BasLon[2*i+1])/2.
        newAz[i]     = (self.ScanParam.Az[2*i]     + self.ScanParam.Az[2*i+1])/2.
        newEl[i]     = (self.ScanParam.El[2*i]     + self.ScanParam.El[2*i+1])/2.
        newAzOff[i]  = (self.ScanParam.AzOff[2*i]  + self.ScanParam.AzOff[2*i+1])/2.
        newElOff[i]  = (self.ScanParam.ElOff[2*i]  + self.ScanParam.ElOff[2*i+1])/2.
        newLatOff[i] = (self.ScanParam.LatOff[2*i] + self.ScanParam.LatOff[2*i+1])/2.
        newLonOff[i] = (self.ScanParam.LonOff[2*i] + self.ScanParam.LonOff[2*i+1])/2.
        newRot[i]    = (self.ScanParam.Rot[2*i]    + self.ScanParam.Rot[2*i+1])/2.
        newPar[i]    = (self.ScanParam.ParAngle[2*i]+ self.ScanParam.ParAngle[2*i+1])/2.
        newLST[i]    = (self.ScanParam.LST[2*i]    + self.ScanParam.LST[2*i+1])/2.
        newMJD[i]    = (self.ScanParam.MJD[2*i]    + self.ScanParam.MJD[2*i+1])/2.
        newRA[i]     = (self.ScanParam.RA[2*i]     + self.ScanParam.RA[2*i+1])/2.
        newDec[i]    = (self.ScanParam.Dec[2*i]    + self.ScanParam.Dec[2*i+1])/2.
        newRAOff[i]  = (self.ScanParam.RAOff[2*i]  + self.ScanParam.RAOff[2*i+1])/2.
        newDecOff[i] = (self.ScanParam.DecOff[2*i] + self.ScanParam.DecOff[2*i+1])/2.
        newFocX[i]   = (self.ScanParam.FocX[2*i]   + self.ScanParam.FocX[2*i+1])/2.
        newFocY[i]   = (self.ScanParam.FocY[2*i]   + self.ScanParam.FocY[2*i+1])/2.
        newFocZ[i]   = (self.ScanParam.FocZ[2*i]   + self.ScanParam.FocZ[2*i+1])/2.

    self.Data        = fUtilities.as_column_major_storage(newData)
    self.DataWeights = fUtilities.as_column_major_storage(newWeig)
    self.FlagHandler = BoaFlagHandler.createFlagHandler(newFlag)

    self.ScanParam.BasLat = fUtilities.as_column_major_storage(newBasLat)
    self.ScanParam.BasLon = fUtilities.as_column_major_storage(newBasLon)
    self.ScanParam.Az     = fUtilities.as_column_major_storage(newAz)
    self.ScanParam.El     = fUtilities.as_column_major_storage(newEl)
    self.ScanParam.AzOff  = fUtilities.as_column_major_storage(newAzOff)
    self.ScanParam.ElOff  = fUtilities.as_column_major_storage(newElOff)
    self.ScanParam.LatOff = fUtilities.as_column_major_storage(newLatOff)
    self.ScanParam.LonOff = fUtilities.as_column_major_storage(newLonOff)
    self.ScanParam.Rot    = fUtilities.as_column_major_storage(newRot)
    self.ScanParam.ParAngle= fUtilities.as_column_major_storage(newPar)
    self.ScanParam.LST    = fUtilities.as_column_major_storage(newLST)
    self.ScanParam.MJD    = fUtilities.as_column_major_storage(newMJD)
    self.ScanParam.RA     = fUtilities.as_column_major_storage(newRA)
    self.ScanParam.Dec    = fUtilities.as_column_major_storage(newDec)
    self.ScanParam.RAOff  = fUtilities.as_column_major_storage(newRAOff)
    self.ScanParam.DecOff = fUtilities.as_column_major_storage(newDecOff)
    self.ScanParam.FocX   = fUtilities.as_column_major_storage(newFocX)
    self.ScanParam.FocY   = fUtilities.as_column_major_storage(newFocY)
    self.ScanParam.FocZ   = fUtilities.as_column_major_storage(newFocZ)
    self.ScanParam.FlagHandler = BoaFlagHandler.createFlagHandler(newFlags)


    self.ScanParam.NInt = NInt/2    
    for i in range(len(self.ScanParam.SubscanNum)):
        self.ScanParam.SubscanIndex[0][i] = self.ScanParam.SubscanIndex[0][i]/2
        self.ScanParam.SubscanIndex[1][i] = self.ScanParam.SubscanIndex[1][i]/2

    self._DataAna__resetStatistics()


  # -------------------------------------------------------------------
  # computeWeight: fill DataWeights attribute
  # -------------------------------------------------------------------

  
      
  
  def computeWeight(self,method='rms',subscan=0, lolim=0.1, hilim=10.0):
    """
    DES: compute weights and store them in DataWeights attribute
    INP: (str)   method : type of weighting (default='rms')
                          'rms' : use 1/rms^2
                          'pow' : use 1/pow^2, where pow is the mean power
                                  between frequencies lolim and hilim
                                  (in Hz)
                 lolim  : low frequency limit for 'pow' method
                          defualt=0.1 Hz
                 hilim  : high freq. limit for 'pow' method
                          defualt=10.0 Hz
                          hilim and lolim are ignored unless method='pow'
         (bool) subscan : compute weight by subscan? default: no
                          ignored if method='pow'
    """
    if not self.__statisticsDone:
      self.__statistics()

    if not method in ['rms','pow']:
        self.MessHand.error("Unknown weighting method - no weight computed")

    if method == 'rms':
        if subscan:
            self.__statistics()
            subnum=self.ScanParam.SubscanNum
            subrms=self.ChanRms_s
            subin=self.ScanParam.SubscanIndex
            for s in range(len(subnum)):
                r = subrms[:,s]
                weight = zeros(shape(r),'f')
                for ch in range(len(r)):
                    if (r[ch] > 0.0):
                        weight[ch] = 1./r[ch]**2  #

                for i in range(subin[0,s],subin[1,s]):
                    self.DataWeights[i,::] = weight.astype(Float32)
                
                    
        else:
            r = self.ChanRms
            weight = zeros(shape(r),'f')
            for i in range(len(r)):
                if str(r[i]) != str(float('nan')):
                    weight[i] = 1./r[i]**2  # NaNs will have zero weight

            for i in range(0,self.ScanParam.NInt):
                self.DataWeights[i,::] = weight.astype(Float32)

    if method == 'pow':
        (c,x,y)=self.plotFFT('all',returnSpectrum=1,plot=0)
        chanListAll=self.BolometerArray.UsedChannels
        weights = zeros(shape(chanListAll),'f')
        ind=0
        for i in range(len(weights)):
            if chanListAll[i] in c:
                xc=x[ind][::]
                yc=y[ind][::]
                mask=where( bitwise_and((xc > lolim), \
                                        (xc < hilim)), 1,0 )
                pows=fStat.f_mean(compress(mask,yc))
                weights[i]=1./(pows**2)
                ind=ind+1
                
        for i in range(0,self.ScanParam.NInt):
            self.DataWeights[i,::] = weights.astype(Float32)
    
    return
    
    
  def slidingWeight(self,chanList=[],nbInteg=50):
      """
      DES: compute weights using 1/rms^2, where rms is computed in
           sliding windows of size nbInteg
      INP: (i list) chanList = the list of channels to compute
           (i)       nbInteg = size of windows (default: 20)
      """
      chanList = self.BolometerArray.checkChanList(chanList)
      if len(chanList)<1: 
          self.MessHand.error("no valid channel")
          return

      for c in chanList:
          index = self.BolometerArray.getChanIndex(c)[0]
          rms = fStat.slidingrms(self.Data[:,index],
                                 self.FlagHandler._aFlags[:,index],
                                 nbInteg)
          wei = 1./rms**2
          self.DataWeights[:,index] = wei.astype(self.DataWeights.typecode())

  # -------------------------------------------------------------------
  # ---- plotting methods ---------------------------------------------
  # -------------------------------------------------------------------

  def plotMean(self,chanList=[], \
               channelFlag=[], plotFlaggedChannels=0, \
               dataFlag=[], plotFlaggedData=0, \
               limitsX=[],limitsY=[], \
               style='l', ci=1, overplot=0, map=0):
    """
    DES: plot mean flux value vs. subscan number
    TODO: flag handling not implemented yet
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
         (logical)       map = plot as a 2D map?
    """

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
    
    dataX = self.getChanListData('subscan',chanList, \
                                 channelFlag=channelFlag, getFlaggedChannels=plotFlaggedChannels, \
                                 dataFlag=dataFlag, getFlaggedData=plotFlaggedData)
    dataY = self.getChanListData('mean_s',chanList, \
                                 channelFlag=channelFlag, getFlaggedChannels=plotFlaggedChannels, \
                                 dataFlag=dataFlag, getFlaggedData=plotFlaggedData)

    xLabel = "Subscan #"
    yLabel = "Mean flux [arb.u.]"
     
    self.MessHand.info("plotting Mean values per subscan")
    if not self.__statisticsDone:
      self.MessHand.warning(" plotting outdated statistics")
     
    if map:
      Plot.draw(transpose(dataY), \
                sizeX=[min(dataX[0]),max(dataX[0])],sizeY=[1,len(chanList)],\
                limitsX=limitsX,\
                labelX=xLabel,labelY='Channel #', \
                caption=self.ScanParam.caption(), wedge=1)
    else:
      MultiPlot.plot(chanList,dataX,dataY,\
                     limitsX = limitsX, limitsY = limitsY, \
                     labelX = xLabel, labelY = yLabel, caption=self.ScanParam.caption(), \
                     style=style,ci=ci,overplot=overplot)
 
  #----------------------------------------------------------------------------
  def plotRms(self,chanList=[], \
              channelFlag=[], plotFlaggedChannels=0, \
              dataFlag=[], plotFlaggedData=0, \
              limitsX=[],limitsY=[], \
              style='l', ci=1, overplot=0, map=0):
    """
    DES: plot flux r.m.s. vs. subscan number
    TODO: flag handling not implemented yet
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
         (logical)       map = plot as a 2D map?
    """

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

    if not self.__statisticsDone:
      self.__statistics()
    
    dataX = self.getChanListData('subscan',chanList, \
                                 channelFlag=channelFlag, getFlaggedChannels=plotFlaggedChannels, \
                                 dataFlag=dataFlag, getFlaggedData=plotFlaggedData)
    dataY = self.getChanListData('rms_s',chanList, \
                                 channelFlag=channelFlag, getFlaggedChannels=plotFlaggedChannels, \
                                 dataFlag=dataFlag, getFlaggedData=plotFlaggedData)
    
    xLabel = "Subscan #"
    yLabel = "Flux r.m.s. [arb.u.]"
    
    self.MessHand.info("plotting r.m.s. per subscan")
    if not map:
            
      MultiPlot.plot(chanList,dataX,dataY,\
                     limitsX = limitsX, limitsY = limitsY, \
                     labelX = xLabel, labelY = yLabel, caption=self.ScanParam.caption(), \
                     style=style,ci=ci,overplot=overplot)
    else:
      Plot.draw(transpose(dataY), \
                sizeX=[min(dataX[0]),max(dataX[0])],sizeY=[1,len(chanList)],\
                limitsX=limitsX,\
                labelX=xLabel,labelY='Channel #', caption=self.ScanParam.caption(), wedge=1)
            

  #----------------------------------------------------------------------------
  def plotMeanChan(self,chanList=[], \
                   channelFlag=[], plotFlaggedChannels=0, \
                   dataFlag=[], plotFlaggedData=0, \
                   limitsX=[],limitsY=[], \
                   style='p', ci=1, overplot=0):
    """
    DES: PLotting the MEAN value for each subscan against channel number.
    """

    if plotFlaggedChannels:
        dataFlag = self._removeReservedFlagValues(dataFlag, self.rflags['CHANNEL FLAGGED'])
        if dataFlag==None:
            self.MessHand.error("no valid flags")
            return
        
    chanList = self.BolometerArray.checkChanList(chanList,
                                                 flag=channelFlag,
                                                 getFlagged=plotFlaggedChannels)
    if len(chanList)<1: 
        self.MessHand.error("no valid channel")
        return

    if not self.__statisticsDone:
        self.__statistics()

    self.MessHand.info("plotting Mean values per channel")
    dataY = self.getChanListData('mean',chanList,
                                 channelFlag=channelFlag,
                                 getFlaggedChannels=plotFlaggedChannels,
                                 dataFlag=dataFlag, getFlaggedData=plotFlaggedData)
    dataY = array(dataY)
    dataX = chanList
    Plot.plot(dataX,ravel(dataY),overplot=overplot,ci=ci,\
              limitsX=limitsX,limitsY=limitsY,\
              labelX='Channel Number',labelY='MEAN value for each subscan',\
              caption=self.ScanParam.caption(),)
        
  #----------------------------------------------------------------------------
  def plotRmsChan(self,chanList=[], \
                  channelFlag=[], plotFlaggedChannels=0, \
                  dataFlag=[], plotFlaggedData=0, \
                  limitsX=[],limitsY=[], \
                  style='p', ci=1, overplot=0, subscan = 0, logY=0):
    """
    DES: PLotting the RMS value for each subscan against channel number.
    INP: (logical) subscan: if 0, plot rms of the complete scan, if 1,
                            plot for each subscan and each channel
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
    """

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

    if not self.__statisticsDone:
      self.__statistics()

    if subscan:
      dataY = self.getChanListData('rms_s',chanList, \
                                   channelFlag=channelFlag, getFlaggedChannels=plotFlaggedChannels, \
                                   dataFlag=dataFlag, getFlaggedData=plotFlaggedData)
      labY = 'RMS value for each subscan'
      nbSubscan = shape(dataY)[1]
    else:
      dataY = self.getChanListData('rms',chanList, \
                                   channelFlag=channelFlag, getFlaggedChannels=plotFlaggedChannels, \
                                   dataFlag=dataFlag, getFlaggedData=plotFlaggedData)
      labY = 'RMS value for the complete scan'
      nbSubscan = 1
    self.MessHand.info("plotting r.m.s. per channel")

    arrayX = ones((nbSubscan),Float32)
    dataX = []
    for n in chanList:
      dataX.extend(n*arrayX)
    Plot.plot(dataX,ravel(dataY),overplot=overplot,ci=ci,\
              limitsX=limitsX,limitsY=limitsY,\
              labelX='Channel Number',labelY=labY,\
              caption=self.ScanParam.caption(),logY=logY)

  #----------------------------------------------------------------------------
  def plotFFT(self,chanList=[], \
              channelFlag=[], plotFlaggedChannels=0, \
              dataFlag=[], plotFlaggedData=0, \
              limitsX=[],limitsY=[], \
              style='l', ci=1, overplot=0, plot=1,logX=1,logY=1,\
              windowSize=0,windowing=3,returnSpectrum=0):
    """
    DES: plot FFT of signal
    INP: (i list) chanList : list of channels
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
         limits, style, ci...: plot parameters (see MultiPlot.plot)
    """
    
    chanList = self.BolometerArray.checkChanList(chanList, \
                                                 flag=channelFlag,getFlagged=plotFlaggedChannels)
    nChan = len(chanList)

    if plotFlaggedChannels:
        dataFlag = self._removeReservedFlagValues(dataFlag, self.rflags['CHANNEL FLAGGED'])
        if dataFlag==None:
            self.MessHand.error("no valid flags")
            return
        
    # We use only MJD values for fft. They are the correct time stamps.
    mjd = self.getChanData('mjd',chanList[0],flag='None')

    freqFFT = []
    modFFT  = []

    # TODO: compute by subscan
    for n in range(nChan):
        flux  = self.getChanData('flux',chanList[n],flag='None')
        num = self.BolometerArray.getChanIndex(chanList[n])[0]
        if not flux.iscontiguous():
            flux = flux.copy()

        if dataFlag in ['','None']:
            maskOk = nonzero(ones(shape=flux.shape, typecode='1'))
        else:
            if plotFlaggedData:
                maskOk = nonzero(self.FlagHandler.isSetMask(dataFlag, dim=1, index=num))
            else:
                maskOk = nonzero(self.FlagHandler.isUnsetMask(dataFlag, dim=1, index=num))

        theFluxOk = zeros(shape(flux),'f')
        put(theFluxOk,maskOk,take(flux,maskOk))

        oneFFT = FilterFFT(mjd,theFluxOk)
        oneFFT.doFFT(windowSize=windowSize,windowing=windowing)
        freq = copy.copy(oneFFT.Freq)
        dens = copy.copy(oneFFT.Power)

        # Remove the Zero frequency to use the log scale
        if logX:
            freq = freq[1::]
            dens = dens[1::]

        freqFFT.append(freq)
        # We plot the PSD (power) in amplitude units (that's V_rms / sqrt(Hz)
        modFFT.append(sqrt(dens))  # sqrt(Power)

    xLabel = "Frequency [Hz]"
    yLabel = "sqrt(PSD) [rms / sqrt(Hz)]"

    if plot:
        self.MessHand.info("plotting FFT(s)")
        
        MultiPlot.plot(chanList,freqFFT,modFFT,\
                       limitsX = limitsX, limitsY = limitsY,\
                       labelX = xLabel, labelY = yLabel, caption=self.ScanParam.caption(),\
                       style=style,ci=ci,overplot=overplot,logX=logX,logY=logY)
        
    if returnSpectrum:
        return(chanList,freqFFT,modFFT)


  #----------------------------------------------------------------------------
  def plotDataGram(self,chanNum=-1, \
                   flag=[], plotFlagged=0, \
                   n=512,limitsZ=[]):
    """
    DES: plot FFT of signal
    INP: (i)  chanNum : channel number to plot
         (integer list) flag : plot data flagged or unflagged accordingly
         (log)   plotFlagged : flag revers to flagged/unflagged data
                               flag   | plotFlagged | Plot..
                               'None' |  0          | all data
                               []     |  0          | unflagged data (default)
                               []     |  1          | data with at least one flag set
                               1      |  0          | data with flag 1 not set
                               1      |  1          | data with flag 1 set
                               [1,2]  |  0          | data with neither flag 1 nor flag 2 set
                               [1,2]  |  1          | data with either flag 1 or flag 2 set
         (i)        n : Number of points for the ffts
         (2f) limitsZ : limits for the color scale
    """

    if chanNum == -1:
        chanNum = self.BolometerArray.RefChannel
    time = self.getChanData('mjd',chanNum,flag=flag,getFlagged=plotFlagged)
    flux = self.getChanData('flux',chanNum,flag=flag,getFlagged=plotFlagged)

    oneFFT = FilterFFT(time,flux)
    oneFFT.plotDataGram(n=n,limitsZ=limitsZ)

  #----------------------------------------------------------------------------
  def bandRms(self,chanList=[],low=1.,high=10., \
              channelFlag=[], getFlaggedChannels=0, \
              dataFlag=[], getFlaggedData=0, \
              windowSize=0,windowing=3):
    """
    DES: compute rms in some spectral range
    INP: (i list) chanList : list of channels
         (f)     low, high : range limits (in Hz)
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
         (f)    windowSize : optional window size to compute FFTs
         (i)     windowing : function type for windowing (see applyWindow)
    """

    if getFlaggedChannels:
        dataFlag = self._removeReservedFlagValues(dataFlag, self.rflags['CHANNEL FLAGGED'])
        if dataFlag==None:
            self.MessHand.error("no valid flags")
            return
        
    chanList = self.BolometerArray.checkChanList(chanList, \
                                                 flag=channelFlag,getFlagged=getFlaggedChannels)
    
    # We use only MJD values for fft. They are the correct time stamps.
    time = self.getChanListData('mjd',chanList, \
                                 channelFlag=channelFlag, getFlaggedChannels=getFlaggedChannels, \
                                 dataFlag=dataFlag, getFlaggedData=getFlaggedData)
    flux = self.getChanListData('flux',chanList, \
                                 channelFlag=channelFlag, getFlaggedChannels=getFlaggedChannels, \
                                 dataFlag=dataFlag, getFlaggedData=getFlaggedData)
    
    for n in range(len(chanList)):      
        oneFFT = FilterFFT(time[n],flux[n])
        oneFFT.doFFT(windowSize=windowSize,windowing=windowing)
        mask = nonzero(greater(oneFFT.Freq,low) and less(oneFFT.Freq,high))
        total = 0.
        for k in mask:
            total += oneFFT.Power[k]
        total = sqrt(total * oneFFT.SamplFreq / oneFFT.N)
        self.MessHand.info("Chan %i: rms in [%f,%f Hz] = %g"%(chanList[n],low,high,total))
    



  #----------------------------------------------------------------------------
  # time shifting routines


  def timeShiftChan(self,chan,step,shiftFlags=1):
    """
    DES: time shift channel by step 
    INP: (i)       chan    : channel number
         (i)       step    : number of time stamps
         (bool) shiftFlags : also shift flags? default yes      
                  
    """
    
    chan          = self.BolometerArray.checkChanList(chan)
    chanListIndex = self.BolometerArray.getChanIndex(chan)
    
    nSamp=shape(self.Data)[0]
    
    if (step > 0):
        # shift data
        temp=copy.deepcopy(self.Data[0:nSamp-step,chanListIndex])       
        self.Data[step:nSamp,chanListIndex] = temp
        # shift flags
        if shiftFlags:
            temp=copy.deepcopy(self.FlagHandler._aFlags[0:nSamp-step,chanListIndex])
            self.FlagHandler._aFlags[step:nSamp,chanListIndex] = temp
        # flag beginning
        mask = zeros(nSamp, typecode='1')
        mask[0:step] = 1
        self.FlagHandler.setOnMask(mask,1,dim=1,index=chanListIndex)
    else:
        if (step < 0):
            # shift data
            self.Data[0:nSamp-abs(step),chanListIndex] = \
                     self.Data[abs(step):nSamp,chanListIndex]
            # shift flags
            if shiftFlags:
                self.FlagHandler._aFlags[0:nSamp-abs(step),chanListIndex] = \
                                         self.FlagHandler._aFlags[abs(step):nSamp,chanListIndex]
            # flag end
            mask = zeros(nSamp, typecode='1')
            mask[nSamp-abs(step):nSamp] = 1
            self.FlagHandler.setOnMask(mask,1,dim=1,index=chanListIndex)


  def timeShiftChanList(self,chanList,steps,shiftFlags=1):
    """
    DES: time shift list of channels by list of steps 
    INP: (i list)  chan : channel list
         (i list) steps : list of number of time stamps
         (bool) shiftFlags : also shift flags? default yes    
                  
    """

    chanList        = self.BolometerArray.checkChanList(chanList)
    chanListIndexes = self.BolometerArray.getChanIndex(chanList)
    
    for i in range(len(chanListIndexes)):
        self.timeShiftChan(chanList[i],steps[i],shiftFlags)


  def printListOfOddNames(self, extra=0):

      names=['womba womba','rolf','somrig ostsas']
      if extra:
          names.extend(['shanti roney'])

      self.MessHand.info('Printing list of odd names...')
      for n in names:
          print n
      self.MessHand.warning('note that some of these may not be real names, but names of sauces etc.')
  
        

  def computeCorTimeShift(self,shiftAz,shiftEl,chanList=[],refChan=-1,distlim=-1.):
    """
    DES: computes mean of absolute correlations for all channel pairs with mutual
         distance smaller than distlim, given time shifts in azimuth and elevation
         directions. To be used by timeshiftAzEl. 
    INP: (f)       shiftAz : time shift in azimuth. unit: milliseconds per arcsecond
         (f)       shiftEl : time shift in elevation. unit: milliseconds per arcsecond
         (i list) chanList : the list of channels to consider
         (i)       refChan : reference channel (timeshift 0)
         (f)       distlim : consider only correlations on bolometer separations
                             smaller than this value (arcseconds)
                  
    """  

    # create a copy of the data array
    chanList        = self.BolometerArray.checkChanList(chanList)
    chanListIndexes = self.BolometerArray.getChanIndex(chanList)
    refChanIndex = self.BolometerArray.getChanIndex(refChan)
    tempData=copy.deepcopy(self)

    # get the channel separations
    offAz = take(self.BolometerArray.Offsets[0,::],chanListIndexes) - \
            self.BolometerArray.RefOffX
    offEl = take(self.BolometerArray.Offsets[1,::],chanListIndexes) - \
            self.BolometerArray.RefOffY

    timestep=average(self.ScanParam.LST[1:20]-self.ScanParam.LST[0:19])*1000.

    # project channel separations onto direction of time shift
    shiftsteps=[]
    for i in range(len(offAz)):
        timeshift=dot(array([offAz[i],offEl[i]]),array([shiftAz,shiftEl]))
        if (timeshift > 0):
            shiftsteps.append(int(0.5+timeshift/timestep))
        else:
            shiftsteps.append(int(timeshift/timestep-0.5))

    tempData.timeShiftChanList(chanList,shiftsteps)
        
    # compute correlations
    #tempData.__resetStatistics()
    #tempData.computeWeights(chanList=[])
    dist, corr = tempData.plotCorDist(chanList=chanList,upperlim=distlim,plot=0)

    # return mean of absolute values of the correlations
    return average(abs(corr))


  
  def timeshiftAzEl(self,chanList=[],refChan=-1,check=1,distlim=300.,shiftmax=10.):
    """
    DES: computes time shifts of all channels, with respect to a reference
         channel, which MAXIMIZES the correlated noise across the array
    
    INP: (i list) chanList : the list of channels to consider
         (i)       refChan : reference channel (will get timeshift=0)
         (l)         check : check the chanList first ( default: yes )
         (f)       distlim : consider only correlations on bolometer separations
                             smaller than this value (arcseconds)
         (f)      shiftmax : maximum timeshift (absolute value) in milliseconds per arcsecond
         
    """

    chanList        = self.BolometerArray.checkChanList(chanList)
