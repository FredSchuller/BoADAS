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
NAM: BoaSkydip.py (module)
DES: contains the BoA skydip class
"""
__version__=  '$Revision: 1726 $'
__date__=     '$Date: 2006-09-25 19:38:54 +0200 (Mon, 25 Sep 2006) $'

# ---------------------------------------------------------------------
# ---- Import ---------------------------------------------------------
# ---------------------------------------------------------------------
from Numeric import *
from boa           import BoaPointing, BoaConfig
from boa.Bogli     import Plot, MultiPlot
from boa.Utilities import fitParabola, modelparabola
from boa.fortran   import fStat,fFit

# ---------------------------------------------------------------------
# ---- Boa Skydip Class ------------------------------------------------
# ---------------------------------------------------------------------

class Skydip(BoaPointing.Point):
  """
  NAM: Skydip (class)
  DES: An object of this class is responsible for the reduction
       of skydips and provides zenithal opacity
  """

  def __init__(self):
    """
    DES: Initialise an instance
    """
    BoaPointing.Point.__init__(self)
    # The following two attributes are for communication with apexCalibrator
    self.ResultAirmass = []
    self.ResultFluxes  = []
    # Calibration V to K - this has  been measured for LABOCA
    self.VtoK     = 1.E6 / 10.6
    # These parameters come into the skydip equation - default are reasonable values
    self.TCabin   = 283.  # +10 deg C
    self.TAtm     = 270.  # -3 deg C
    self.Coupling = 0.8   # 80% optical coupling
    self.Feff     = 0.97  # APEX forward efficiency at 350 GHz
    self.Tau      = 0.3   # zenith opacity, that's what we want to determine

  # -------------------------------------------------------------------
  # ---- public methods -----------------------------------------------
  # -------------------------------------------------------------------
  def solveSkydip(self):
    """
    DES: compute the zenithal opacity
    """

    refChan = self.BolometerArray.RefChannel
    flux    = self.getChanData('flux',refChan)

    nSubscan = self.ScanParam.NObs
    subscanIndex = self.ScanParam.SubscanIndex

    self.MessHand.warning("Coming soon...")
    return
  
  def reduce(self,datasetName='',obstoProc=[],blind=317):
    """
    DES: Process a skydip scan - this method is called by the apexCalibrator
    INP: (string) datasetName: path to the dataset to be reduced
         (i list)   obstoProc: list of subscans to consider (default: all)
         (i)            blind: channel number of the ref. (blind) bolo
                           !!! SPECIFIC TO LABOCA !!!
    """
    self.read(inFile=datasetName,subscans=obstoProc)
            
    self.medianBaseline(subscan=0)
    # Continuous mode: try to find subscans
    if self.ScanParam.NObs == 1:
      self.findElSubscan()
    #refChan = self.BolometerArray.RefChannel
    # for now, use channel 110
    refChan = 110
    flux = self.getChanData('flux',refChan)
    el = self.getChanData('el')

    # Number of subscans read (or found):
    nObs = len(self.ScanParam.SubscanNum)
    resultAir  = []
    resultFlux = []
    for i in range(nObs):
      ind1 = self.ScanParam.SubscanIndex[0,i]
      ind2 = self.ScanParam.SubscanIndex[1,i]
      airmass = 1./sin(el[ind1:ind2]*pi/180.)
      resultAir.append(fStat.f_median(airmass))
      # Compute difference chan - gain*blind
      if blind:
        refFlux = self.getChanData('flux',blind)
        # Fit the relative gain
        f1 = flux[ind1:ind2]
        f2 = refFlux[ind1:ind2]
        #relG = fFit.fit_slope(f2,f1)
        #print "subscan %i - relG = %f"%(i,relG)
        #ratio = fStat.f_median(relG)
        #self.MessHand.info('Computed relative gain %i/%i = %f'%(blind,refChan,ratio))
        ratio = 0.7
        diff = f1 - ratio*f2
      else:
        refFlux = 0
        diff = flux[ind1:ind2] - flux[0]
      
      resultFlux.append(fStat.f_median(diff))
      #resultFlux.append(fStat.f_median(diff[ind1:ind2]))
      
    self.ResultAirmass = resultAir
    self.ResultFluxes = resultFlux

  def findElAccSubscan(self):
    """
    determine subscans by looking at elevation acceleration
    """
    ae = self.getChanData('elacc',flag='None')  # get ALL in order to compute subscan indices
    ind1 = [0]
    ind2 = []
    nb = len(ae)
    i = 0
    while i < nb-3:
      while((ae[i] > -2000 or ae[i+1] > -2000 or ae[i+2] > -2000) and i < nb-3):
        i += 1
      ind2.append(i)
      while((ae[i] < 2000 or ae[i+1] < 2000 or ae[i+2] < 2000) and i < nb-3):
        i += 1
      while((ae[i] > 2000 or ae[i+1] > 2000 or ae[i+2] > 2000) and i < nb-3):
        i += 1
      ind1.append(i)
    ind1 = ind1[:-1]

    # Reinitialise Subscan related attributes
    self.ScanParam.SubscanIndex = []
    self.ScanParam.SubscanNum   = []
    self.ScanParam.SubscanType  = []
    nb = len(ind1)
    for i in range(nb):
      self.ScanParam.SubscanIndex.append([ind1[i],ind2[i]])
      self.ScanParam.SubscanNum.append(i+1)
      self.ScanParam.SubscanType.append('ON')

    self.ScanParam.SubscanIndex = transpose(self.ScanParam.SubscanIndex)
    self.ScanParam.NObs = nb
    self.MessHand.info(str("Found %i subscans"%(nb)))


  def findElSubscan(self,tolerance=0.2,delta=5.e-4):
    """
    DES: Determine subscans by looking at elevation
    INP: (f) tolerance: tolerance (in deg) to assume that we're still in the same subscan
         (f)    delta : minimum variation of El (in deg) for the steps between subscans
    """
    el = self.getChanData('el')
    dEl = abs(el[1::] - el[:-1])
    ind1 = []
    ind2 = []
    i = 0
    nb = len(dEl)
    while i < nb-1:
      # The following is true when between two subscans
      while dEl[i] > delta and i < nb-1:
        i += 1
      currEl = el[i]
      ind1.append(i)
      # The following stays true during one subscan
      while abs(el[i]-currEl) < tolerance and i < nb-1:
        i += 1
      ind2.append(i)
    
    # Reinitialise Subscan related attributes
    self.ScanParam.SubscanIndex = []
    self.ScanParam.SubscanNum   = []
    self.ScanParam.SubscanType  = []
    nb = len(ind1)
    for i in range(nb):
      self.ScanParam.SubscanIndex.append([ind1[i],ind2[i]])
      self.ScanParam.SubscanNum.append(i+1)
      self.ScanParam.SubscanType.append('ON')

    self.ScanParam.SubscanIndex = transpose(self.ScanParam.SubscanIndex)
    self.ScanParam.NObs = nb
    self.MessHand.info(str("Found %i subscans"%(nb)))
      
