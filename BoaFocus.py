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
NAM: BoaFocus.py (module)
DES: contains the BoA focus class
"""
__version__=  '$Revision: 2794 $'
__date__=     '$Date: 2015-03-02 15:30:03 +0100 (Mon, 02 Mar 2015) $'

# ---------------------------------------------------------------------
# ---- Import ---------------------------------------------------------
# ---------------------------------------------------------------------
from Numeric import *
from boa           import BoaPointing
from boa.Bogli     import Plot, MultiPlot, BogliConfig
from boa.Utilities import fitParabola, modelparabola
from boa.fortran   import fStat

# ---------------------------------------------------------------------
# ---- Boa Focus Class ------------------------------------------------
# ---------------------------------------------------------------------

class Focus(BoaPointing.Point):
  """
  NAM: Focus (class)
  DES: An object of this class is responsible for the focus reduction
       of single or multiple scans and provides the offsets.
  """

  def __init__(self):
    """
    DES: Initialise an instance
    """
    BoaPointing.Point.__init__(self)
    self.FocusOffset = 0.
    self.FocusError  = 0.

  # -------------------------------------------------------------------
  # ---- public methods -----------------------------------------------
  # -------------------------------------------------------------------
  def solveFocus(self,noerase=0,caption=''):
    """
    DES: compute the optimal focus position
    """

    refChan = self.BolometerArray.RefChannel
    modType = self.ScanParam.ScanType
    nSubscan = len(self.ScanParam.SubscanNum)  # subscans read in

    iflux  = zeros((nSubscan),Float32)
    isdev  = zeros((nSubscan),Float32)
    ifocus = zeros((nSubscan),Float32)

    for i in range(nSubscan):
      flux    = self.getChanData('flux',refChan,subscans=[self.ScanParam.SubscanNum[i]])
      focus   = self.getChanData(modType,refChan,subscans=[self.ScanParam.SubscanNum[i]])
      iflux[i] = fStat.f_median(flux)
      isdev[i] = fStat.f_rms(flux,iflux[i])
      ifocus[i] = fStat.f_median(focus)
      
    if string.find(modType.upper(),'TILT') > -1:
      ifocus = ifocus*3600.
      xLabel = modType +' Offsets [arcsec]'
    else:
      xLabel = modType +' Offsets [mm]'
    yLabel = 'flux density [arb u.] Channel :'+str(refChan)

    result = fitParabola(ifocus,iflux,isdev)

    x = arange(101.)/100.*(max(ifocus)-min(ifocus))+min(ifocus)
    y = modelparabola(result.params,x)
    
    if not caption:
      caption = self.ScanParam.caption()
    Plot.plot(ifocus,iflux,style='p',\
              labelX=xLabel,labelY=yLabel,\
              caption=caption,noerase=noerase)
    
    # TODO : need plot with error bar.
    #Plot.plot(ifocus,iflux,style='p',ci=1,overplot=1)
    Plot.plot(x,y,overplot=1,style='l')
    if result.params[2]:
      self.FocusOffset = result.params[1]/(-2.0*result.params[2])
      if result.perror:
        self.FocusError  = result.perror[1]/(-2.0*result.params[2])
      else:
        self.FocusError  = 0.
      self.MessHand.info("Offset : %5.2f +- %5.2f" \
                         %(self.FocusOffset,self.FocusError))
      self.MessHand.longinfo("Parabola: a + b x + c x^2, "+
                          "a = %7.4f b = %7.4f c = %7.4f"
                          %(result.params[0],result.params[1],result.params[2]))
    else:
      self.MessHand.warning("no focus offset found")
      self.FocusOffset = 0.
      self.FocusError  = 0.

  def reduce(self,datasetName='',obstoProc=[],febe='',baseband=1):
    """
    DES: Process a Focus scan - this method is called by the apexCalibrator
    INP: (string) datasetName: path to the dataset to be reduced
         (i list)   obstoProc: list of subscans to consider (default: all)
    """
    if len(obstoProc)==1:
      if type(obstoProc[0]) == type([]): # e.g. obstoProc == [range(4,8)]
        self.read(inFile=datasetName,subscans=obstoProc[0],
                  febe=febe,baseband=baseband)
      else:
        # cannot work subscan by subscan
        self.read(inFile=datasetName,subscans=range(1,obstoProc[0]+1),
                  febe=febe,baseband=baseband)
    else:
      self.read(inFile=datasetName,subscans=obstoProc,
                febe=febe,baseband=baseband)

    # If chopped data, then compute phase diff
    if self.ScanParam.WobUsed:
      self.ScanParam.computeOnOff()
      self._phaseDiff()
            
    self.zeroStart(subscan=0)

    Plot.panels(2,2)
    Plot.nextpage()  # start a new page
    textColBak = BogliConfig.xyouttext['color']
    BogliConfig.xyouttext['color'] = 3
    # First plot = signal before skynoise removal
    four = self.BolometerArray.fourpixels()
    self.signal(four,noerase=1,caption=self.ScanParam.caption()+' - Raw signal')
    Plot.nextpage()
    self.solveFocus(noerase=1,caption=' ')
    y1,y2 = Plot.yAxis['limits']
    Plot.xyout(0.,(y1+y2)/2.,
               str("%5.2f +- %4.2f"%(self.FocusOffset,self.FocusError)),size=2.5)
    Plot.plot([self.FocusOffset,self.FocusOffset],y1+array([0.7,1.],'f')*(y2-y1),
         overplot=1,style='l',ci=2)
    x1,x2 = Plot.xAxis['limits']
    Plot.xyout(x1+0.85*(x2-x1),y1+0.9*(y2-y1),'Raw data',size=2)
    
    # again after skynoise removal
    self.medianNoiseRemoval(computeFF=0,chanRef=-1,nbloop=3)
    Plot.nextpage()
    self.signal(four,caption='Skynoise subtracted',noerase=1)
    Plot.nextpage()
    self.solveFocus(noerase=1,caption=' ')
    y1,y2 = Plot.yAxis['limits']
    Plot.xyout(0.,(y1+y2)/2.,
               str("%5.2f +- %4.2f"%(self.FocusOffset,self.FocusError)),size=2.5)
    Plot.plot([self.FocusOffset,self.FocusOffset],y1+array([0.7,1.],'f')*(y2-y1),
         overplot=1,style='l',ci=2)
    modType = self.ScanParam.ScanType
    Plot.xyout(0.,y1+0.2*(y2-y1),modType,size=3)
    x1,x2 = Plot.xAxis['limits']
    Plot.xyout(x1+0.82*(x2-x1),y1+0.91*(y2-y1),'Skynoise',size=2)
    Plot.xyout(x1+0.82*(x2-x1),y1+0.84*(y2-y1),'subtracted',size=2)
    Plot.panels(1,1)
    BogliConfig.xyouttext['color'] = textColBak
