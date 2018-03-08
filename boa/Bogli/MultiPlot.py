# Copyright (C) 2002-2010
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
NAM: BogliMultiChanPlot.py (module)
DES: contains Bogli multi-channel plot class 
"""
__version__=  '$Revision: 2722 $'
__date__=     '$Date: 2010-04-29 11:43:32 +0200 (Thu, 29 Apr 2010) $'

#------------------------------------------------------------------------------
#----- Import -----------------------------------------------------------------
#------------------------------------------------------------------------------
from ppgplot import *
from Numeric import *
import Plot
from boa import BoaMessageHandler
from boa.Utilities import Timing
from boa.fortran.fUtilities import mreplacenan, masknan

from boa.fortran.fStat import minmax

import DeviceHandler

#------------------------------------------------------------------------------
#----- MultiChan --------------------------------------------------------------
#------------------------------------------------------------------------------

# Description of the Channel Number, same definition as Plot.__labelX
chanNum={'charheight'     : 1.0,  \
         'color'          : 2,    \
         'displacement'   : -1.5, \
         'coordinate'     : 0.1, \
         'justification'  : 0.0,  \
         'side'           : 't',  \
         'text'           : ' ',  \
         'draw_it'        : 1}


# The individuals viewPoint size, saved for overplot
viewPointStep = {'X':1.0, 'Y': 1.0} 

# The private message Handler
__MessHand = BoaMessageHandler.MessHand('MultiPlot')

#----------------------------------------------------------------------------
#----- Methods --------------------------------------------------------------
#----------------------------------------------------------------------------

def detSubDivView(numPlot,wedge=0):
  """
  DES: Determine sub-division of plot page into viewpoints/boxes.
  INP: (int) numPlot   : number of panels in the plot
  OPT: (logical) wedge : if a wedge is present (no by default) 
  OUT: (ints) (numPlotWinX,numPlotWinY) the number of viewport in both directions
  """

  global viewPointStep, chanNum

  # Divide in a square
  numPlotWinX = int(floor(sqrt(numPlot)))
  numPlotWinY = int(floor(numPlot/numPlotWinX))
  
  # In case the is more window than a perfect square
  if numPlotWinX*numPlotWinY < numPlot:
    numPlotWinY += 1

  if wedge == 1 or Plot.wasWedge == 1:
    # If wedge is present reserve space for it
    viewPointStep['X'] = (Plot.globalViewPoint['marging_right'] - \
                          Plot.globalViewPoint['marging_left'] - \
                          Plot.globalViewPoint['wedge_size'])/ \
                          numPlotWinX
  else:
    viewPointStep['X'] = (Plot.globalViewPoint['marging_right'] - \
                          Plot.globalViewPoint['marging_left'])/ \
                          numPlotWinX 

  viewPointStep['Y'] = (Plot.globalViewPoint['marging_top'] - \
                        Plot.globalViewPoint['marging_bottom'])/ \
                        numPlotWinY

    
  # Reserve room for the wedge if needed
  if wedge:
    viewPointStep['X'] -= Plot.globalViewPoint['wedge_size']/numPlotWinX

  # Change the character height of the plot according to the
  # computed viewPoint size
#  box['charheight']     = 1.0/sqrt(0.7*numPlotWinY)
  chanNum['charheight'] = 1.0/sqrt(0.3*numPlotWinY)
  
  # TODO: somehow the Point['size'] should be change
  # Point['size'] = 4.0/sqrt(0.7*numPlotWinY)

  return (numPlotWinX,numPlotWinY)

#----------------------------------------------------------------------------
def plotBox(numPlot,x,y):
  """
  DES: Draw box and labels.
  INP: (int) numPlot : number of panels in the plot
       (int) x,y     : the indices of the panel
  """

  myTiming = Timing()
  
  pgsch(Plot.box['charheight']/sqrt(0.5*sqrt(numPlot)))                            # set character height 
  pgsci(Plot.box['color'])                                 # set colour index 
  pgsls(Plot.box['linestyle'])                             # set line style
  pgslw(Plot.box['linewidth'])                             # set line width 
  lastRow=0

  if (1+y)*floor(sqrt(numPlot)) >= numPlot:
    lastRow=1

  if Plot.xAxis['log']:
    codeX = 'BCTSL'
  else:
    codeX = 'BCTS'
  if Plot.yAxis['log']:
    codeY = 'BCTSL'
  else:
    codeY = 'BCTS'
  
  if x==0 and lastRow:
    pgbox(codeX+'N',0.0,0,codeY+'NV',0.0,0)              # corner: label x and y
  elif x==0:
    pgbox(codeX,0.0,0,codeY+'NV',0.0,0)                  # label y
  elif lastRow:
    pgbox(codeX+'N',0.0,0,codeY,0.0,0)                   # label x 
  else:
    pgbox(codeX,0.0,0,codeY,0.0,0)                       # label nothing

  __MessHand.debug("plotBox done in "+str(myTiming))

#--------------------------------------------------------------------------------
def setLimits(dataX,dataY,limitsX=[],limitsY=[]):
  """
  DES: compute and/or set the limits for the multiplot
  INP: (list arrays)   dataX/Y : the array to be plotted
  OPT: (2elts array) limitsX/Y : limits to use in X/Y for the plot
  """

  myTiming = Timing()
  
  xLimits = []
  yLimits = []
  
  # compute the limits for each plot and ... 
  for i in range(len(dataX)):
    if len(dataX[i])>0 :
      # Make sure that X and Y have the same length
      if len(dataX[i]) > len(dataY[i]):
        dataX[i] = dataX[i][:len(dataY[i])]
      elif len(dataY[i]) > len(dataX[i]):
        dataY[i] = dataY[i][:len(dataX[i])]
        
      # Remove all NaN from X and Y arrays
      maskX = masknan(dataX[i])
      maskY = masknan(dataY[i])
      mNan = bitwise_or(maskX,maskY)
      bad = nonzero(not_equal(mNan,0))
      if len(bad):
        # at least one NaN somewhere
        __MessHand.warning("%i data points NaN were found"%(len(bad)))
        good = nonzero(equal(mNan,0))
        if len(good):
          dataX[i] = take(dataX[i],good)
          dataY[i] = take(dataY[i],good)
        else:
          dataX[i] = array([0],'f')
          dataY[i] = array([0],'f')      

      Plot.setLimits(dataX[i],dataY[i],limitsX=limitsX,limitsY=limitsY)
      xLimits.append(Plot.xAxis['limits'])
      yLimits.append(Plot.yAxis['limits'])

  xLimits = ravel(xLimits)
  yLimits = ravel(yLimits)

  # .. take the maximum of all
  Plot.xAxis['limits'] = minmax(xLimits)
  Plot.yAxis['limits'] = minmax(yLimits)

  __MessHand.debug("setLimits time: "+myTiming.__str__())
  

#--------------------------------------------------------------------------------
def drawChanNum(c):
  """
  DES: Draw channel number.
  INT: (int) c : the channel number
  """

  global chanNum

  chanNum['text'] = str(c)
  Plot.drawLabel(chanNum)

#--------------------------------------------------------------------------------
def setMultiViewPoint( x=0, y=0, wedge=0):
  """
  DES: Determine and set view points.
  INP: (ints)    x/y      : position of the viewpoint from 0 to numPlotWinX/Y
  """

  global viewPointStep

  xleftIndVP  = Plot.globalViewPoint['marging_left'] +  x*viewPointStep['X']
  xrightIndVP = Plot.globalViewPoint['marging_left'] + (x+1.0)*viewPointStep['X']
  ytopIndVP   = Plot.globalViewPoint['marging_top']  -  y*viewPointStep['Y']
  ybotIndVP   = Plot.globalViewPoint['marging_top']  - (y+1.0)*viewPointStep['Y']

  # Set the viewpoint...
  pgsvp(xleftIndVP, xrightIndVP, ybotIndVP, ytopIndVP)
  
  # ... and set the limits for the plots
  pgswin(Plot.xAxis['limits'][0],Plot.xAxis['limits'][1],\
         Plot.yAxis['limits'][0],Plot.yAxis['limits'][1])

  # TODO : Do we need an option for aspect ratio in multichan plots ?
#--------------------------------------------------------------------------------
def gloLabelling(wedge=0):
  """
  DES: Label x, y, caption and channel number.
  OPT: (logical) wedge : if a wedge is present
  """


  # Use the ViewPoint as if there were only one plot
  Plot.setViewPoint(wedge=wedge)
  # And label it.
  Plot.labelling(wedge=wedge)

#--------------------------------------------------------------------------------

def plot(chanList,dataX,dataY, \
         limitsX=[], limitsY=[], \
         labelX='x',labelY = 'y',caption=' ',\
         style='p', ci=1, overplot=0, logX=0, logY=0,
         nan=0, noerase=0):
  """
  DES: do a multi channel plot
  INP: (int list) chanList = list of channels
       (array list)  dataX = values to plot along X
       (array list)  dataY = values to plot along Y
  OPT: (2elts array) limitsX/Y = limits to use in X/Y for the plot
       (string)      labelX  = x label (default 'x')
       (string)      labelY  = y label (default 'y')
       (string)      caption = the caption of the plot (default ' ')
       (char)        style   = the style used for the plot ('l': line,
                               'p': point (default), 'b': histogram)
       (int)            ci   = color index (default 1)
       (logical)    overplot = are we overplotting ? (default no)
       (logical)      logX/Y = do we use log scale ? (default no)
       (logical)     noerase = do not clear the window
  """
  myTiming = Timing()

  # Check for an open device first
  if not DeviceHandler.CurrentDev:
    __MessHand.warning('no device open')
    DeviceHandler.openDev()
  
  if overplot:
    if len(chanList) != Plot.numPlot:
      __MessHand.error('cannot overplot, wrong number of channels')
      return
  else:
    # If we are not overplotting clear the graph and set the limits/label
    if not noerase:
      Plot.clear()
    if logX:
      Plot.xAxis['log'] = 1
    if logY:
      Plot.yAxis['log'] = 1
    setLimits(dataX,dataY,limitsX=limitsX,limitsY=limitsY)
    Plot.setLabels(labelX=labelX,labelY=labelY,caption=caption)
    Plot.numPlot = len(chanList)

  #Divide the screen properly
  (numPlotWinX,numPlotWinY) = detSubDivView(Plot.numPlot)

  # pgbbuf() buffers the drawing on the screen until the pgebuf() command.
  # There is a gain of up to 10% in time when using this on a
  # multichanplot of 400 channels on full screen resolution.

  pgbbuf() 

  # Now do the plots
  for y in range(numPlotWinY): 
    for x in range(numPlotWinX):

      # counter of plotted channel
      cc = y*numPlotWinX+x

      # terminate when last chan to plot is reached
      if cc == Plot.numPlot:
        break
      
      # determine and set view point
      setMultiViewPoint(x,y)    

      if not overplot and Plot.box['draw_it']:   # draw box and labels
        plotBox(Plot.numPlot,x,y)

      # Do the plot
      Plot.plotDataXY(dataX[cc],dataY[cc],style=style,ci=ci)

      if not overplot and chanNum['draw_it']: # draw the channel number
         drawChanNum(chanList[cc])
      
  if not overplot:
    # do not forget the label everything if we are not overplotting
    gloLabelling()

  pgebuf()

  __MessHand.debug("plot done in "+str(myTiming))

#--------------------------------------------------------------------------------

def draw(chanList,map_arrays,sizeX=[], sizeY=[], WCS = [], \
         limitsX=[], limitsY=[], limitsZ=[], \
         nan = 0,\
         labelX='x',labelY = 'y',caption=' ',\
         style='g2r', contrast=1.0, brightness=0.5, \
         wedge=0, overplot=0 ):
  
  """
  DES: do a multi channel image drawing
  INP: (int list) chanList = list of channels
       (map_arrays)  lits of map to display
  OPT: (2elts arrays)   sizeX/Y = the 'physical' size of the array (default pixel numbers)
       (2elts arrays) limitsX/Y = limits to use in X/Y for the plot
       (string)         labelX  = x label (default 'x')
       (string)         labelY  = y label (default 'y')
       (string)         caption = the caption of the plot (default ' ')
       (char)           style   = the color used for the plot (default 'g2r'
                                  see Plot.Plot.setImaCol())
       (logical)        wedge   = shall we draw a wedge ? (default no)
       (logical)       overplot = are we overplotting ? (default no)
  """

  myTiming = Timing()

  # Check for an open device first
  if not DeviceHandler.CurrentDev:
    __MessHand.warning('no device open')
    DeviceHandler.openDev()

  lmap_arrays = array(copy.copy(map_arrays))

  # points in the map that were not observed are at NaN,
  # retrieved them... ugly
  if nan:
    new = mreplacenan(lmap_arrays)
    for i in range(len(lmap_arrays)):
      lmap_arrays[i] = new[i]

  # If the size of the map have not been given use a integer
  # range, use the lmap_arrays[0], i.e. all the maps have to
  # have the same dimension for the moment
  # TODO: fix that
  
  if sizeX == []:
    sizeX = array([0,lmap_arrays[0].shape[0]-1],Float)
  else:
    # Just to be sure that sizeX is an array of Float
    sizeX = array(sizeX,Float)
    
  # Same for sizeY
  if sizeY == []:
    sizeY = array([0,lmap_arrays[0].shape[1]-1],Float)
  else:
    sizeY = array(sizeY,Float)

  # Compute the transformation matrix for this image even in case of overplot
  
  Plot.setMapTransformation(lmap_arrays, sizeX=sizeX, sizeY=sizeY, WCS=WCS)

  if not overplot:
    # Since we are note overplotting, clear the plot, set
    # size/limits/labels/color
    
    Plot.clear()

    Plot.setMapLimits(lmap_arrays, limitsX=limitsX, limitsY=limitsY, limitsZ=limitsZ)
    Plot.setLabels(labelX,labelY,caption)
    Plot.setImaCol(style=style, contrast=contrast, brightness=brightness)
    Plot.numPlot = len(chanList)
    # Keep in min that they were a wedge
    Plot.wasWedge = wedge

  if overplot:
    if len(chanList) != numPlot:
      __MessHand.error('cannot overplot, wrong number of channels')
      return

  #Divide the screen properly
  (numPlotWinX,numPlotWinY) = detSubDivView(Plot.numPlot,wedge=wedge)

  # pgbbuf() buffers the drawing on the screen til the pgebuf()
  # command.
  pgbbuf() 

  # Now do the plots
  for y in range(numPlotWinY): 
    for x in range(numPlotWinX):

      # Channel counter
      cc = y*numPlotWinX+x

      # terminate when last chan to plot is reached
      if cc == Plot.numPlot:
        break
      
      # determine and set view point
      setMultiViewPoint(x,y)    

      if Plot.box['draw_it']:   # draw box and labels
        plotBox(Plot.numPlot,x,y)

      # Draw the image
      pgimag(transpose(lmap_arrays[cc]),lmap_arrays[cc].shape[0],lmap_arrays[cc].shape[1],\
             0, lmap_arrays[cc].shape[0]-1,0,lmap_arrays[cc].shape[1]-1,\
             Plot.zAxis['limits'][0],Plot.zAxis['limits'][1],\
             Plot.transformationMatrix)


      if chanNum['draw_it']: # draw the channel number
         drawChanNum(chanList[cc])  # 0=black
      
  if not overplot:
    # Do not forget to label everything with a wedge if needed
    gloLabelling(wedge=wedge)

  pgebuf()
    
  __MessHand.debug("draw done in "+str(myTiming))


## @}
