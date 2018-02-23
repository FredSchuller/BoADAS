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
NAM: BogliPlot.py (file)
DES: contains Bogli command handler class
"""
__version__=  '$Revision: 2722 $'
__date__=     '$Date: 2010-04-29 11:43:32 +0200 (Thu, 29 Apr 2010) $'
#------------------------------------------x----------------------------------------
#----- Import ---------------------------------------------------------------------
#----------------------------------------------------------------------------------
from ppgplot import *
from Numeric import *
import copy, cPickle
from boa import BoaMessageHandler
from boa.Utilities import solvePoly, Timing

from BogliConfig import *
import DeviceHandler
from boa.fortran.fUtilities import replacenan, masknan
from boa.fortran.fStat import minmax

# Keep in mind if a wedge was plotted or not (used to overplot on
# image) or that the aspect keyword was used

wasWedge  = 0
wasAspect = 0

# The number of plot to do (used to check for overplot)
numPlot = 1

# Used for maps : 
transformationMatrix = array([0,1,0,0,0,1], Float)

# Add a dedicated message handler as attribute
__MessHand = BoaMessageHandler.MessHand('Plot')


#--------------------------------------------------------------------------------
#----- methods ------------------------------------------------------------------
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
def erase():
  """
  DES: Erase any existing graphics.
  """
  pgeras()
  
#--------------------------------------------------------------------------------
def clear():
  """
  DES: clear plot
  """

  global __MessHand, xAxis, yAxis, wasWedge, wasAspect

  # Do we really need this ? Couldn't erase be enough ?
  current_device = pgqid()
  __MessHand.debug("clearing plot of device "+ `current_device`)
  pgeras()

  # Resize the graph so that the resize command is not needed anymore only for /xwin
  if not pgqinf('dev/type'):
    pgpage()

  # reinitialise log-scale flags - FS050515
  # TODO Do we need that ? if yes all the other attribute should be also set 
  xAxis['log'] = 0
  yAxis['log'] = 0

  # reinitialise the wedge and aspect memory
  wasWedge = 0
  wasAspect = 0
  
#--------------------------------------------------------------------------------
def setViewPoint(wedge=0,aspect=0):
  """
  DES: set the  view point
  OPT: (logical) wegde  : if wedge is present (default no)
       (logical) aspect : should we keep the aspect ratio (in 'physical' unit)
                          of the graph (default no)
  """

  global wasWedge, wasAspect, globalViewPoint, xAxis, yAxis
  
  # set the viewPoint
  if wedge == 1 or wasWedge == 1:
    # If wedge is present reserve space for it
    pgsvp(globalViewPoint['marging_left'], \
          globalViewPoint['marging_right']-globalViewPoint['wedge_size'], \
          globalViewPoint['marging_bottom'],globalViewPoint['marging_top'])
  else: 
    pgsvp(globalViewPoint['marging_left'], globalViewPoint['marging_right'], \
          globalViewPoint['marging_bottom'], globalViewPoint['marging_top'])

  if aspect == 0 and wasAspect == 0:
    # and set the limits for the plots
    pgswin(xAxis['limits'][0],xAxis['limits'][1],\
           yAxis['limits'][0],yAxis['limits'][1])
  else:
    pgwnad(xAxis['limits'][0],xAxis['limits'][1],\
           yAxis['limits'][0],yAxis['limits'][1])
    
#--------------------------------------------------------------------------------
def setLabels(labelX, labelY, caption):
  """
  DES: check and set labels 
  INP: (strings) labelX, labelX, caption : the X and Y labels and the caption
  """

  global xAxis, yAxis, cLabel
  
  if isinstance(labelX,str):
    xAxis['text'] = labelX
  else:
    xAxis['text'] = 'x'
    __MessHand.warning("invalid x-label: "+`labelX`)

  if isinstance(labelY,str):
    yAxis['text'] = labelY
  else:
    yAxis['text'] = 'y'
    __MessHand.warning("invalid y-label: "+`labelY`)

  if isinstance(caption,str):
    cLabel['text'] = caption
  else:
    cLabel['text'] = ''
    __MessHand.warning("invalid caption: "+`caption`)

#--------------------------------------------------------------------------------
def labelling( wedge=0):
  """
  DES: Label x, y, caption and channel number.
  OPT: (logical) wedge : should we draw a wedge ? (default no)
  """

  global xAxis, yAxis, cLabel
  
  if xAxis['draw_it']: drawLabel(xAxis)
  if yAxis['draw_it']: drawLabel(yAxis)
  if cLabel['draw_it']: drawLabel(cLabel)

  if wedge:
    drawWedge()

#--------------------------------------------------------------------------------
def drawWedge():
  """
  DES: Draw a wedge.
  """

  global zAxis

  # Maintained as a separate function but could be merged with
  # BogliPlot.Plot.labelling()...
  # TODO : think about that and fix
  
  if zAxis['draw_it']: 
    
    pgsch(zAxis['charheight']) # set character height
    pgsci(zAxis['color'])      # set colour index
    pgsls(zAxis['linestyle'])  # set line style
    pgslw(zAxis['linewidth'])  # set line width
    pgwedg(zAxis['side'],\
           zAxis['displacement'],\
           zAxis['width'],\
           zAxis['limits'][1],\
           zAxis['limits'][0],\
           zAxis['text'])

#--------------------------------------------------------------------------------
def drawLabel(label):
  """
  DES: generic function to draw labels/captions
  INP: (dict) label : a 'label' attribute like labelX
  """
  
  pgsch(label['charheight']) # set character height
  pgsci(label['color'])      # set colour index
  pgmtxt(label['side'],\
         label['displacement'],\
         label['coordinate'],\
         label['justification'],\
         label['text'])

#--------------------------------------------------------------------------------
def xyout(X,Y,text,size=0):
  """
  DES: generic function to overplot text
  INP: (int) X,Y : pixel position
       (str) text: text to output
       (flt) size: charheight (optional, default = use BogliConfig)
  """

  if size:
    pgsch(size)
  else:
    pgsch(xyouttext['charheight']) # set character height
  pgsci(xyouttext['color'])      # set colour index
  pgptxt(X,Y,\
         xyouttext['angle'],\
         xyouttext['justification'],\
         text)


#--------------------------------------------------------------------------------
def plotDataXY(dataX,dataY,style='p',ci=1,width=0,ls=1):
  """
  DES: Plot x y data.
  INP: (array) dataX/Y : the array to be plotted (same dimension)
  OPT: (string) style : the style used for the plot ('l': line,
                        'p': point (default), 'b': histogram)
       (int)    ci    : the color index (default 1)
       (int)   width  : linewidth (def. 0 = use previous)
       (int)    ls    : linestyle (def. 1 = solid line)
  """

  myTiming = Timing()
  
  global xAxis, yAxis, point, line

  if xAxis['log']:
    dataX = log(array(dataX))/log(10.)
  else:
    dataX = array(dataX)
  if yAxis['log']:
    dataY = log(array(dataY))/log(10.)
  else:
    dataY = array(dataY)

  pgsci(ci)                                # set colour index

  if   style == 'p':                       # points
    pgsch(point['size'])              # set character height
    pgpt(dataX,dataY,point['symbol']) # plot points
    
  elif style == 'l':                       # lines
    if ls:
      pgsls(ls)
    else:
      pgsls(line['linestyle'])          # set line style 
    if width:
      pgslw(width)                         # set line width
    else:
      pgslw(line['linewidth'])        # set line width
    pgline(dataX,dataY)                    # plot lines
    
  elif style == 'b':                       # bins
    pgsls(line['linestyle'])          # set line style 
    pgslw(line['linewidth'])          # set line width
    pgbin(dataX,dataY)

  __MessHand.debug("plotDataXY done in "+str(myTiming))
#----------------------------------------------------------------------------
def plotBox():
  """
  DES: Draw box and labels.
  """
  
  myTiming = Timing()

  global box, xAxis, yAxis, numPlot

  pgsch(box['charheight'])           # set character height 
  pgsci(box['color'])                # set colour index 
  pgsls(box['linestyle'])            # set line style
  pgslw(box['linewidth'])            # set line width 
  
  if (xAxis['log']):
    if (yAxis['log']):
      pgbox('BCTSNL',0.0,0,'BCTSNVL',0.0,0)
    else:
      pgbox('BCTSNL',0.0,0,'BCTSNV',0.0,0)
  else:
    if (yAxis['log']):
      pgbox('BCTSN',0.0,0,'BCTSNVL',0.0,0)
    else:
      pgbox('BCTSN',0.0,0,'BCTSNV',0.0,0)
      
  __MessHand.debug("plotBox done in "+str(myTiming))

#--------------------------------------------------------------------------------
def setLimits(dataX,dataY,limitsX=[],limitsY=[]):
  """
  DES: compute and/or set the limits for the graph
  INP: (arrays)        dataX/Y : the array to be plotted
  OPT: (2elts array) limitsX/Y : limits to use in X/Y for the plot
  """
  
  myTiming = Timing()

  global xAxis, yAxis

  if limitsX != []:
    # limitsX is set so use it
    xAxis['limits'] = array(limitsX)
  elif limitsY != []:
    #limitsX is not set but limitsY is set so compute the extrema
    #for dataX in the range defined by limitsY
    
    X_insideY = take(dataX,nonzero(greater_equal(dataY,min(limitsY)) & \
                                    greater_equal(max(limitsY),dataY)))

    xAxis['limits'] = minmax(X_insideY)

  else:
    # nothing is set, compute the extrema
    xAxis['limits'] = minmax(dataX)
    
  if xAxis['log']:
    xAxis['limits'] = log(xAxis['limits'])/log(10.)

  # Same for Y
  if limitsY != []:
    yAxis['limits'] = array(limitsY)
  elif limitsX != []:
    Y_insideX = take(dataY,nonzero(greater_equal(dataX,min(limitsX)) & \
                                    greater(max(limitsX),dataX)))
    
    yAxis['limits'] = minmax(Y_insideX)
  else:
    yAxis['limits'] = minmax(dataY)

  if yAxis['log']:
    yAxis['limits'] = log(yAxis['limits'])/log(10.)

  # If a limit is infinite, set it to 'nan'
  if str(xAxis['limits'][0]) == '-inf' or str(xAxis['limits'][0]) == 'inf':
    xAxis['limits'][0] = float('nan')
  if str(xAxis['limits'][1]) == '-inf' or str(xAxis['limits'][1]) == 'inf':
    xAxis['limits'][1] = float('nan')
  if str(yAxis['limits'][0]) == '-inf' or str(yAxis['limits'][0]) == 'inf':
    yAxis['limits'][0] = float('nan')
  if str(yAxis['limits'][1]) == '-inf' or str(yAxis['limits'][1]) == 'inf':
    yAxis['limits'][1] = float('nan')
  
  # Sometimes... an error occurs... and you ask for blank or
  # constant plots so...
  
  if xAxis['limits'][0] == xAxis['limits'][1]:
    xAxis['limits'][0] = xAxis['limits'][0] - 0.5 # this is completly arbitrary... but
    xAxis['limits'][1] = xAxis['limits'][1] + 0.5 # we are plotting constant so...

  if yAxis['limits'][0] == yAxis['limits'][1]:
    yAxis['limits'][0] = yAxis['limits'][0] - 0.5 
    yAxis['limits'][1] = yAxis['limits'][1] + 0.5

  # add a margin to clearify the plot
  dist=max(xAxis['limits'])-min(xAxis['limits'])
  # crude way to determine the order of the limits
  if str(dist) != str(float('nan')):
    if xAxis['limits'][0] < xAxis['limits'][1]:
      xAxis['limits'] = xAxis['limits'] + array([-1.,1.])*addFracLim*dist
    else:
      xAxis['limits'] = xAxis['limits'] - array([-1.,1.])*addFracLim*dist
    
  dist=max(yAxis['limits'])-min(yAxis['limits'])
  if str(dist) != str(float('nan')):
    if yAxis['limits'][0] < yAxis['limits'][1]:
      yAxis['limits'] = yAxis['limits'] + array([-1.,1.])*addFracLim*dist
    else:
      yAxis['limits'] = yAxis['limits'] - array([-1.,1.])*addFracLim*dist

  __MessHand.debug("setLimits done in "+str(myTiming))


#--------------------------------------------------------------------------------
def setMapTransformation( map_array, sizeX = [0.,1.], sizeY = [0.,1.], WCS = []):
  """
  DES: compute transformation matrix
  """
  global transformationMatrix

  if len(map_array.shape) == 2:
    nX = map_array.shape[0]
    nY = map_array.shape[1]
  elif len(map_array.shape) == 3:
    nX = map_array[0].shape[0]
    nY = map_array[0].shape[1]

  if WCS != []:
    if WCS['CUNIT1'] == 'arcsec':
      cosY = cos(WCS['CRVAL2']*pi/180./3600.)
    else:
      cosY = cos(WCS['CRVAL2']*pi/180.)
    cdeltCos = WCS['CDELT1'] / cosY
    transformationMatrix[0] = WCS['CRVAL1']-(cdeltCos)*(WCS['CRPIX1']+1)
    transformationMatrix[1] = cdeltCos
    transformationMatrix[2] = 0
    transformationMatrix[3] = WCS['CRVAL2']-WCS['CDELT2']*(WCS['CRPIX2']+1)
    transformationMatrix[4] = 0
    transformationMatrix[5] = WCS['CDELT2']
    
  else:
  
    # Making the assumption that the sizeX/sizeY array are sorted..
    # It could be useful to plot things in reverse order...

    x_pixel_size = (sizeX[1]-sizeX[0])/(nX-1)
    y_pixel_size = (sizeY[1]-sizeY[0])/(nY-1)
  
    transformationMatrix[0] = sizeX[0]-x_pixel_size # to be centered
    transformationMatrix[1] = x_pixel_size
    transformationMatrix[3] = sizeY[0]-y_pixel_size
    transformationMatrix[5] = y_pixel_size
      

#--------------------------------------------------------------------------------
def setMapLimits( map_array, limitsX=[], limitsY=[], limitsZ=[]):

  """
  DES: compute and/or set the limits for the map
  INP: (2D array)    map_array : the map array
  OPT: (2elts array) limitsX/Y : limits to use in X/Y for the plot
       (2elts array)  limitsZ  : the plotted color range in unit of the map_array       
  """

  myTiming = Timing()
  
  global xAxis, yAxis, zAxis, transformationMatrix
  
  # Compute the map limits in X and Y

  if len(map_array.shape) == 2:
    nX = map_array.shape[0]
    nY = map_array.shape[1]
  elif len(map_array.shape) == 3:
    nX = map_array[0].shape[0]
    nY = map_array[0].shape[1]

  # if the limits on X/Y are given use them, or use all the map
  if limitsX != []:
    xAxis['limits'] = limitsX
  else:
    xAxis['limits'] = transformationMatrix[0]+(transformationMatrix[1])*(array([0,nX])+0.5)
    
  if limitsY != []:
    yAxis['limits'] = limitsY
  else:
    yAxis['limits'] = transformationMatrix[3]+(transformationMatrix[5])*(array([0,nY])+0.5)

  # if limits on Z are given use them
  # else compute the minmax of the map, restricted to the limitsX/Y  
  if limitsZ != []:
    zAxis['limits'] = limitsZ
  else:
    # compute pixel positions corresponding to limitsX/Y
    tm = transformationMatrix
    x1 = int((xAxis['limits'][0] - tm[0])/tm[1] + 0.5)
    x2 = int((xAxis['limits'][1] - tm[0])/tm[1] + 0.5)
    y1 = int((yAxis['limits'][0] - tm[3])/tm[5] + 0.5)
    y2 = int((yAxis['limits'][1] - tm[3])/tm[5] + 0.5)
    lmap_array = ravel(array(map_array[x1:x2,y1:y2]))
    zAxis['limits'] = minmax(lmap_array)
  
  __MessHand.debug("setMapLimits done in "+str(myTiming))

#--------------------------------------------------------------------------------
def setImaCol( style = 'g2r', contrast=1, brightness=0.5, transferFunction=0,nan=0):
  """
  NAM: setImaCol (method)
  DES: Set image colours.
  OPT: (string) style    : the style (default 'g2r' : green to red)
                           also defined 'r2g' red to green, 'b2r' blue to red,
                           'r2b' red to blue, 'blue' or any lut file define in
                           the BogliConfig.lutDir variable
       (float)  contrast   : the contrast of the plot (default 1)
       (float)  brightness : the brightness of the plot (default 0.5)
       (int) transferFunction : the transfer funtion (linear/log/sqrt)
  """

  (iCiLo, iCiHi) = pgqcir() # color index used for the images
  pgscir(iCiLo,iCiHi)       # set colour index range 

  ncolor = 255
  ramp_intensity = arrayrange(ncolor*1.0)/ncolor # set the color ramp

  
  if style is 'g2r':
    red   = 1.0-ramp_intensity
    green = ramp_intensity
    blue  = ramp_intensity*0.0
  elif style is 'r2g':
    red   = ramp_intensity
    green = 1.0-ramp_intensity
    blue  = ramp_intensity*0.0
  elif style is 'b2r':
    red   = 1.0-ramp_intensity
    green = ramp_intensity*0.0
    blue  = ramp_intensity
  elif style is 'r2b':
    red   = ramp_intensity
    green = ramp_intensity*0.0
    blue  = 1.0-ramp_intensity
  elif style is 'blue':
    red   = ramp_intensity*0.0
    green = ramp_intensity*0.0
    blue  = ramp_intensity
  elif style in ['grey','gray']:
    red = 1.0-ramp_intensity
    green = 1.0-ramp_intensity
    blue = 1.0-ramp_intensity
  elif style is 'inverse':
    red = ramp_intensity
    green = ramp_intensity
    blue = ramp_intensity

  else:
    red,green,blue = readLut(style)
    ncolor = len(red)
    ramp_intensity = arrayrange(ncolor*1.0)/ncolor # set the color ramp
    
  pgctab(ramp_intensity,red,green,blue,ncolor,contrast,brightness)

  if nan:
    # Blank the lowest value (background) I do not understand why...
    pgscr(iCiHi,0,0,0)

  pgsitf(transferFunction) # set image transfer function
                           # ITF = 0 : linear
                           # ITF = 1 : logarithmic
                           # ITF = 2 : square-root

#--------------------------------------------------------------------------------
def readLut(lutFile):
  """
  NAM: readLut (method)
  DES: read a LUT file
  INP: (string) lutFile : the name of the input lut file
  """

  try:
    f = file(lutDir+lutFile+'.lut')
  except IOError:
    __MessHand.setMess(1," E: could not open file %s"%(lutFile))
    return
  
  # read and process LUT file
  lines = f.readlines()
  f.close()
  red, green, blue = [], [], []
  
  for i in range(len(lines)):
    tmp = string.split(lines[i])
    if tmp[0] != '!':			# skip comments
      red.append(string.atof(tmp[0]))
      green.append(string.atof(tmp[1]))
      blue.append(string.atof(tmp[2]))

  # Lut are defined in the inverse order so
  red.reverse()
  green.reverse()
  blue.reverse()

  # we need array to output
  return array(red), array(green), array(blue)

#--------------------------------------------------------------------------------
def plot(dataX, dataY = [], limitsX=[], limitsY=[], labelX = 'x', labelY = 'y',
         caption='', style='p', ci=1, width=0, ls=1, overplot=0, aspect=0,
         logX=0, logY=0, nodata=0, noerase=0):
  """
  DES: do a plot
  INP: (array)  dataX = values to plot along X
       (array)  dataY = values to plot along Y (optional - default:
                                     plot dataX vs. running number)
  OPT: (2elts array) limitsX/Y = limits to use in X/Y for the plot
       (string)      labelX  = x label (default 'x')
       (string)      labelY  = y label (default 'y')
       (string)      caption = the caption of the plot (default ' ')
       (char)        style   = the style used for the plot ('l': line,
                               'p': point (default), 'b': histogram)
       (int)            ci   = color index (default 1)
       (int)            ls   = line style (def 1) FB240307
       (int)          width  = linewidth (defaut 0 = use previous)
       (logical)      aspect = keep the aspect ratio in 'physical' unit
       (logical)    overplot = are we overplotting ? (default no)
       (logical)      logX/Y = do we use log scale ? (default no)
       (logical)      nodata = do not plot the data
       (logical)     noerase = do not clear the window
  """

  myTiming = Timing()

  global xAxis, yAxis, wasWedge, wasAspect, numPlot

  if dataY != [] :
    dataX = array(dataX)
    dataY = array(dataY)
    # Make sure that X and Y have the same length
    if len(dataX) > len(dataY):
      dataX = dataX[:len(dataY)]
    elif len(dataY) > len(dataX):
      dataY = dataY[:len(dataX)]
      
  else:
    dataY = dataX
    dataX = array(range(len(dataY)))
      
  # Remove all NaN from X and Y arrays
  maskX = masknan(dataX)
  maskY = masknan(dataY)
  mNan = bitwise_or(maskX,maskY)
  bad = nonzero(not_equal(mNan,0))
  if len(bad):
    # at least one NaN somewhere
    __MessHand.warning("%i data points NaN were found"%(len(bad)))
    good = nonzero(equal(mNan,0))
    if len(good):
      dataX = take(dataX,good)
      dataY = take(dataY,good)
    else:
      dataX = array([0],'f')
      dataY = array([0],'f')      

  if not overplot:
    # Check for an open device first
    if not DeviceHandler.CurrentDev:
      __MessHand.warning('no device open')
      DeviceHandler.openDev()
    else:
      # to overcome the closing-without-closing feature of apexCalibrator
      DeviceHandler.__checkReopen()

    if not noerase:
      clear()

    if logX:
      xAxis['log'] = 1
    if logY:
      yAxis['log'] = 1
    setLimits(dataX,dataY,limitsX=limitsX,limitsY=limitsY)
    setLabels(labelX=labelX,labelY=labelY,caption=caption)

    wasAspect = aspect
    numPlot = 1
    setViewPoint(aspect=aspect)

  if overplot:
    if numPlot != 1:
      __MessHand.error('cannot overplot single Plot on multiPlot')
      return

  if not nodata:
    plotDataXY(dataX,dataY,style,ci,width,ls)      
    
  if not overplot:
    if box['draw_it']:
      plotBox()
    labelling()

  __MessHand.debug("plot done in "+str(myTiming))

#--------------------------------------------------------------------------------
# Utilitie for one-D plots
#
def getPixel(order=0):
  """
  DES: allow user to get pixel values using mouse
  INP: (int)  order = for polynomial interpolation
  """
  
  __MessHand.info("Click left to get one pixel, right (or 'e' or 'x') to exit")
  again = 1
  x,y = 0,0
  nbData = 0
  posX = []  # store positions to perform polynomial interpolation
  posY = []
  __MessHand.info(str("(%0.10s, %0.10s)\n"%(str("%10s"%xAxis['text']),\
                                            str("%10s"%yAxis['text']))))
  while(again):
    x,y,char = pgband(7,0,x,y)
    if (char in ['X','x','e','E']):
      again = 0
    else:
      if xAxis['log']:
        xValue = 10.**x
      else:
        xValue = x
      if yAxis['log']:
        yValue = 10.**y
      else:
        yValue = y
      # display position
      __MessHand.info(str("(%10.3f, %10.3f)"%(xValue,yValue)))
      if order:
        #if (nbData > order):  # remove first point
        #  old=posX.pop(0)
        #  old=posY.pop(0)
        posX.append(x)        # append new one
        posY.append(y)
        if nbData:            # if at least two points
          coeff = solvePoly(min([order,nbData]),posX,posY)
          out = str(" Interpolate %1i last points: Y = " % (order+1))
          for n in range(order+1):
            if (coeff[n] > 0):
              out += str("%+5.2f X^%1i " % (coeff[n],order-n))
            else:
              out += str("%6.2f X^%1i " % (coeff[n],order-n))
          __MessHand.info(out)
          
          #### TODO: overplot the polynom
        nbData += 1


def draw( map_array, sizeX=[], sizeY=[], WCS = [], \
         limitsX= [], limitsY= [], limitsZ=[], \
         nan = 1, \
         labelX = 'x', labelY = 'y', caption='', \
         style='g2r', contrast=1.0, brightness=0.5, \
         wedge=0, overplot=0, aspect=0, \
         doContour=0, levels=[], labelContour=0, noerase=0  ):
  
  """
  DES: do a image drawing
  INP: (map_array)  map to display
  OPT: (2elts arrays)   sizeX/Y = the 'physical' size of the array (default pixel numbers)
                                  defined by the center of the two extreme pixels !
       (2elts arrays) limitsX/Y = limits to use in X/Y for the plot
       (logical)           nan  = set if NaN are present in the array
       (string)         labelX  = x label (default 'x')
       (string)         labelY  = y label (default 'y')
       (string)         caption = the caption of the plot (default ' ')
       (char)           style   = the color used for the plot (default 'g2r'
                                  see BogliPlot.Plot.setImaCol())
       (logical)        wedge   = shall we draw a wedge ? (default no)
       (logical)         aspect = keep the aspect ratio in 'physical' unit
       (logical)       overplot = are we overplotting ? (default no)
       (logical)     doContour  = draw contour instead of map (default no)
       (array)         levels   = the levels for the contours (default nContour
                                  withing plotLimitsZ)
       (logical)   labelContour = label the contours (default no)
       (logical)        noerase = do not clear the window
  """

  myTiming=Timing()
  
  global wasWedge, wasAspect, transformationMatrix, zAxis, contour, numPlot

  # Check for an open device first
  if not DeviceHandler.CurrentDev:
    __MessHand.warning('no device open')
    DeviceHandler.openDev()
  else:
    # to overcome the closing-without-closing feature of apexCalibrator
    DeviceHandler.__checkReopen()

  lmap_array = array(copy.copy(map_array))

  # points in the map that were not observed are at NaN,
  # retrieved them... ugly
  if nan:
    lmap_array = replacenan(lmap_array)
    
  # If the size of the map are not gived use a integer range starting
  # at 0

  if sizeX == []:
    sizeX = array([0.5,lmap_array.shape[0]-0.5],Float)
  else:
    # Just to be sure that sizeX is an array of Float
    sizeX = array(sizeX,Float)

  # Same for sizeY
  if sizeY == []:
    sizeY = array([0.5,lmap_array.shape[1]-0.5],Float)
  else:
    sizeY = array(sizeY,Float)

  # Compute the transformation matrix for this image even in case of overplot
  
  setMapTransformation(lmap_array, sizeX=sizeX, sizeY=sizeY, WCS=WCS)

  if not overplot:
    if not noerase:
      clear()

    setMapLimits(lmap_array, limitsX=limitsX, limitsY=limitsY, limitsZ=limitsZ)
    setLabels(labelX=labelX,labelY=labelY,caption=caption)
    setImaCol(style=style, contrast=contrast, brightness=brightness,nan=nan)

    # Keep in mind if there was a wedge
    wasWedge = wedge
    wasAspect = aspect
    numPlot = 1

    setViewPoint(aspect=aspect,wedge=wedge)

    
  if not doContour:
    pgimag(transpose(lmap_array),lmap_array.shape[0],lmap_array.shape[1],\
           0, lmap_array.shape[0]-1,0,lmap_array.shape[1]-1,\
           zAxis['limits'][0],zAxis['limits'][1],\
           transformationMatrix)
  else:
    # set the line style first
    pgsch(contour['charheight']) # set character height
    pgsci(contour['color'])       # set colour index
    pgslw(contour['linewidth'])   # set line width


    # If no levels have been set compute them
    if levels == []:
      nContour = 10
      levels = arrayrange(nContour*1.0)/(nContour+1)*\
               (zAxis['limits'][1]-zAxis['limits'][0])+zAxis['limits'][0]

    # Draw the contour
    pgcont(transpose(lmap_array),lmap_array.shape[0],lmap_array.shape[1],\
           0, lmap_array.shape[0]-1,0,lmap_array.shape[1]-1,\
           array(levels),len(levels),\
           transformationMatrix)

    if labelContour == 1: 
      # Label the contour
      for level in array(levels):
        pgconl(transpose(lmap_array),lmap_array.shape[0],lmap_array.shape[1],\
               0, lmap_array.shape[0]-1,0,lmap_array.shape[1]-1,\
               level,transformationMatrix, str("%7.3f"%level),20,10)
        
  if not overplot:
    plotBox()
    labelling(wedge=wedge)

  __MessHand.debug("draw done in "+str(myTiming))
    
#--------------------------------------------------------------------------------
def panels(nX,nY):
  """
  DES: divide the plot window in nX x nY sub-panels
  INP: (int) nX,nY = number of panels in X and Y
  """
  pgsubp(nX,nY)
  
#--------------------------------------------------------------------------------
def nextpage():
  """
  DES: start a new page or go to next panel
  """
  pgpage()

#--------------------------------------------------------------------------------
def getpix(x,y):
  """
  DES: get pixel value using the mouse pointer
  """
  x,y,char = pgband(7,0,x,y)
  return x,y,char
#--------------------------------------------------------------------------------
# Polygon methods
#--------------------------------------------------------------------------------
def defPolygon():
  """
  DES: define a polygon
  """

  # Check for an open device first
  if not DeviceHandler.CurrentDev:
    __MessHand.warning('no device open')
    DeviceHandler.openDev()

  poly=[]
  cursPos=pgband(7,0,0,0)
  startPos=cursPos
  poly.append((cursPos[0],cursPos[1]))
  __MessHand.info(' x = '+`cursPos[0]`+' y = '+`cursPos[1]`)
  if cursPos[2]=='X':
    __MessHand.warning('not enough vertices')
  else:
    pgmove(startPos[0],startPos[1])
    cursPos=pgband(1,0,cursPos[0],cursPos[1])
    poly.append((cursPos[0],cursPos[1]))
    __MessHand.info(' x = '+`cursPos[0]`+' y = '+`cursPos[1]`)
    if cursPos[2]=='X':
      __MessHand.warning('not enough vertices')
    else:
      pgdraw(cursPos[0],cursPos[1])
      while 1:
        cursPos=pgband(1,0,cursPos[0],cursPos[1])
        if cursPos[2]=='X':
          pgdraw(startPos[0],startPos[1])
          break
        poly.append((cursPos[0],cursPos[1]))
        __MessHand.info(' x = '+`cursPos[0]`+' y = '+`cursPos[1]`)
        pgdraw(cursPos[0],cursPos[1])
        pgmove(cursPos[0],cursPos[1])
  poly=array(poly)
  return poly
#--------------------------------------------------------------------------------
def savePolygon(poly=zeros((1,2)),outFile='polygon1.pol'):
  """
  DES: save polygon to a file
  """
  try: 
    f=file(outFile,'w')      
  except IOError:
    __MessHand.error('could not open file %s in write mode'%(outFile)) 
    return 
  try: 
    cPickle.dump(poly,f)
    f.close()
  except: 
    __MessHand.error('cannot save polygon')
  return

#--------------------------------------------------------------------------------
def loadPolygon(inFile='polygon1.pol'):
  """
  DES: load polygon from a file
  """
  try: 
    f=file(inFile,'r')      
  except IOError: 
    __MessHand.error('could not open file %s in read mode'%(inFile))
    return
  try:
    poly=cPickle.load(f)
    f.close()
  except:
    __MessHand.error('cannot load polygon')
    return
  return poly
#--------------------------------------------------------------------------------
def plotPolygon(poly=zeros((1,2))):
  """
  DES: plot a polygon
  """
  pgline(poly[:,0],poly[:,1])
  pgdraw(poly[0,0],poly[0,1])

## @}
