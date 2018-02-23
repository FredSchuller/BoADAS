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
NAM: DeviceHandler.py (module)
"""
__version__=  '$Revision: 2381 $'
__date__='$Date: 2007-08-29 10:06:23 +0200 (Wed, 29 Aug 2007) $'
#----------------------------------------------------------------------------------
#----- Import ---------------------------------------------------------------------
#----------------------------------------------------------------------------------
import ppgplot
import string
from boa import BoaMessageHandler
#----------------------------------------------------------------------------------
#----- BoGLi Device Handler -------------------------------------------------------
#----------------------------------------------------------------------------------

# attributes:
CurrentDev=0    
DevList=[]
__MessHand = BoaMessageHandler.MessHand('DeviceHandler')
        
#--------------------------------------------------------------------------------
#----- public methods -----------------------------------------------------------
#--------------------------------------------------------------------------------
def openDev(type='/XWINDOW'):
  """
  DES: open a device, return the device id
  INP: (string) type = pgplot device type 
  """ 

  global CurrentDev, DevList
  
  if isinstance(type,str):
    CurrentDev=ppgplot.pgopen(type)
    ppgplot.pgask(0)
    DevList.append(CurrentDev)
  else:
    __MessHand.setMess(0, "openDev: string argument required ")
  __queryDev()
#--------------------------------------------------------------------------------
def selectDev(devID = ''):
  """
  DES: select an open device
  INP: (int) device ID
  """

  global CurrentDev, DevList

  if devID in DevList:
    ppgplot.pgslct(devID)             # select the deviceID
    CurrentDev=ppgplot.pgqid()        # store the actual pgplot deviceID
  else: 
    __MessHand.error("this is not an open device: "+`devID`)
  __queryDev()


#--------------------------------------------------------------------------------
def closeDev(devID='current'):
  """
  DES: close selected device
  INP: (int) device ID, 'all','current' (default)
  """

  global CurrentDev, DevList

  if devID in ['current','curren','curre','curr','cur']:
    
    __MessHand.longinfo("closing the current device")
    ppgplot.pgclos();
    DevList.remove(CurrentDev)
    if not DevList==[]:
      ppgplot.pgslct(DevList[0])
      CurrentDev=ppgplot.pgqid()        # store the actual deviceID
    else:
      CurrentDev = 0

  elif devID in ['all']:
    
    __MessHand.longinfo("closing all devices")
    for device in DevList:
      ppgplot.pgslct(device)
      ppgplot.pgclos();

    DevList=[]
    CurrentDev=0
    
  elif devID in DevList:
    
    __MessHand.longinfo("closing device "+`devID`)
    ppgplot.pgslct(devID)
    ppgplot.pgclos()
    DevList.remove(devID)
    if devID == CurrentDev and DevList != []:
      ppgplot.pgslct(DevList[0])
      CurrentDev=ppgplot.pgqid()        # store the actual deviceID
    else:
      ppgplot.pgslct(CurrentDev)
  else: 
    __MessHand.error("this is not an open device: "+`devID`)
  __queryDev()
  
#--------------------------------------------------------------------------------
def resizeDev():
  """
  DES: resize plot area after resizing window with mouse
  ABB: resize
  """
  __MessHand.info("resizing plot page")
  ppgplot.pgpage()


#--------------------------------------------------------------------------------
def resizeDev():
  """
  Des: resize plot area after resizing window with mouse
  ABB: resize
  """
  __MessHand.info("resizing plot page")
  ppgplot.pgpage()

  
#--------------------------------------------------------------------------------
#----- private methods ----------------------------------------------------------
#--------------------------------------------------------------------------------
def __queryDev():
  """
  DES: query open and current devices
  INP: (int) a = message weight
  """

  global CurrentDev, DevList
  
  DevList.sort()
  CurrentDev=ppgplot.pgqid()        # store the actual deviceID

  if CurrentDev:
    __MessHand.info("current device = "+`CurrentDev`)
    if not CurrentDev in DevList:
      # for some reason a device is open and not in the DevList
      DevList.append(CurrentDev)
  else:
    __MessHand.warning("current device = none")

  if len(DevList) > 0:
    __MessHand.longinfo("open devices   = "+`DevList`)
  else:
    __MessHand.warning("open devices   = none")

def __checkReopen():
  """
  DES: check if the current device is actually open and usable, or if it
       has to be reopenned (e.g. after a plot produced by apexCalibrator)
  """

  global CurrentDev, DevList

  ppgId = ppgplot.pgqid()
  if ppgId != CurrentDev:
    CurrentDev=ppgId
    ppgplot.pgask(0)
    DevList.append(CurrentDev)
