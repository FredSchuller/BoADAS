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
NAM: BoaMessageHandler.py (module)  
DES: contains the BoA message handler class
"""
__version__=  '$Revision: 2794 $'
__date__=     '$Date: 2015-03-02 15:30:03 +0100 (Mon, 02 Mar 2015) $'

#----------------------------------------------------------------------------------
#----- Import ---------------------------------------------------------------------
#----------------------------------------------------------------------------------
import os
import string
import sys
import time
from boa import BoaConfig

# Keeps all the messages in memory
_history = {} # for test

def printHistory():
  allkeys = _history.keys()
  allkeys.sort()
  
  for key in allkeys:
    print key, _history[key]

#----------------------------------------------------------------------------------
#----- BoA Message Handler --------------------------------------------------------
#----------------------------------------------------------------------------------

class MessHand:
  """
  NAM: MessHand (class)
  DES: An object of this class is responsible for the management of output 
       messages as well as the creation of message files. 
  """ 


  def __init__(self, logName='Unknown'):
    """
    DES: initialise an instance 
    """

    # parameter attributes:
 
    # default value for max. weight of messages to be printed
    self.maxWeight = BoaConfig.maxMessHandWeight
    
    # Private attributes:

    self.__logName = '' # Name of the class using the Message Handler
                        # 10 char long left hand stripped.
    self.__WeightList = {1: 'E:', 2: 'W:',3: 'I:', 4: 'L:', 5: 'D:'}  # list of allowed weights 
    self.__WeightDescription = ['errors', 'warnings', \
                                'short info', 'extended info', \
                                'debug']

    self.__messFileName="boa.mes"                                   # default name of message file
#    self.__prompt = chr(1) + "\033[1;32m" + chr(2) + 'boa>' + chr(1) + "\033[0m" + chr(2)
    self.__prompt = 'boa<'

    # methods to run at initialisation 
    self.setLogName(logName)


  def __del__(self):
    self.debug("closing message handler '"+self.__logName.strip()+"'")
    
  #--------------------------------------------------------------------------------
  #----- methods ------------------------------------------------------------------
  #--------------------------------------------------------------------------------
  def Welcome(self):
    """
    DES: print welcome message
    OUT: screen output
    """
    str = ""
    print
    print " Welcome to The BOlometer data Analysis project ! "
    print "   -------------------------------------------    "

    if os.environ.get('USER'):
      str += " User: "+os.environ.get('USER')+" "
    str += "("+ time.ctime()+ ")\n"
    #   str += "Running on " +sys.platform+" with Python "+sys.version
    print str

  #--------------------------------------------------------------------------------
  def setMaxWeight(self,weight='2'):
    """
    DES: Set the maximum weight of messages to be printed.
    INP: (int) weight = maximum weight
    
                  1: errors, queries
                  2: warnings
                  3: short info
                  4: extended info
                  5: debug
    """
    if weight in self.__WeightList: 
      self.maxWeight=weight
      self.info("max weight of messages = "+`self.maxWeight`)
    else:
      self.error("invalid max weight of messages: "+`weight`)
        
  #--------------------------------------------------------------------------------
  def setLogName(self,logName='Unknown'):
    self.__logName = logName.ljust(15)    # Name of the class using
    self.__logName = self.__logName[0:15] # the Message Handler 10 char  
                                          # long left hand stripped

  #--------------------------------------------------------------------------------
  def pause(self,message= ''):
    """
    DES: allow to make a pause in the program
    OPT: (sring) : a message to display
    """
    if not message:
      message='Please press enter to continue'

    self.ask(message)


  #--------------------------------------------------------------------------------
  def yesno(self,message= ''):
    """
    DES: ask the user a question with yes/no answer type
    INP: (sring)  : the question
    OUT: (l) :    : the answer
    """

    yes = ['y','Y','yes','Yes','YES']
    no  = ['n','N','no','No','NO']
    answer = ''
    while answer not in yes and answer not in no:
      answer = self.ask(message+' (y/n): ')

    if answer in yes :
      return 1
    else:
      return 0
    
                                          
  #--------------------------------------------------------------------------------
  def ask(self,message= ''):
    """
    DES: ask the user 
    INP: (sring)  : the question
    OUT: (string) : the answer
    """
    return raw_input(self.__prompt + ' ?: ' + message)
                                          
  #--------------------------------------------------------------------------------
  def error(self,message= ''):
    """
    DES: to print an error message
    INP: (sring) message
    """
    self.setMess(1,message)

  #--------------------------------------------------------------------------------
  def warning(self,message=''):
    """
    DES: to print an warning message
    INP: (sring) message
    """
    self.setMess(2,message)
    
  #--------------------------------------------------------------------------------
  def info(self,message=''):
    """
    DES: to print an info message
    INP: (sring) message
    """
    self.setMess(3,message)
    
  #--------------------------------------------------------------------------------
  def longinfo(self,message=''):
    """
    DES: to print an long info message
    INP: (sring) message
    """
    self.setMess(4,message)
    
  #--------------------------------------------------------------------------------
  def debug(self,message=''):
    """
    DES: to print an debug message
    INP: (sring) message
    """
    self.setMess(5,message)

  #--------------------------------------------------------------------------------
  def line(self):
    """
    DES: to print a line
    """
    message = '-'*(80-10)
    self.setMess(4,message)

  #--------------------------------------------------------------------------------
  def setMess(self,weight=1,message=' '):
    """
    DES: deposit messages for screen output and message files
    INP: (int) weight = weight of transferred message (see setMaxWeight)
         (string) message = message to be printed and added to message file 
    """
    if weight in self.__WeightList:
      message = message.strip() 
      short_prefix = " " +self.__WeightList[weight] + " "
      long_prefix  = time.ctime() + " " + self.__logName + short_prefix
      
      if weight<=self.maxWeight:                         # Print if asked ...

        # ... in file
        if (self.__messFileName != ""):
          messFile=open(self.__messFileName,'a')
          for line in message.split('\n'):
            messFile.write(long_prefix + line + "\n")
          messFile.close()
          del messFile
        else:
          
          print self.__prompt + " " + self.__WeightList[1] + " " + \
                "the log file for this message handler is closed or undefined"
          
        
        # ... on screen
        for line in message.split('\n'):
          print self.__prompt + short_prefix + line

        # and add this to history
        _history[long_prefix] = message

        del short_prefix, long_prefix
      
  #--------------------------------------------------------------------------------
  def initMessFile(self, filename="boa.mes"):
    """
    DES: set & initialise new message file
    OUT: screen output  
    """
    if self.__messFileName=="":
      self.__messFileName=BoaConfig.outDir+filename.strip()

    try:
      
      if os.path.isfile(self.__messFileName):
        # File name already exist so move it to a new location
        # First find a available file name as self.__messFileName+i
        i=1
        while(os.path.isfile(self.__messFileName+str(i))):
              i+=1
        newFileName = self.__messFileName+str(i)
        self.info('old log file renamed to '+newFileName)
        os.rename(self.__messFileName,newFileName)
        
      messFile=open(self.__messFileName,'w')
      messFile.write("Logfile for The BOlometer data Analysis project. ")
      messFile.write("created on "+ time.ctime() + "\n")
      messFile.close()
      del messFile
      
    except IOError:
      self.error("cannot open message file "+self.__messFileName)

  #--------------------------------------------------------------------------------
  def closeMessFile(self):
    """
    DES: set self.__existMessFile to 0 and file name to "" 
    """
    self.debug("closing message file "+self.__messFileName)
    self.__messFileName=""

#----------------------------------------------------------------------------------
#----- For compatibility with the CalibratorLog class -----------------------------
#----------------------------------------------------------------------------------
class Logger:
  """
  NAM: Logger (class)
  DES: for compatiliby with the CalibratorLog.Logger class
  """
  
  def __init__(self, logType='ACS'):
    """
    DES: Initiabise an instance
    """
    
    self.logger = printLogger()
    
class printLogger(MessHand):
  """
  NAM: printLogger (class)
  DES: for compatibility with the CalibratorLog.printLogger class
  """
  def __init__(self):

    MessHand.__init__(self)
    self.setLogName(logName='CalibratorLog')

  #--------------------------------------------------------------------------------
  # Compatibility with the CalibratorLog class
  def logInfo(self,message):
    self.info(message)

  def logError(self,message):
    self.error(message)

  def logWarning(self,message):
    self.warning(message)
  
  def logDebug(self,message):
    self.debug(message)

    
