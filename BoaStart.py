#!/usr/bin/env python
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
NAM: BoaStart.py (file)
DES: BoA starting procedure 
"""

__version__=  '$Revision: 2794 $'
__date__=     '$Date: 2015-03-02 15:30:03 +0100 (Mon, 02 Mar 2015) $'


# ---------------------------------------------------------------------
# ---- Import ---------------------------------------------------------
# ---------------------------------------------------------------------

import os, sys
import string
from Numeric import *          # to be able to define new arrays and so on
import cPickle

import readline, rlcompleter
readline.parse_and_bind("tab: complete")

from boa.Bogli import *
from boa import BoaConfig
from boa import BoaMessageHandler, BoaCommandHistory, BoaFlagHandler
from boa import BoaMapping ,BoaPointing, BoaFocus
from boa.Utilities import getTau,getCalCorr


# ---------------------------------------------------------------------
# ---- Print welcome message ------------------------------------------
# ---------------------------------------------------------------------

print
print "BoA - the Bolometer Array data Analysis Software"
print
print "Versions: BoA        :", os.getenv('BOA_VERSION')
print "          BoA Library:", os.getenv('BOA_LIB_VERSION')
print

# ---------------------------------------------------------------------
# ---- Initialisation of the Message File -----------------------------
# ---------------------------------------------------------------------
messages = BoaMessageHandler.MessHand('User')
messages.initMessFile()

# ---------------------------------------------------------------------
# ---- Initialisation of the Command History --------------------------
# ---------------------------------------------------------------------
BoaCommandHistory.initHistory()
    
# ---------------------------------------------------------------------
# ---- Initialisation of the plotting routine -------------------------
# ---------------------------------------------------------------------

plot  = Plot.plot
draw  = Plot.draw
mplot = MultiPlot.plot
mdraw = MultiPlot.draw        

# ---------------------------------------------------------------------
# ---- Change the prompt if online mode  -----------------------------
# if you dont like the color, choose a different one from the menue below
# ---------------------------------------------------------------------

if BoaConfig.online:
    messages.info("BoA - Online mode")
else:
    messages.info("BoA - Offline mode")
    
sys.ps1 =  'boa> '
sys.ps2 =  '.... '


# ---------------------------------------------------------------------
# ---- Main program ---------------------------------------------------
# ---------------------------------------------------------------------

"""
DES: You have to define an environment variable MBFITSXML,
     whose value is the path to the file MBFits.xml,
     including the file name itself 
     (ex: setenv MBFITSXML ./MBFits.xml )
"""

# TODO make an external read() routine which will return the
# appropriate object and take care of memory allocation

#data = BoaMapping.Map()
#data = BoaPointing.Point()
data = BoaFocus.Focus()

try:
    execfile(os.getenv('BOA_HOME_BOA')+'/BoaShortcut.py')
except:
    messages.warning("BoaShortcut not found, or error in execution")


# ---------------------------------------------------------------------
# ---- Try to execute the interactive script --------------------------
# ---- The variable PYTHONSTARTUP should point to interactive.py ------
# ---------------------------------------------------------------------
try:
    if os.environ.get('PYTHONSTARTUP') \
           and os.path.isfile(os.environ['PYTHONSTARTUP']):
        execfile(os.environ['PYTHONSTARTUP'])
        
except:
    messages.warning("interactive.py could not be loaded")

# ---------------------------------------------------------------------
# Function definitions
# ---------------------------------------------------------------------
import gc
def test_memory():
    attrDic = vars(data)
    attrName = attrDic.keys()
    nbAttr = len(attrName)
    attrName.sort()

    for badattribute in ['Bogli','MessHand','addLatWT','addLonWT',\
                         'flagValue','fwhm','numDataPoints','phase','scanNum']:
        if badattribute in attrName:
            attrName.remove(badattribute)
    totsize = 0
    for a in attrName:
        totsize += len(attrDic[a])
        print a, " : ", len(attrDic[a])
        attrDic[a] = [0]
            
    print "Total :",totsize

def newRestoreData(fileName='BoaData.sav'):
    """ 
    DES: restore a DataEntity object previously saved in a file, and
    set it as the currData attribute of BoaB
    INP: (string) fileName: name of the input file
    optional - default value = 'BoaData.sav'
    """
    tmp = restoreFile(fileName)
    
    if hasattr(tmp,'DataFlags'):
        tmp.FlagHandler = BoaFlagHandler.createFlagHandler(tmp.DataFlags.astype(Int8))
        tmp.DataFlags = None
    
        tmp.ScanParam.FlagHandler = BoaFlagHandler.createFlagHandler(tmp.ScanParam.Flags.astype(Int32))
        tmp.ScanParam.Flags = None
    
        tmp.BolometerArray.FlagHandler = BoaFlagHandler.createFlagHandler(tmp.BolometerArray.Flags.astype(Int32))
        tmp.BolometerArray.Flags = None

    return tmp

def restoreFile(fileName):
    """
    DES: read an object that was stored using cPickle (e.g. using dumpData or dumpMap)
    """
    try:
        f = file(fileName)
    except IOError:
        messages.error(" could not open file %s"%(fileName))
        return 0
    tmp = cPickle.load(f)
    f.close()
    return tmp

# ---------------------------------------------------------------------
# Finally, execute BoaOnOff.py which contains functions to process
# ON-OFF scans observed with LABOCA
execfile(os.getenv('BOA_HOME_BOA')+'/BoaOnOff.py')
