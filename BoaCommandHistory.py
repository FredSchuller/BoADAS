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
NAM: BoaCommandHistory.py (module)  
DES: provides methods to access the command history
"""
__version__=  '$Revision: 1730 $'
__date__=     '$Date: 2006-09-26 22:16:08 +0200 (Tue, 26 Sep 2006) $'

#----------------------------------------------------------------------------------
#----- Import ---------------------------------------------------------------------
#----------------------------------------------------------------------------------
import readline

from boa import BoaMessageHandler

#----------------------------------------------------------------------------------
#----- _CommandHistory Class ------------------------------------------------------
#----------------------------------------------------------------------------------

class _CommandHistory:
    # NAM: _CommandHistory (class)
    # DES: Encapsulate calls to the readline module to access tho command history

    def __init__(self):
        self._messHand=BoaMessageHandler.MessHand(self.__module__)

        self._hasGetHistoryItem = 0

        try:
            readline.read_history_file()
        except IOError:
            self._messHand.warning("Not able to retrieve last session commands")
            return
        
        if not hasattr(readline,"get_history_item"):
            # For Python earlier than Version 2.3
            self._messHand.warning("Command history not avaiable for this Python version")
            return
        else:
            self._hasGetHistoryItem = 1

            self._tags = []

            try:
                curHistLen = readline.get_current_history_length()
            except:
                curHistLen = 0

            if curHistLen:
                self._tags.append(curHistLen+1)
            else:
                self._tags.append(1)


    def __del__(self):
        try:
            readline.write_history_file()
        except IOError:
            self._messHand.warning("Not able to save session commands")

    def addTag(self):
        if self._hasGetHistoryItem:
            try:
                self._tags.append(readline.get_current_history_length())
            except:
                self._messHand.warning("Not able to add tag to command history")

    def getCommands(self, itag=0):
        ret = []
        if self._hasGetHistoryItem:
            try:
                iStop = readline.get_current_history_length()
                for i in xrange(self._tags[itag], iStop+1):
                    ret.append(readline.get_history_item(i))
            except:
                self._messHand.warning("Not able to get command history")
        return ret


#----------------------------------------------------------------------------------
#----- The only _CommandHistory object --------------------------------------------
#----------------------------------------------------------------------------------

_theHistory = None

#----------------------------------------------------------------------------------
#----- The interface of the module ------------------------------------------------
#----------------------------------------------------------------------------------

def initHistory():
    """
    DES: Initialize the command history from the .history file
    """
    global _theHistory
    if not _theHistory:
        _theHistory = _CommandHistory()

def tagHistory():
    """
    DES: Adds a tag to the command history.
         Tags are integer numbers starting at 0. 
    """
    global _theHistory
    if _theHistory:
        _theHistory.addTag()

def getHistory(itag=0):
    """
    DES: Returns the command history starting at tag itag
    INP: (int) itag: Tag from whereon commands are returned.
                     For itag=0, all command since the call of
                     initHistory are returned; for itag=-1, the
                     command since the last call of tagHistory are
                     returned.
    OUT: (list str): Commands in the history
    """
    global _theHistory
    if _theHistory:
        return _theHistory.getCommands(itag)
    else:
        return None

def writeHistoryToFile(filename, itag=0):
    """
    DES: Write command history to the specified file
    INP: (str filename): The full filename
    """
    global _theHistory
    if _theHistory:
        of = open(filename, 'w')
        commands = getHistory(itag)
        for command in commands:
            of.write("%s\n" % command)
        of.close()
    else:
        self._messHand.warning("Not able to write history to file")
        



