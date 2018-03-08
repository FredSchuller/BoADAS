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
NAM: BoaFlagHandler.py (file)
DES: Contains classes for bitwise manipulation of flags.
     The module serves as an object-oriented interface for the Fortran
     module fFlagHandler.
"""                    
# ---------------------------------------------------------------------
# ---- Import ---------------------------------------------------------
# ---------------------------------------------------------------------

import copy
import Numeric
from fortran import fFlagHandler, fUtilities


# ---------------------------------------------------------------------
# ---- Function createFlagHandler B-------------------------------------
# ---------------------------------------------------------------------

def createFlagHandler(aFlags):
    """
    NAM: createFlagHandler (Function)
    DES: Creates an object of the correct subtype of FlagHandler to
         handle the flag array aFlags.
         Presently, 1-dim flag arrays of type Int32 and 2-dim flag arrays
         of type Int8 can be handled.
    INP: aFlags (array) : Flag array
    OUT: (object)       : FlagHandler object
    """
    typecode = aFlags.typecode()
    nDim = len(aFlags.shape)

    if nDim==1 and typecode==Numeric.Int32:
        return FlagHandler1d32b(aFlags)
    elif nDim==2 and typecode==Numeric.Int8:
        return FlagHandler2d8b(aFlags)
    else:
        raise "Invalid Flag array: dim=%d, typecode=%s" % (nDim, typecode)


# ---------------------------------------------------------------------
# ---- Class FlagHandler ----------------------------------------------
# ---------------------------------------------------------------------

class FlagHandler:
    """
    NAM: FlagHandler (Class)
    DES: Provides methods to manipulate and query flag arrays bitwise.
         Here, a flag array is a Numeric.array, each element of which
         represents n independent flags, with n the number of bits of
         each array element. The n independent flags are enumerated with
         the flag value, ranging from 1 (the rightmost bit) to n (the
         leftmost bit).

         Throughout the class, the nomenclature is as follows:

         - self._aFlags (Numeric.array):
             The flag array to be manipulated or queried.

             Permitted shapes and datatypes of self._aFlags depend on the
             exact subtype of FlagHandler.

         - iFlags (integer, list of integers, or empty list):
             The flag values used for manipulating or queriing the flag
             array.

             Valid flag values are from 1 to n, with n the number of bits
             of each element of aFlags.

             If iFlags is the empty list, the list [1,2,...,n] is assumed.
    """

    # -----------------------------------------------------------------

    def __init__(self, aFlags):
        self._aFlags = aFlags

        self._validFlagValues = range(1, self._aFlags.itemsize()*8+1)

    # -----------------------------------------------------------------

    def getFlags(self):
        """
        NAM: FlagHandler.getFlags (Method)
        DES: Returns the flag array self._aFlags
        """
        return self._aFlags

    # -----------------------------------------------------------------

    def getValidFlagValues(self):
        """
        NAM: FlagHandler.getValidFlagValues (Method)
        DES: Returns a list of all valid flag values for
             flag array self._aFlags
        OUT: (int list) : All valid flag values for aFlags
        """
        return copy.copy(self._validFlagValues)

    # -----------------------------------------------------------------

    def _getBitmask(self, iFlags):
        # DES: Determine the bitmask corresponding to iFlags.
        # INP: iFlags (int list): List of flag values that specify bits.
        # OUT: Numeric.array of the same typecode as self._aFlags,
        #      shape(1,), and value such that bits corresponding
        #      to iFlags are 1, all others 0
        typecode = self._aFlags.typecode()
        zero = Numeric.zeros(1,typecode=typecode)
        one  = Numeric.ones(1,typecode=typecode)

        bitmask = zero

        if iFlags==[]:
            iFlags = self._validFlagValues
        if type(iFlags)==type(1):
            iFlags = [iFlags]

        for iFlag in iFlags:
            if not iFlag in self._validFlagValues:
                raise "Invalid flag %d" % iFlag
            bitmask |= (one<<(iFlag-1)).astype(typecode)

        return bitmask


class FlagHandler1d32b(FlagHandler):
    """
    NAM: FlagHandler1d32b (Class)
    DES: Bitwise manipulation of 1-dim flag arrays of type Int32
    """

    # -----------------------------------------------------------------

    def __init__(self, aFlags):
        FlagHandler.__init__(self, aFlags)

    # -----------------------------------------------------------------

    def setAll(self, iFlags=[]):
        """
        NAM: FlagHandler1d32b.setAll (Method)
        DES: Sets the flag values iFlags for all elements of self._aFlags
        INP: iFlags (int list) : Flag values (see also doc string of class FlagHandler)
        """
        bitmask = self._getBitmask(iFlags)
        fFlagHandler.flags1d32b.setall(self._aFlags, bitmask)

    # -----------------------------------------------------------------

    def unsetAll(self, iFlags=[]):
        """
        NAM: FlagHandler1d32b.unsetAll (Method)
        DES: Unsets the flag values iFlags for all elements of self._aFlags
        INP: iFlags (int list) : Flag values (see also doc string of class FlagHandler)
        """
        bitmask = self._getBitmask(iFlags)
        fFlagHandler.flags1d32b.unsetall(self._aFlags, bitmask)

    # -----------------------------------------------------------------

    def nSet(self, iFlags=[]):
        """
        NAM: FlagHandler1d32b.nSet (Method)
        DES: Returns the number of elements of self._aFlags for which
             at least one flag value specified by iFlags is set.
        INP: iFlags (int list) : Flag values (see also doc string of class FlagHandler)
        OUT: (int) : Number of elements of self._aFlags for which at least
                     one flag value specified by iFlags is set.
        """
        bitmask = self._getBitmask(iFlags)
        return fFlagHandler.flags1d32b.nset(self._aFlags, bitmask)

    # -----------------------------------------------------------------

    def nUnset(self, iFlags=[]):
        """
        NAM: FlagHandler1d32b.nUnset (Method)
        DES: Returns the number of elements of self._aFlags for which
             none of the flag values specified by iFlags is set.
        INP: iFlags (int list) : Flag values (see also doc string of class FlagHandler)
        OUT: (int) : Number of elements of self._aFlags for which none
                     of the flag values specified by iFlags is set.
        """
        bitmask = self._getBitmask(iFlags)
        return fFlagHandler.flags1d32b.nunset(self._aFlags, bitmask)

    # -----------------------------------------------------------------

    def setOnIndex(self, index, iFlags=[]):
        """
        NAM: FlagHandler1d32b.setOnIndex (Method)
        DES: Sets the flag values iFlags for a single element of
             flag array self._aFlags
        INP: index (int)       : Index of the element of self._aFlags to be set
             iFlags (int list) : Flag values (see also doc string of class FlagHandler)
        """
        bitmask = self._getBitmask(iFlags)
        fFlagHandler.flags1d32b.setonindex(self._aFlags, index+1, bitmask)

    # -----------------------------------------------------------------

    def unsetOnIndex(self, index, iFlags=[]):
        """
        NAM: FlagHandler1d32b.unsetOnIndex (Method)
        DES: Unsets the flag values iFlags for a single element of
             flag array self._aFlags
        INP: index (int)       : Index of the element of self._aFlags to be unset
             iFlags (int list) : Flag values (see also doc string of class FlagHandler)
        """
        bitmask = self._getBitmask(iFlags)
        fFlagHandler.flags1d32b.unsetonindex(self._aFlags, index+1, bitmask)

    # -----------------------------------------------------------------

    def isSetOnIndex(self, index, iFlags=[]):
        """
        NAM: FlagHandler1d32b.isSetOnIndex (Method)
        DES: Returns 1 if at least one flag value specified by iFlags
             is set for a single element of flag array self._aFlags.
        INP: index (int)       : Index of the element of self._aFlags to be set
             iFlags (int list) : Flag values (see also doc string of class FlagHandler)
        OUT: (int) : 1 if at least one flag value specified by iFlags
                     is set, 0 else.
        """
        bitmask = self._getBitmask(iFlags)
        return fFlagHandler.flags1d32b.issetonindex(self._aFlags, \
                                                    index+1, \
                                                    bitmask)
    
    # -----------------------------------------------------------------

    def isUnsetOnIndex(self, index, iFlags=[]):
        """
        NAM: FlagHandler1d32b.isSetOnIndex (Method)
        DES: Returns 1 if none of the flag values specified by iFlags
             are set for a single element of flag array self._aFlags.
        INP: index (int)       : Index of the element of self._aFlags to be set
             iFlags (int list) : Flag values (see also doc string of class FlagHandler)
        OUT: (int) : 1 if none of the flag values specified by iFlags
                     is set, 0 else.
        """
        bitmask = self._getBitmask(iFlags)
        return fFlagHandler.flags1d32b.isunsetonindex(self._aFlags, \
                                                      index+1, \
                                                      bitmask)

    # -----------------------------------------------------------------

    def setOnMask(self, aMask, iFlags=[]):
        """
        NAM: FlagHandler1d32b.setOnMask (Method)
        DES: Sets the flag values iFlags for all elements of
             flag array self._aFlags specified by aMask
        INP: aMask (Numeric.array) : Mask specifiing the elements of
                                     self._aFlags to be manipulated.
                                     The shape of aMask must be the shape
                                     of self._aFlags.
                                     Only elements of self._aFlags, for which aMask
                                     is not 0, will be manipulated.
             iFlags (int list) : Flag values (see also doc string of class FlagHandler)
        """
        if aMask.shape != self._aFlags.shape:
            raise "Invalid mask"

        bitmask = self._getBitmask(iFlags)
        fFlagHandler.flags1d32b.setonmask(self._aFlags, aMask, bitmask)

    # -----------------------------------------------------------------

    def unsetOnMask(self, aMask, iFlags=[]):
        """
        NAM: FlagHandler1d32b.setOnMask (Method)
        DES: Unsets the flag values iFlags for all elements of
             flag array self._aFlags specified by aMask
        INP: aMask (Numeric.array) : Mask specifiing the elements of
                                     self._aFlags to be manipulated.
                                     The shape of aMask must be the shape
                                     of self._aFlags.
                                     Only elements of self._aFlags, for which aMask
                                     is not 0, will be manipulated.
             iFlags (int list) : Flag values (see also doc string of class FlagHandler)
        """
        if aMask.shape != self._aFlags.shape:
            raise "Invalid mask"

        bitmask = self._getBitmask(iFlags)
        fFlagHandler.flags1d32b.unsetonmask(self._aFlags, aMask, bitmask)

    # -----------------------------------------------------------------

    def isSetMask(self, iFlags=[]):
        """
        NAM: FlagHandler1d32b.isSetMask (Method)
        DES: Returns a mask that indicates for which elements
             of self._aFlags at least one flag value specified by
             iFlags is set.
        INP: iFlags (int list) : Flag values (see also doc string of class FlagHandler)
        OUT: (Numeric.array) : Mask indicating that at least one flag value
                               specified by iFlags is set for the corresponding
                               element of self._aFlags.
        """
        bitmask = self._getBitmask(iFlags)
        return fFlagHandler.flags1d32b.issetmask(self._aFlags, bitmask)

    # -----------------------------------------------------------------

    def isUnsetMask(self, iFlags=[]):
        """
        NAM: FlagHandler1d32b.isUnsetMask (Method)
        DES: Returns a mask that indicates for which elements
             of self._aFlags none of the flag values specified by
             iFlags is set.
        INP: iFlags (int list) : Flag values (see also doc string of class FlagHandler)
        OUT: (Numeric.array) : Mask indicating that none of the flag values
                               specified by iFlags is set for the corresponding
                               element of self._aFlags.
        """
        bitmask = self._getBitmask(iFlags)
        return fFlagHandler.flags1d32b.isunsetmask(self._aFlags, bitmask)




class FlagHandler2d8b(FlagHandler):
    """
    NAM: FlagHandler2d8b (Class)
    DES: Bitwise manipulation of 2-dim flag arrays of type Int8
    """

    # -----------------------------------------------------------------

    def __init__(self, aFlags):
        FlagHandler.__init__(self, aFlags)
        self._aFlags = fUtilities.as_column_major_storage(self._aFlags)

    # -----------------------------------------------------------------

    def setAll(self, iFlags=[], dim=None, index=None):
        """
        NAM: FlagHandler2d8b.setAll (Method)
        DES: Sets the flag values iFlags for all elements of self._aFlags
        INP: iFlags (int list) : Flag values (see also doc string of class FlagHandler)
             dim, index (int)  : Specify a slice of self._aFlags on which this method operates.
                                 If None (default): Method operates on complete flag array.
        """
        bitmask = self._getBitmask(iFlags)

        if dim is None:
            fFlagHandler.flags2d8b.setall(self._aFlags, bitmask)
        elif dim==0:
            fFlagHandler.flags2d8b.setall2(self._aFlags, index+1, bitmask)
        elif dim==1:
            fFlagHandler.flags2d8b.setall1(self._aFlags, index+1, bitmask)
        else:
            raise "Invalid argument 'dim'"

    # -----------------------------------------------------------------

    def unsetAll(self, iFlags=[], dim=None, index=None):
        """
        NAM: FlagHandler2d8b.unsetAll (Method)
        DES: Unsets the flag values iFlags for all elements of self._aFlags
        INP: iFlags (int list) : Flag values (see also doc string of class FlagHandler)
             dim, index (int)  : Specify a slice of self._aFlags on which this method operates.
                                 If None (default): Method operates on complete flag array.
        """
        bitmask = self._getBitmask(iFlags)

        if dim is None:
            fFlagHandler.flags2d8b.unsetall(self._aFlags, bitmask)
        elif dim==0:
            fFlagHandler.flags2d8b.unsetall2(self._aFlags, index+1, bitmask)
        elif dim==1:
            fFlagHandler.flags2d8b.unsetall1(self._aFlags, index+1, bitmask)
        else:
            raise "Invalid argument 'dim'"

    # -----------------------------------------------------------------

    def nSet(self, iFlags=[], dim=None, index=None):
        """
        NAM: FlagHandler2d8b.nSet (Method)
        DES: Returns the number of elements of self._aFlags for which
             at least one flag value specified by iFlags is set.
        INP: iFlags (int list) : Flag values (see also doc string of class FlagHandler)
             dim, index (int)  : Specify a slice of self._aFlags on which this method operates.
                                 If None (default): Method operates on complete flag array.
        OUT: (int) : Number of elements of self._aFlags for which at least
                     one flag value specified by iFlags is set.
        """
        bitmask = self._getBitmask(iFlags)

        if dim is None:
            n = fFlagHandler.flags2d8b.nset(self._aFlags, bitmask)
        elif dim==0:
            n = fFlagHandler.flags2d8b.nset2(self._aFlags, index+1, bitmask)
        elif dim==1:
            n = fFlagHandler.flags2d8b.nset1(self._aFlags, index+1, bitmask)
        else:
            raise "Invalid argument 'dim'"

        return n

    # -----------------------------------------------------------------

    def nUnset(self, iFlags=[], dim=None, index=None):
        """
        NAM: FlagHandler2d8b.nUnset (Method)
        DES: Returns the number of elements of self._aFlags for which
             none of the flag values specified by iFlags is set.
        INP: iFlags (int list) : Flag values (see also doc string of class FlagHandler)
             dim, index (int)  : Specify a slice of self._aFlags on which this method operates.
                                 If None (default): Method operates on complete flag array.
        OUT: (int) : Number of elements of self._aFlags for which none
                     of the flag values specified by iFlags is set.
        """
        bitmask = self._getBitmask(iFlags)

        if dim is None:
            n = fFlagHandler.flags2d8b.nunset(self._aFlags, bitmask)
        elif dim==0:
            n = fFlagHandler.flags2d8b.nunset2(self._aFlags, index+1, bitmask)
        elif dim==1:
            n = fFlagHandler.flags2d8b.nunset1(self._aFlags, index+1, bitmask)
        else:
            raise "Invalid argument 'dim'"

        return n

    # -----------------------------------------------------------------

    def setOnIndex(self, index, iFlags=[]):
        """
        NAM: FlagHandler2d8b.setOnIndex (Method)
        DES: Sets the flag values iFlags for a single element of
             flag array self._aFlags
        INP: index (int)       : Index of the element of self._aFlags to be set
             iFlags (int list) : Flag values (see also doc string of class FlagHandler)
        """
        bitmask = self._getBitmask(iFlags)
        fFlagHandler.flags2d8b.setonindex(self._aFlags,
                                          index[0]+1, index[1]+1,
                                          bitmask)

    # -----------------------------------------------------------------

    def unsetOnIndex(self, index, iFlags=[]):
        """
        NAM: FlagHandler2d8b.unsetOnIndex (Method)
        DES: Sets the flag values iFlags for a single element of
             flag array self._aFlags
        INP: index (int)       : Index of the element of self._aFlags to be unset
             iFlags (int list) : Flag values (see also doc string of class FlagHandler)
        """
        bitmask = self._getBitmask(iFlags)
        fFlagHandler.flags2d8b.unsetonindex(self._aFlags,
                                            index[0]+1, index[1]+1,
                                            bitmask)

    # -----------------------------------------------------------------

    def isSetOnIndex(self, index, iFlags=[]):
        """
        NAM: FlagHandler2d8b.isSetOnIndex (Method)
        DES: Returns 1 if at least one flag value specified by iFlags
             is set for a single element of flag array self._aFlags.
        INP: index  (int)      : Index of the element of self._aFlags
             iFlags (int list) : Flag values (see also doc string of class FlagHandler)
        OUT: (int) : 1 if at least one flag value specified by iFlags
                     is set, 0 else.
        """
        bitmask = self._getBitmask(iFlags)
        return fFlagHandler.flags2d8b.issetonindex(self._aFlags,
                                                   index[0]+1, index[1]+1,
                                                   bitmask)
    
    # -----------------------------------------------------------------

    def isUnsetOnIndex(self, index, iFlags=[]):
        """
        NAM: FlagHandler2d8b.isSetOnIndex (Method)
        DES: Returns 1 if none of the flag values specified by iFlags
             are set for a single element of flag array self._aFlags.
        INP: index  (int)      : Index of the element of self._aFlags
             iFlags (int list) : Flag values (see also doc string of class FlagHandler)
        OUT: (int) : 1 if none of the flag values specified by iFlags
                     is set, 0 else.
        """
        bitmask = self._getBitmask(iFlags)
        return fFlagHandler.flags2d8b.isunsetonindex(self._aFlags,
                                                     index[0]+1, index[1]+1,
                                                     bitmask)

    # -----------------------------------------------------------------

    def setOnMask(self, aMask, iFlags=[], dim=None, index=None):
        """
        NAM: FlagHandler2d8b.setOnMask (Method)
        DES: Sets the flag values iFlags for all elements of
             flag array self._aFlags specified by aMask
        INP: aMask (Numeric.array) : Mask specifiing the elements of
                                     self._aFlags to be manipulated.
                                     The shape of aMask must be the shape
                                     of self._aFlags.
                                     Only elements of self._aFlags, for which aMask
                                     is not 0, will be manipulated.
             iFlags (int list) : Flag values (see also doc string of class FlagHandler)
             dim, index (int)  : Specify a slice of self._aFlags on which this method operates.
                                 If None (default): Method operates on complete flag array.
        """
        bitmask = self._getBitmask(iFlags)
        
        if dim is None:
            fFlagHandler.flags2d8b.setonmask(self._aFlags, aMask, bitmask)
        elif dim==0:
            fFlagHandler.flags2d8b.setonmask2(self._aFlags, aMask, index+1, bitmask)
        elif dim==1:
            fFlagHandler.flags2d8b.setonmask1(self._aFlags, aMask, index+1, bitmask)
        else:
            raise "Invalid argument 'dim'"
            

    # -----------------------------------------------------------------

    def unsetOnMask(self, aMask, iFlags=[], dim=None, index=None):
        """
        NAM: FlagHandler2d8b.setOnMask (Method)
        DES: Sets the flag values iFlags for all elements of
             flag array self._aFlags specified by aMask
        INP: aMask (Numeric.array) : Mask specifiing the elements of
                                     self._aFlags to be manipulated.
                                     The shape of aMask must be the shape
                                     of self._aFlags.
                                     Only elements of self._aFlags, for which aMask
                                     is not 0, will be manipulated.
             iFlags (int list) : Flag values (see also doc string of class FlagHandler)
             dim, index (int)  : Specify a slice of self._aFlags on which this method operates.
                                 If None (default): Method operates on complete flag array.
        """
        bitmask = self._getBitmask(iFlags)
        
        if dim is None:
            fFlagHandler.flags2d8b.unsetonmask(self._aFlags, aMask, bitmask)
        elif dim==0:
            fFlagHandler.flags2d8b.unsetonmask2(self._aFlags, aMask, index+1, bitmask)
        elif dim==1:
            fFlagHandler.flags2d8b.unsetonmask1(self._aFlags, aMask, index+1, bitmask)
        else:
            raise "Invalid argument 'dim'"

    # -----------------------------------------------------------------

    def isSetMask(self, iFlags=[], dim=None, index=None):
        """
        NAM: FlagHandler2d8b.isSetMask (Method)
        DES: Returns a mask that indicates for which elements
             of self._aFlags at least one flag value specified by
             iFlags is set.
        INP: iFlags (int list) : Flag values (see also doc string of class FlagHandler)
             dim, index (int)  : Specify a slice of self._aFlags on which this method operates.
                                 If None (default): Method operates on complete flag array.
        OUT: (Numeric.array) : Mask indicating that at least one flag value
                               specified by iFlags is set for the corresponding
                               element of self._aFlags.
        """
        bitmask = self._getBitmask(iFlags)

        if dim is None:
            aMask = fFlagHandler.flags2d8b.issetmask(self._aFlags, bitmask)
        elif dim==0:
            aMask = fFlagHandler.flags2d8b.issetmask2(self._aFlags, index+1, bitmask)

        elif dim==1:
            aMask = fFlagHandler.flags2d8b.issetmask1(self._aFlags, index+1, bitmask)
        else:
            raise "Invalid argument 'dim'"

        return aMask

    # -----------------------------------------------------------------

    def isUnsetMask(self, iFlags=[], dim=None, index=None):
        """
        NAM: FlagHandler2d8b.isUnsetMask (Method)
        DES: Returns a mask that indicates for which elements
             of self._aFlags none of the flag values specified by
             iFlags is set.
        INP: iFlags (int list) : Flag values (see also doc string of class FlagHandler)
             dim, index (int)  : Specify a slice of self._aFlags on which this method operates.
                                 If None (default): Method operates on complete flag array.
        OUT: (Numeric.array) : Mask indicating that none of the flag values
                               specified by iFlags is set for the corresponding
                               element of self._aFlags.
        """
        bitmask = self._getBitmask(iFlags)

        if dim is None:
            aMask = fFlagHandler.flags2d8b.isunsetmask(self._aFlags, bitmask)
        elif dim==0:
            aMask = fFlagHandler.flags2d8b.isunsetmask2(self._aFlags, index+1, bitmask)

        elif dim==1:
            aMask = fFlagHandler.flags2d8b.isunsetmask1(self._aFlags, index+1, bitmask)
        else:
            raise "Invalid argument 'dim'"

        return aMask
