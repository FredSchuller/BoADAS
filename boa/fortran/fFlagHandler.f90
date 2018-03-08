!"""
!NAM: fFlagHandler.f90 (file)
!DES: Contains modules with subroutines for fast bitwise manipulation of flags.
!     Throughout the file, the nomenclature is as follows:
!
!     - aFlags(n_aFlags):
!       The flag array to be manipulated or queried.
!
!       Permitted shapes and datatypes of aFlags depend on the
!       the module
!
!     - bitmask:
!       Specifies the flag values used for manipulating or queriing the flag array.
!
!       The datatype corresponds to the datatype of aFlags and depends on the module.
!"""
module flags1d32b
!"""
!NAM: flags1d32b (module)
!DES: Bitwise manipulation of 1-dim flag arrays of type integer*4
!"""
contains
  subroutine setAll(aFlags, bitmask, n_aFlags)
    !"""
    !NAM: flags1d32b.setAll (Subroutine)
    !DES: Sets the flag values specified by bitmask for all elements of aFlags
    !INP: aFlags (integer*4(n_aFlags)) : Flag array
    !     bitmask (integer*4)          : Specifies the flag values
    !"""
    integer*4                                   :: n_aFlags
    integer*4, dimension(n_aFlags)              :: aFlags
    integer*4                                   :: bitmask
    !f2py integer*4 intent(inplace)             :: aFlags
    !f2py integer*4 intent(in)                  :: bitmask
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags = len(aFlags)

    aFlags = ior(aFlags, bitmask)
  end subroutine setAll

  subroutine unsetAll(aFlags, bitmask, n_aFlags)
    !"""
    !NAM: flags1d32b.unsetAll (Subroutine)
    !DES: Unsets the flag values specified by bitmask for all elements of aFlags
    !INP: aFlags (integer*4(n_aFlags)) : Flag array
    !     bitmask (integer*4)          : Specifies the flag values
    !"""
    integer*4                                   :: n_aFlags
    integer*4, dimension(n_aFlags)              :: aFlags
    integer*4                                   :: bitmask
    !f2py integer*4 intent(inplace)             :: aFlags
    !f2py integer*4 intent(in)                  :: bitmask
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags = len(aFlags)

    aFlags = iand(aFlags, not(bitmask))
  end subroutine unsetAll

  subroutine nSet(aFlags, bitmask, n_aFlags, countSet)
    !"""
    !NAM: flags1d32b.nSet (Subroutine)
    !DES: Returns the number of elements of aFlags for which
    !     at least one flag value specified by bitmask is set.
    !INP: aFlags (integer*4(n_aFlags)) : Flag array
    !     bitmask (integer*4)          : Specifies the flag values
    !OUT: countSet(integer*4)          : Number of elements of aFlags for which at least
    !                                    one flag value specified by bitmask is set.
    !"""
    integer*4                                   :: n_aFlags
    integer*4, dimension(n_aFlags)              :: aFlags
    integer*4                                   :: bitmask
    integer*4                                   :: countSet
    !f2py integer*4 intent(inplace)             :: aFlags
    !f2py integer*4 intent(in)                  :: bitmask
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags = len(aFlags)
    !f2py integer*4 intent(out)                 :: countSet

    countSet = count( iand(aFlags, bitmask) /= 0 )
  end subroutine nSet

  subroutine nUnset(aFlags, bitmask, n_aFlags, countUnset)
    !"""
    !NAM: flags1d32b.nUnset (Subroutine)
    !DES: Returns the number of elements of aFlags for which
    !     none of the flag values specified by bitmask is set.
    !INP: aFlags (integer*4(n_aFlags)) : Flag array
    !     bitmask (integer*4)          : Specifies the flag values
    !OUT: countUnset(integer*4)        : Number of elements of aFlags for which none of
    !                                    the flag values specified by bitmask is set.
    !"""
    integer*4                                   :: n_aFlags
    integer*4, dimension(n_aFlags)              :: aFlags
    integer*4                                   :: bitmask
    integer*4                                   :: countUnset
    !f2py integer*4 intent(inplace)             :: aFlags
    !f2py integer*4 intent(in)                  :: bitmask
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags = len(aFlags)
    !f2py integer*4 intent(out)                 :: countUnset

    countUnset = count( iand(aFlags, bitmask) == 0 )
  end subroutine nUnset

  subroutine setOnIndex(aFlags, index, bitmask, n_aFlags)
    !"""
    !NAM: flags1d32b.setOnIndex (Subroutine)
    !DES: Sets the flag values specified by bitmask for a single element of
    !     the flag array aFlags
    !INP: aFlags (integer*4(n_aFlags)) : Flag array
    !     index (integer*4)            : Index of the element of aFlags to be set
    !     bitmask (integer*4)          : Specifies the flag values
    !"""
    integer*4                                   :: n_aFlags
    integer*4, dimension(n_aFlags)              :: aFlags
    integer*4                                   :: index
    integer*4                                   :: bitmask
    !f2py integer*4 intent(inplace)             :: aFlags
    !f2py integer*4 intent(in)                  :: index
    !f2py integer*4 intent(in)                  :: bitmask
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags = len(aFlags)

    aFlags(index) = ior(aFlags(index), bitmask)
  end subroutine setOnIndex

  subroutine unsetOnIndex(aFlags, index, bitmask, n_aFlags)
    !"""
    !NAM: flags1d32b.unsetOnIndex (Subroutine)
    !DES: Unsets the flag values specified by bitmask for a single element of
    !     the flag array aFlags
    !INP: aFlags (integer*4(n_aFlags)) : Flag array
    !     index (integer*4)            : Index of the element of aFlags to be unset
    !     bitmask (integer*4)          : Specifies the flag values
    !"""
    integer*4                                   :: n_aFlags
    integer*4, dimension(n_aFlags)              :: aFlags
    integer*4                                   :: index
    integer*4                                   :: bitmask
    !f2py integer*4 intent(inplace)             :: aFlags
    !f2py integer*4 intent(in)                  :: index
    !f2py integer*4 intent(in)                  :: bitmask
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags = len(aFlags)

    aFlags(index) = iand(aFlags(index), not(bitmask))
  end subroutine unsetOnIndex

  subroutine isSetOnIndex(aFlags, index, bitmask, n_aFlags, isSet)
    !"""
    !NAM: flags1d32b.isSetOnIndex (Subroutine)
    !DES: Returns 1 if at least one flag value specified by bitmask
    !     is set for a single element of flag array aFlags.
    !INP: aFlags (integer*4(n_aFlags)) : Flag array
    !     index (integer*4)            : Index of the element of aFlags to be unset
    !     bitmask (integer*4)          : Specifies the flag values
    !OUT: isSet (integer*1)            : 1 if at least one flag value specified by bitmask
    !                                    is set, 0 else.
    !"""
    integer*4                                   :: n_aFlags
    integer*4, dimension(n_aFlags)              :: aFlags
    integer*4                                   :: index
    integer*4                                   :: bitmask
    integer*1                                   :: isSet
    !f2py integer*4 intent(inplace)             :: aFlags
    !f2py integer*4 intent(in)                  :: index
    !f2py integer*4 intent(in)                  :: bitmask
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags = len(aFlags)
    !f2py integer*1 intent(out)                 :: isSet

    if (iand(aFlags(index), bitmask) /= 0) then
       isSet = 1
    else
       isSet = 0
    end if
  end subroutine isSetOnIndex

  subroutine isUnsetOnIndex(aFlags, index, bitmask, n_aFlags, isUnset)
    !"""
    !NAM: flags1d32b.isUnsetOnIndex (Subroutine)
    !DES: Returns 1 if none of the flag values specified by bitmask
    !     is set for a single element of flag array aFlags.
    !INP: aFlags (integer*4(n_aFlags)) : Flag array
    !     index (integer*4)            : Index of the element of aFlags to be unset
    !     bitmask (integer*4)          : Specifies the flag values
    !OUT: isUnset (integer*1)          : 1 if none of the flag values specified by bitmask
    !                                    is set, 0 else.
    !"""
    integer*4                                   :: n_aFlags
    integer*4, dimension(n_aFlags)              :: aFlags
    integer*4                                   :: index
    integer*4                                   :: bitmask
    integer*1                                   :: isUnset
    !f2py integer*4 intent(inplace)             :: aFlags
    !f2py integer*4 intent(in)                  :: index
    !f2py integer*4 intent(in)                  :: bitmask
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags = len(aFlags)
    !f2py integer*1 intent(out)                 :: isUnset

    if (iand(aFlags(index), bitmask) /= 0) then
       isUnset = 0
    else
       isUnset = 1
    end if
  end subroutine isUnsetOnIndex

  subroutine setOnMask(aFlags, aMask, bitmask, n_aFlags)
    !"""
    !NAM: flags1d32b.setOnMask (Subroutine)
    !DES: Sets the flag values specified by bitmask for all elements of
    !     the flag array aFlags specified by aMask
    !INP: aFlags (integer*4(n_aFlags)) : Flag array
    !     aMask (integer*1(n_aFlags))  : Mask specifiing teh elements of
    !                                    aFlags to be manipulated.
    !     index (integer*4)            : Index of the element of aFlags to be set
    !     bitmask (integer*4)          : Specifies the flag values
    !"""
    integer*4                                   :: n_aFlags
    integer*4, dimension(n_aFlags)              :: aFlags
    integer*1, dimension(n_aFlags)              :: aMask
    integer*4                                   :: bitmask
    !f2py integer*4 intent(inplace)             :: aFlags
    !f2py integer*1 intent(in)                  :: aMask
    !f2py integer*4 intent(in)                  :: bitmask
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags = len(aFlags)

    where( aMask /= 0 )
       aFlags = ior(aFlags, bitmask)
    end where
  end subroutine setOnMask

  subroutine unsetOnMask(aFlags, aMask, bitmask, n_aFlags)
    !"""
    !NAM: flags1d32b.unsetOnMask (Subroutine)
    !DES: Unsets the flag values specified by bitmask for all elements of
    !     the flag array aFlags specified by aMask
    !INP: aFlags (integer*4(n_aFlags)) : Flag array
    !     aMask (integer*1(n_aFlags))  : Mask specifiing teh elements of
    !                                    aFlags to be manipulated.
    !     index (integer*4)            : Index of the element of aFlags to be set
    !     bitmask (integer*4)          : Specifies the flag values
    !"""
    integer*4                                   :: n_aFlags
    integer*4, dimension(n_aFlags)              :: aFlags
    integer*1, dimension(n_aFlags)              :: aMask
    integer*4                                   :: bitmask
    !f2py integer*4 intent(inplace)             :: aFlags
    !f2py integer*1 intent(in)                  :: aMask
    !f2py integer*4 intent(in)                  :: bitmask
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags = len(aFlags)

    where( aMask /= 0 )
       aFlags = iand(aFlags, not(bitmask))
    end where
  end subroutine unsetOnMask

  subroutine isSetMask(aFlags, bitmask, n_aFlags, aMask)
    !"""
    !NAM: flags1d32b.isSetMask (Subroutine)
    !DES: Returns a mask that indicates for which elements
    !     of aFlags at least one of the flag values specified by
    !     bitmask is set.
    !INP: aFlags (integer*4(n_aFlags)) : Flag array
    !     bitmask (integer*4)          : Specifies the flag values
    !OUT: aMask (integer*1(n_aFlags))  : Mask indicating that at least one of the flag values
    !                                    specified by bitmask is set for the corresponding
    !                                    element of aFlags.
    !"""
    integer*4                                   :: n_aFlags
    integer*4, dimension(n_aFlags)              :: aFlags
    integer*4                                   :: bitmask
    integer*1, dimension(n_aFlags)              :: aMask
    !f2py integer*4 intent(in)                  :: aFlags
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags = len(aFlags)
    !f2py integer*4 intent(in)                  :: bitmask
    !f2py intent(out)                           :: aMask

    aMask = (iand(aFlags, bitmask) /= 0)
    ! Necessary since conversion from logical to integer*1 may yield -1
    aMask = abs(aMask)
  end subroutine isSetMask

  subroutine isUnsetMask(aFlags, bitmask, n_aFlags, aMask)
    !"""
    !NAM: flags1d32b.isUnsetMask (Subroutine)
    !DES: Returns a mask that indicates for which elements
    !     of aFlags none of the flag values specified by
    !     bitmask is set.
    !INP: aFlags (integer*4(n_aFlags)) : Flag array
    !     bitmask (integer*4)          : Specifies the flag values
    !OUT: aMask (integer*1(n_aFlags))  : Mask indicating that none of the flag values
    !                                    specified by bitmask is set for the corresponding
    !                                    element of aFlags.
    !"""
    integer*4                                   :: n_aFlags
    integer*4, dimension(n_aFlags)              :: aFlags
    integer*4                                   :: bitmask
    integer*1, dimension(n_aFlags)              :: aMask
    !f2py integer*4 intent(in)                  :: aFlags
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags = len(aFlags)
    !f2py integer*4 intent(in)                  :: bitmask
    !f2py intent(out)                           :: aMask

    aMask = (iand(aFlags, bitmask) == 0)
    ! Necessary since conversion from logical to integer*1 may yield -1
    aMask = abs(aMask)
  end subroutine isUnsetMask

end module flags1d32b

module flags2d8b
!"""
!NAM: flags1d32b (module)
!DES: Bitwise manipulation of 1-dim flag arrays of type integer*4
!"""
contains
  subroutine setAll(aFlags, bitmask, n_aFlags1, n_aFlags2)
    !"""
    !NAM: flags2d8b.setAll (Subroutine)
    !DES: Sets the flag values specified by bitmask for all elements of aFlags
    !INP: aFlags (integer*1(n_aFlags1, n_aFlags2)) : Flag array
    !     bitmask (integer*1)          : Specifies the flag values
    !"""
    integer*4                                   :: n_aFlags1, n_aFlags2
    integer*1, dimension(n_aFlags1, n_aFlags2)  :: aFlags
    integer*1                                   :: bitmask
    !f2py integer*1 intent(inplace)             :: aFlags
    !f2py integer*1 intent(in)                  :: bitmask
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags1 = shape(aFlags,0)
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags2 = shape(aFlags,1)

    aFlags = ior(aFlags, bitmask)
  end subroutine setAll

  subroutine setAll1(aFlags, index, bitmask, n_aFlags1, n_aFlags2)
    !"""
    !NAM: flags2d8b.setAll1 (Subroutine)
    !DES: flags2d8b.setAll restricted to slice in dimension 1
    !INP: index (integer*4) : Specifies slice
    !"""
    integer*4                                   :: n_aFlags1, n_aFlags2
    integer*4                                   :: index
    integer*1, dimension(n_aFlags1, n_aFlags2)  :: aFlags
    integer*1                                   :: bitmask
    !f2py integer*1 intent(inplace)             :: aFlags
    !f2py integer*4 intent(in)                  :: index
    !f2py integer*1 intent(in)                  :: bitmask
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags1 = shape(aFlags,0)
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags2 = shape(aFlags,1)

    aFlags(:,index) = ior(aFlags(:,index), bitmask)
  end subroutine setAll1

  subroutine setAll2(aFlags, index, bitmask, n_aFlags1, n_aFlags2)
    !"""
    !NAM: flags2d8b.setAll2 (Subroutine)
    !DES: flags2d8b.setAll restricted to slice in dimension 2
    !INP: index (integer*4) : Specifies slice
    !"""
    integer*4                                   :: n_aFlags1, n_aFlags2
    integer*4                                   :: index
    integer*1, dimension(n_aFlags1, n_aFlags2)  :: aFlags
    integer*1                                   :: bitmask
    !f2py integer*1 intent(inplace)             :: aFlags
    !f2py integer*4 intent(in)                  :: index
    !f2py integer*1 intent(in)                  :: bitmask
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags1 = shape(aFlags,0)
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags2 = shape(aFlags,1)

    aFlags(index,:) = ior(aFlags(index,:), bitmask)
  end subroutine setAll2

  subroutine unsetAll(aFlags, bitmask, n_aFlags1, n_aFlags2)
    !"""
    !NAM: flags2d8b.unsetAll (Subroutine)
    !DES: Unsets the flag values specified by bitmask for all elements of aFlags
    !INP: aFlags (integer*1(n_aFlags)) : Flag array
    !     bitmask (integer*4)          : Specifies the flag values
    !"""
    integer*4                                   :: n_aFlags1, n_aFlags2
    integer*1, dimension(n_aFlags1, n_aFlags2)  :: aFlags
    integer*1                                   :: bitmask
    !f2py integer*1 intent(inplace)             :: aFlags
    !f2py integer*1 intent(in)                  :: bitmask
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags1 = shape(aFlags,0)
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags2 = shape(aFlags,1)

    aFlags = iand(aFlags, not(bitmask))
  end subroutine unsetAll

  subroutine unsetAll1(aFlags, index, bitmask, n_aFlags1, n_aFlags2)
    !"""
    !NAM: flags2d8b.unsetAll1 (Subroutine)
    !DES: flags2d8b.unsetAll restricted to slice in dimension 1
    !INP: index (integer*4) : Specifies slice
    !"""
    integer*4                                   :: n_aFlags1, n_aFlags2
    integer*4                                   :: index
    integer*1, dimension(n_aFlags1, n_aFlags2)  :: aFlags
    integer*1                                   :: bitmask
    !f2py integer*1 intent(inplace)             :: aFlags
    !f2py integer*4 intent(in)                  :: index
    !f2py integer*1 intent(in)                  :: bitmask
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags1 = shape(aFlags,0)
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags2 = shape(aFlags,1)

    aFlags(:,index) = iand(aFlags(:,index), not(bitmask))
  end subroutine unsetAll1

  subroutine unsetAll2(aFlags, index, bitmask, n_aFlags1, n_aFlags2)
    !"""
    !NAM: flags2d8b.unsetAll1 (Subroutine)
    !DES: flags2d8b.unsetAll restricted to slice in dimension 1
    !INP: index (integer*4) : Specifies slice
    !"""
    integer*4                                   :: n_aFlags1, n_aFlags2
    integer*4                                   :: index
    integer*1, dimension(n_aFlags1, n_aFlags2)  :: aFlags
    integer*1                                   :: bitmask
    !f2py integer*1 intent(inplace)             :: aFlags
    !f2py integer*4 intent(in)                  :: index
    !f2py integer*1 intent(in)                  :: bitmask
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags1 = shape(aFlags,0)
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags2 = shape(aFlags,1)

    aFlags(index,:) = iand(aFlags(index,:), not(bitmask))
  end subroutine unsetAll2

  subroutine nSet(aFlags, bitmask, n_aFlags1, n_aFlags2, countSet)
    !"""
    !NAM: flags2d8b.nSet (Subroutine)
    !DES: Returns the number of elements of aFlags for which
    !     at least one flag value specified by bitmask is set.
    !INP: aFlags (integer*1(n_aFlags1, n_aFlags2)) : Flag array
    !     bitmask (integer*1)          : Specifies the flag values
    !OUT: countSet(integer*4)          : Number of elements of aFlags for which at least
    !                                    one flag value specified by bitmask is set.
    !"""
    integer*4                                   :: n_aFlags1, n_aFlags2
    integer*1, dimension(n_aFlags1, n_aFlags2)  :: aFlags
    integer*4                                   :: bitmask
    integer*4                                   :: countSet
    !f2py integer*1 intent(inplace)             :: aFlags
    !f2py integer*4 intent(in)                  :: bitmask
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags1 = shape(aFlags,0)
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags2 = shape(aFlags,1)
    !f2py integer*4 intent(out)                 :: countSet

    countSet = count( iand(aFlags, bitmask) /= 0 )
  end subroutine nSet

  subroutine nSet1(aFlags, index, bitmask, n_aFlags1, n_aFlags2, countUnset)
    !"""
    !NAM: flags2d8b.nSet1 (Subroutine)
    !DES: flags2d8b.nSet restricted to slice in dimension 1
    !INP: index (integer*4) : Specifies slice
    !"""
    integer*4                                   :: n_aFlags1, n_aFlags2
    integer*1, dimension(n_aFlags1, n_aFlags2)  :: aFlags
    integer*4                                   :: index
    integer*4                                   :: bitmask
    integer*4                                   :: countUnset
    !f2py integer*1 intent(inplace)             :: aFlags
    !f2py integer*4 intent(in)                  :: index
    !f2py integer*4 intent(in)                  :: bitmask
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags1 = shape(aFlags,0)
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags2 = shape(aFlags,1)
    !f2py integer*4 intent(out)                 :: countUnset

    countUnset = count( iand(aFlags(:,index), bitmask) /= 0 )
  end subroutine nSet1

  subroutine nSet2(aFlags, index, bitmask, n_aFlags1, n_aFlags2, countUnset)
    !"""
    !NAM: flags2d8b.nSet2 (Subroutine)
    !DES: flags2d8b.nSet restricted to slice in dimension 2
    !INP: index (integer*4) : Specifies slice
    !"""
    integer*4                                   :: n_aFlags1, n_aFlags2
    integer*1, dimension(n_aFlags1, n_aFlags2)  :: aFlags
    integer*4                                   :: index
    integer*4                                   :: bitmask
    integer*4                                   :: countUnset
    !f2py integer*1 intent(inplace)             :: aFlags
    !f2py integer*4 intent(in)                  :: index
    !f2py integer*4 intent(in)                  :: bitmask
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags1 = shape(aFlags,0)
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags2 = shape(aFlags,1)
    !f2py integer*4 intent(out)                 :: countUnset

    countUnset = count( iand(aFlags(index,:), bitmask) /= 0 )
  end subroutine nSet2


  subroutine nUnset(aFlags, bitmask, n_aFlags1, n_aFlags2, countUnset)
    !"""
    !NAM: flags2d8b.nUnset (Subroutine)
    !DES: Returns the number of elements of aFlags for which
    !     none of the flag values specified by bitmask is set.
    !INP: aFlags (integer*1(n_aFlags1, n_aFlags2)) : Flag array
    !     bitmask (integer*1)          : Specifies the flag values
    !OUT: countUnset(integer*4)        : Number of elements of aFlags for which none of
    !                                    the flag values specified by bitmask is set.
    !"""
    integer*4                                   :: n_aFlags1, n_aFlags2
    integer*1, dimension(n_aFlags1, n_aFlags2)  :: aFlags
    integer*4                                   :: bitmask
    integer*4                                   :: countUnset
    !f2py integer*1 intent(inplace)             :: aFlags
    !f2py integer*4 intent(in)                  :: bitmask
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags1 = shape(aFlags,0)
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags2 = shape(aFlags,1)
    !f2py integer*4 intent(out)                 :: countUnset

    countUnset = count( iand(aFlags, bitmask) == 0 )
  end subroutine nUnset

  subroutine nUnset1(aFlags, index, bitmask, n_aFlags1, n_aFlags2, countUnset)
    !"""
    !NAM: flags2d8b.nUnset1 (Subroutine)
    !DES: flags2d8b.nUnset restricted to slice in dimension 1
    !INP: index (integer*4) : Specifies slice
    !"""
    integer*4                                   :: n_aFlags1, n_aFlags2
    integer*1, dimension(n_aFlags1, n_aFlags2)  :: aFlags
    integer*4                                   :: index
    integer*4                                   :: bitmask
    integer*4                                   :: countUnset
    !f2py integer*1 intent(inplace)             :: aFlags
    !f2py integer*4 intent(in)                  :: index
    !f2py integer*4 intent(in)                  :: bitmask
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags1 = shape(aFlags,0)
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags2 = shape(aFlags,1)
    !f2py integer*4 intent(out)                 :: countUnset

    countUnset = count( iand(aFlags(:,index), bitmask) == 0 )
  end subroutine nUnset1

  subroutine nUnset2(aFlags, index, bitmask, n_aFlags1, n_aFlags2, countUnset)
    !"""
    !NAM: flags2d8b.nUnset2 (Subroutine)
    !DES: flags2d8b.nUnset restricted to slice in dimension 2
    !INP: index (integer*4) : Specifies slice
    !"""
    integer*4                                   :: n_aFlags1, n_aFlags2
    integer*1, dimension(n_aFlags1, n_aFlags2)  :: aFlags
    integer*4                                   :: index
    integer*4                                   :: bitmask
    integer*4                                   :: countUnset
    !f2py integer*1 intent(inplace)             :: aFlags
    !f2py integer*4 intent(in)                  :: index
    !f2py integer*4 intent(in)                  :: bitmask
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags1 = shape(aFlags,0)
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags2 = shape(aFlags,1)
    !f2py integer*4 intent(out)                 :: countUnset

    countUnset = count( iand(aFlags(index,:), bitmask) == 0 )
  end subroutine nUnset2

  subroutine setOnIndex(aFlags, index1, index2, bitmask, n_aFlags1, n_aFlags2)
    !"""
    !NAM: flags2d8b.setOnIndex (Subroutine)
    !DES: Sets the flag values specified by bitmask for a single element of
    !     the flag array aFlags
    !INP: aFlags (integer*1(n_aFlags1, n_aFlags2)) : Flag array
    !     index1, index2 (integer*4)   : Indices of the element of aFlags to be set
    !     bitmask (integer*1)          : Specifies the flag values
    !"""
    integer*4                                   :: n_aFlags1, n_aFlags2
    integer*1, dimension(n_aFlags1, n_aFlags2)  :: aFlags
    integer*4                                   :: index1, index2
    integer*1                                   :: bitmask
    !f2py integer*1 intent(inplace)             :: aFlags
    !f2py integer*4 intent(in)                  :: index1, index2
    !f2py integer*1 intent(in)                  :: bitmask
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags1 = shape(aFlags,0)
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags2 = shape(aFlags,1)

    aFlags(index1, index2) = ior(aFlags(index1, index2), bitmask)
  end subroutine setOnIndex

  subroutine unsetOnIndex(aFlags, index1, index2, bitmask, n_aFlags1, n_aFlags2)
    !"""
    !NAM: flags2d8b.unsetOnIndex (Subroutine)
    !DES: Unsets the flag values specified by bitmask for a single element of
    !     the flag array aFlags
    !INP: aFlags (integer*1(n_aFlags1, n_aFlags2)) : Flag array
    !     index1, index2 (integer*4)   : Indices of the element of aFlags to be set
    !     bitmask (integer*1)          : Specifies the flag values
    !"""
    integer*4                                   :: n_aFlags1, n_aFlags2
    integer*1, dimension(n_aFlags1, n_aFlags2)  :: aFlags
    integer*4                                   :: index1, index2
    integer*1                                   :: bitmask
    !f2py integer*1 intent(inplace)             :: aFlags
    !f2py integer*4 intent(in)                  :: index1, index2
    !f2py integer*1 intent(in)                  :: bitmask
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags1 = shape(aFlags,0)
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags2 = shape(aFlags,1)

    aFlags(index1, index2) = iand(aFlags(index1, index2), not(bitmask))
  end subroutine unsetOnIndex

  subroutine isSetOnIndex(aFlags, index1, index2, bitmask, n_aFlags1, n_aFlags2, isSet)
    !"""
    !NAM: flags2d8b.isSetOnIndex (Subroutine)
    !DES: Returns 1 if at least one flag value specified by bitmask
    !     is set for a single element of flag array aFlags.
    !INP: aFlags (integer*1(n_aFlags1, n_aFlags2)) : Flag array
    !     index1, index2 (integer*4)   : Indices of the element of aFlags to be set
    !     bitmask (integer*1)          : Specifies the flag values
    !OUT: isSet (integer*1)            : 1 if at least one flag value specified by bitmask
    !                                    is set, 0 else.
    !"""
    integer*4                                   :: n_aFlags1, n_aFlags2
    integer*1, dimension(n_aFlags1, n_aFlags2)  :: aFlags
    integer*4                                   :: index1, index2
    integer*1                                   :: bitmask
    integer*1                                   :: isSet
    !f2py integer*1 intent(inplace)             :: aFlags
    !f2py integer*4 intent(in)                  :: index1, index2
    !f2py integer*1 intent(in)                  :: bitmask
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags1 = shape(aFlags,0)
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags2 = shape(aFlags,1)
    !f2py integer*1 intent(out)                 :: isSet

    if (iand(aFlags(index1, index2), bitmask) /= 0) then
       isSet = 1
    else
       isSet = 0
    end if
  end subroutine isSetOnIndex

  subroutine isUnsetOnIndex(aFlags, index1, index2, bitmask, n_aFlags1, n_aFlags2, isUnset)
    !"""
    !NAM: flags2d8b.isUnsetOnIndex (Subroutine)
    !DES: Returns 1 if none of the flag values specified by bitmask
    !     is set for a single element of flag array aFlags.
    !INP: aFlags (integer*1(n_aFlags1, n_aFlags2)) : Flag array
    !     index1, index2 (integer*4)   : Indices of the element of aFlags to be set
    !     bitmask (integer*1)          : Specifies the flag values
    !OUT: isUnset (integer*1)          : 1 if none of the flag values specified by bitmask
    !                                    is set, 0 else.
    !"""
    integer*4                                   :: n_aFlags1, n_aFlags2
    integer*1, dimension(n_aFlags1, n_aFlags2)  :: aFlags
    integer*4                                   :: index1, index2
    integer*1                                   :: bitmask
    integer*1                                   :: isUnset
    !f2py integer*1 intent(inplace)             :: aFlags
    !f2py integer*4 intent(in)                  :: index1, index2
    !f2py integer*1 intent(in)                  :: bitmask
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags1 = shape(aFlags,0)
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags2 = shape(aFlags,1)
    !f2py integer*1 intent(out)                 :: isUnset

    if (iand(aFlags(index1, index2), bitmask) /= 0) then
       isUnset = 0
    else
       isUnset = 1
    end if
  end subroutine isUnsetOnIndex

  subroutine setOnMask(aFlags, aMask, bitmask, n_aFlags1, n_aFlags2)
    !"""
    !NAM: flags2d8b.setOnMask (Subroutine)
    !DES: Sets the flag values specified by bitmask for all elements of
    !     the flag array aFlags specified by aMask
    !INP: aFlags (integer*1(n_aFlags1, n_aFlags2)) : Flag array
    !     aMask (integer*1(n_aFlags1, n_aFlags2))  : Mask specifiing the elements of
    !                                    aFlags to be manipulated.
    !     bitmask (integer*1)          : Specifies the flag values
    !"""
    integer*4                                   :: n_aFlags1, n_aFlags2
    integer*1, dimension(n_aFlags1, n_aFlags2)  :: aFlags
    integer*1, dimension(n_aFlags1, n_aFlags2)  :: aMask
    integer*1                                   :: bitmask
    !f2py integer*1 intent(inplace)             :: aFlags
    !f2py integer*1 intent(in)                  :: aMask
    !f2py integer*1 intent(in)                  :: bitmask
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags1 = shape(aFlags,0)
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags2 = shape(aFlags,1)

    where( aMask /= 0 )
       aFlags = ior(aFlags, bitmask)
    end where
  end subroutine setOnMask

  subroutine setOnMask1(aFlags, aMask, index, bitmask, n_aFlags1, n_aFlags2)
    !"""
    !NAM: flags2d8b.setOnMask1 (Subroutine)
    !DES: flags2d8b.setOnMask restricted to slice in dimension 1
    !INP: index (integer*4) : Specifies slice
    !"""
    integer*4                                   :: n_aFlags1, n_aFlags2
    integer*1, dimension(n_aFlags1, n_aFlags2)  :: aFlags
    integer*1, dimension(n_aFlags1)             :: aMask
    integer*4                                   :: index
    integer*1                                   :: bitmask
    !f2py integer*1 intent(inplace)             :: aFlags
    !f2py integer*1 intent(in)                  :: aMask
    !f2py integer*4 intent(in)                  :: index
    !f2py integer*1 intent(in)                  :: bitmask
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags1 = shape(aFlags,0)
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags2 = shape(aFlags,1)

    where( aMask /= 0 )
       aFlags(:,index) = ior(aFlags(:,index), bitmask)
    end where
  end subroutine setOnMask1

  subroutine setOnMask2(aFlags, aMask, index, bitmask, n_aFlags1, n_aFlags2)
    !"""
    !NAM: flags2d8b.setOnMask2 (Subroutine)
    !DES: flags2d8b.setOnMask restricted to slice in dimension 2
    !INP: index (integer*4) : Specifies slice
    !"""
    integer*4                                   :: n_aFlags1, n_aFlags2
    integer*1, dimension(n_aFlags1, n_aFlags2)  :: aFlags
    integer*1, dimension(n_aFlags2)             :: aMask
    integer*4                                   :: index
    integer*1                                   :: bitmask
    !f2py integer*1 intent(inplace)             :: aFlags
    !f2py integer*1 intent(in)                  :: aMask
    !f2py integer*4 intent(in)                  :: index
    !f2py integer*1 intent(in)                  :: bitmask
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags1 = shape(aFlags,0)
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags2 = shape(aFlags,1)

    where( aMask /= 0 )
       aFlags(index,:) = ior(aFlags(index,:), bitmask)
    end where
  end subroutine setOnMask2

  subroutine unsetOnMask(aFlags, aMask, bitmask, n_aFlags1, n_aFlags2)
    !"""
    !NAM: flags2d8b.unsetOnMask (Subroutine)
    !DES: Unsets the flag values specified by bitmask for all elements of
    !     the flag array aFlags specified by aMask
    !INP: aFlags (integer*1(n_aFlags1, n_aFlags2)) : Flag array
    !     aMask (integer*1(n_aFlags1, n_aFlags2))  : Mask specifiing the elements of
    !                                    aFlags to be manipulated.
    !     bitmask (integer*1)          : Specifies the flag values
    !"""
    integer*4                                   :: n_aFlags1, n_aFlags2
    integer*1, dimension(n_aFlags1, n_aFlags2)  :: aFlags
    integer*1, dimension(n_aFlags1, n_aFlags2)  :: aMask
    integer*1                                   :: bitmask
    !f2py integer*1 intent(inplace)             :: aFlags
    !f2py integer*1 intent(in)                  :: aMask
    !f2py integer*1 intent(in)                  :: bitmask
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags1 = shape(aFlags,0)
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags2 = shape(aFlags,1)

    where( aMask /= 0 )
       aFlags = iand(aFlags, not(bitmask))
    end where
  end subroutine unsetOnMask

  subroutine unsetOnMask1(aFlags, aMask, index, bitmask, n_aFlags1, n_aFlags2)
    !"""
    !NAM: flags2d8b.unsetOnMask1 (Subroutine)
    !DES: flags2d8b.unsetOnMask restricted to slice in dimension 1
    !INP: index (integer*4) : Specifies slice
    !"""
    integer*4                                   :: n_aFlags1, n_aFlags2
    integer*1, dimension(n_aFlags1, n_aFlags2)  :: aFlags
    integer*1, dimension(n_aFlags1)             :: aMask
    integer*4                                   :: index
    integer*1                                   :: bitmask
    !f2py integer*1 intent(inplace)             :: aFlags
    !f2py integer*1 intent(in)                  :: aMask
    !f2py integer*4 intent(in)                  :: index
    !f2py integer*1 intent(in)                  :: bitmask
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags1 = shape(aFlags,0)
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags2 = shape(aFlags,1)

    where( aMask /= 0 )
       aFlags(:,index) = iand(aFlags(:,index), not(bitmask))
    end where
  end subroutine unsetOnMask1
 
  subroutine unsetOnMask2(aFlags, aMask, index, bitmask, n_aFlags1, n_aFlags2)
    !"""
    !NAM: flags2d8b.setOnMask2 (Subroutine)
    !DES: flags2d8b.setOnMask restricted to slice in dimension 2
    !INP: index (integer*4) : Specifies slice
    !"""
    integer*4                                   :: n_aFlags1, n_aFlags2
    integer*1, dimension(n_aFlags1, n_aFlags2)  :: aFlags
    integer*1, dimension(n_aFlags2)             :: aMask
    integer*4                                   :: index
    integer*1                                   :: bitmask
    !f2py integer*1 intent(inplace)             :: aFlags
    !f2py integer*1 intent(in)                  :: aMask
    !f2py integer*4 intent(in)                  :: index
    !f2py integer*1 intent(in)                  :: bitmask
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags1 = shape(aFlags,0)
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags2 = shape(aFlags,1)

    where( aMask /= 0 )
       aFlags(index,:) = iand(aFlags(index,:), not(bitmask))
    end where
  end subroutine unsetOnMask2

  subroutine isSetMask(aFlags, bitmask, n_aFlags1, n_aFlags2, aMask)
    !"""
    !NAM: flags2d8b.isSetMask (Subroutine)
    !DES: Returns a mask that indicates for which elements
    !     of aFlags at least one of the flag values specified by
    !     bitmask is set.
    !INP: aFlags (integer*1(n_aFlags1, n_aFlags2)) : Flag array
    !     bitmask (integer*1)          : Specifies the flag values
    !OUT: aMask (integer*1(n_aFlags1, n_aFlags2)) : Mask indicating that at least one of
    !                                    the flag values
    !                                    specified by bitmask is set for the corresponding
    !                                    element of aFlags.
    !"""
    integer*4                                   :: n_aFlags1, n_aFlags2
    integer*1, dimension(n_aFlags1, n_aFlags2)  :: aFlags
    integer*1                                   :: bitmask
    integer*1, dimension(n_aFlags1, n_aFlags2)  :: aMask
    !f2py integer*1 intent(in)                  :: aFlags
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags1 = shape(aFlags,0)
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags2 = shape(aFlags,1)
    !f2py integer*1 intent(in)                  :: bitmask
    !f2py intent(out)                           :: aMask

    aMask = (iand(aFlags, bitmask) /= 0)
    ! Necessary since conversion from logical to integer*1 may yield -1
    aMask = abs(aMask)
  end subroutine isSetMask

  subroutine isSetMask1(aFlags, index, bitmask, n_aFlags1, n_aFlags2, aMask)
    !"""
    !NAM: flags2d8b.isSetMask1 (Subroutine)
    !DES: flags2d8b.isSetMask restricted to slice in dimension 1
    !INP: index (integer*4) : Specifies slice
    !"""
    integer*4                                   :: n_aFlags1, n_aFlags2
    integer*1, dimension(n_aFlags1, n_aFlags2)  :: aFlags
    integer*1                                   :: bitmask
    integer*4                                   :: index
    integer*1, dimension(n_aFlags1)             :: aMask
    !f2py integer*1 intent(in)                  :: aFlags
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags1 = shape(aFlags,0)
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags2 = shape(aFlags,1)
    !f2py integer*4 intent(in)                  :: index
    !f2py integer*1 intent(in)                  :: bitmask
    !f2py intent(out)                           :: aMask

    aMask = (iand(aFlags(:,index), bitmask) /= 0)
    ! Necessary since conversion from logical to integer*1 may yield -1
    aMask = abs(aMask)
  end subroutine isSetMask1

  subroutine isSetMask2(aFlags, index, bitmask, n_aFlags1, n_aFlags2, aMask)
    !"""
    !NAM: flags2d8b.isSetMask2 (Subroutine)
    !DES: flags2d8b.isSetMask restricted to slice in dimension 2
    !INP: index (integer*4) : Specifies slice
    !"""
    integer*4                                   :: n_aFlags1, n_aFlags2
    integer*1, dimension(n_aFlags1, n_aFlags2)  :: aFlags
    integer*1                                   :: bitmask
    integer*4                                   :: index
    integer*1, dimension(n_aFlags2)             :: aMask
    !f2py integer*1 intent(in)                  :: aFlags
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags1 = shape(aFlags,0)
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags2 = shape(aFlags,1)
    !f2py integer*4 intent(in)                  :: index
    !f2py integer*1 intent(in)                  :: bitmask
    !f2py intent(out)                           :: aMask

    aMask = (iand(aFlags(index,:), bitmask) /= 0)
    ! Necessary since conversion from logical to integer*1 may yield -1
    aMask = abs(aMask)
  end subroutine isSetMask2

  subroutine isUnsetMask(aFlags, bitmask, n_aFlags1, n_aFlags2, aMask)
    !"""
    !NAM: flags2d8b.isSetMask (Subroutine)
    !DES: Returns a mask that indicates for which elements
    !     of aFlags at least one of the flag values specified by
    !     bitmask is set.
    !INP: aFlags (integer*1(n_aFlags1, n_aFlags2)) : Flag array
    !     bitmask (integer*1)          : Specifies the flag values
    !OUT: aMask (integer*1(n_aFlags1, n_aFlags2)) : Mask indicating that at least one of
    !                                    the flag values
    !                                    specified by bitmask is set for the corresponding
    !                                    element of aFlags.
    !"""
    integer*4                                   :: n_aFlags1, n_aFlags2
    integer*1, dimension(n_aFlags1, n_aFlags2)  :: aFlags
    integer*1                                   :: bitmask
    integer*1, dimension(n_aFlags1, n_aFlags2)  :: aMask
    !f2py integer*1 intent(in)                  :: aFlags
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags1 = shape(aFlags,0)
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags2 = shape(aFlags,1)
    !f2py integer*1 intent(in)                  :: bitmask
    !f2py intent(out)                           :: aMask

    aMask = (iand(aFlags, bitmask) == 0)
    ! Necessary since conversion from logical to integer*1 may yield -1
    aMask = abs(aMask)
  end subroutine isUnsetMask

  subroutine isUnsetMask1(aFlags, index, bitmask, n_aFlags1, n_aFlags2, aMask)
    !"""
    !NAM: flags2d8b.isUnsetMask1 (Subroutine)
    !DES: flags2d8b.isUnsetMask restricted to slice in dimension 1
    !INP: index (integer*4) : Specifies slice
    !"""
    integer*4                                   :: n_aFlags1, n_aFlags2
    integer*1, dimension(n_aFlags1, n_aFlags2)  :: aFlags
    integer*1                                   :: bitmask
    integer*4                                   :: index
    integer*1, dimension(n_aFlags1)             :: aMask
    !f2py integer*1 intent(in)                  :: aFlags
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags1 = shape(aFlags,0)
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags2 = shape(aFlags,1)
    !f2py integer*4 intent(in)                  :: index
    !f2py integer*1 intent(in)                  :: bitmask
    !f2py intent(out)                           :: aMask

    aMask = (iand(aFlags(:,index), bitmask) == 0)
    ! Necessary since conversion from logical to integer*1 may yield -1
    aMask = abs(aMask)
  end subroutine isUnsetMask1

  subroutine isUnsetMask2(aFlags, index, bitmask, n_aFlags1, n_aFlags2, aMask)
    !"""
    !NAM: flags2d8b.isUnsetMask2 (Subroutine)
    !DES: flags2d8b.isUnsetMask restricted to slice in dimension 2
    !INP: index (integer*4) : Specifies slice
    !"""
    integer*4                                   :: n_aFlags1, n_aFlags2
    integer*1, dimension(n_aFlags1, n_aFlags2)  :: aFlags
    integer*1                                   :: bitmask
    integer*4                                   :: index
    integer*1, dimension(n_aFlags2)             :: aMask
    !f2py integer*1 intent(in)                  :: aFlags
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags1 = shape(aFlags,0)
    !f2py integer*4 intent(hide),depend(aFlags) :: n_aFlags2 = shape(aFlags,1)
    !f2py integer*4 intent(in)                  :: index
    !f2py integer*1 intent(in)                  :: bitmask
    !f2py intent(out)                           :: aMask

    aMask = (iand(aFlags(index,:), bitmask) == 0)
    ! Necessary since conversion from logical to integer*1 may yield -1
    aMask = abs(aMask)
  end subroutine isUnsetMask2

end module flags2d8b

