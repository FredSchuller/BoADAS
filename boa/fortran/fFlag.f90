! Copyright (C) 2002-2006
! Max-Planck-Institut fuer Radioastronomie Bonn
! Argelander Institut fuer Astronomie
! Astronomisches Institut der Ruhr-Universität Bochum
!
! Produced for the LABOCA project
!
! This library is free software; you can redistribute it and/or modify it under
! the terms of the GNU Library General Public License as published by the Free
! Software Foundation; either version 2 of the License, or (at your option) any
! later version.
!
! This library is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU Library General Public License for more
! details.
!
! You should have received a copy of the GNU Library General Public License
! along with this library; if not, write to the Free Software Foundation, Inc.,
! 675 Massachusetts Ave, Cambridge, MA 02139, USA.  
!

subroutine flagInTime(ChanList,DataFlag,Time,above,below,flagIn,flagOut,nFlags,nChanList,nChan,nInt)
! """
! DES: flag the data w. time
! INP:
!      ChanList = flag only this channels indexes
!      Dataflag = entire flag array
!          time = time stream  
!   below,above = the limits
!        flagIn = the flag which should be tested
!       flagOut = the value of the flag to put
!     nChanList = the lenght of ChanList
!         nChan = nb. channels
!          nInt = number of time steps
! OUT: DataFlag = updated DataFlag array
!         nFlag = number of timestamp flagged
! """
!INOUT ----------------------------
  integer                             :: nFlags, nInt, nChan, nChanList
  integer, dimension(nChanList)       :: ChanList
  integer, dimension(nInt,nChan)      :: DataFlag
  real, dimension(nInt)               :: Time
  real                                :: above, below
  integer                             :: flagIn,flagOut
  !f2py intent(hide),depend(DataFlag) :: nInt=shape(DataFlag,0),nChan=shape(DataFlag,1)
  !f2py intent(hide),depend(chanList) :: nChanList=shape(ChanList)
  !f2py intent(in)                    :: ChanList, Time, below, above, flag
  !f2py intent(in,out)                :: DataFlag
  !f2py intent(hide,out)              :: nFlags
  !LOCAL ----------------------------
  integer                             :: iChannel, index, nToFlag
  logical, dimension(nInt)            :: time_mask, flag_mask
  ! ---------------------------------
  

  time_mask = (Time.ge.above).and.(Time.le.below)

  nFlags = 0
  DO iChannel=1,nChanList
     ! index in the data/flag array (remember fortran indexing)
     index = ChanList(iChannel)+1
     
     ! find the value we want to flag i.e. the flagIn within time_mask
     flag_mask = (DataFlag(:,index).eq.flagIn).and.time_mask
     nToFlag = COUNT(flag_mask)
     
     IF (nToFlag>0) THEN
        nFlags = nFlags+nToFlag
        where(flag_mask) DataFlag(:,index) = flagOut
     END IF
     
  END DO
  
END SUBROUTINE flagInTime

subroutine flagPosition(ChanList,DataFlag, AzEl, Offsets, AzElRef, radius, flag, nChanList, nChan, nInt, nFlags)
! """
! DES: flag the data around a position in sky
! INP:
!      ChanList = flag only this channels indexes
!      Dataflag = entire flag array
!         AzEl  = Reference Position of the array with...
!       Offsets = ... channel Offsets for the Used Channel only
!       AzElRef = Position  to flag...
!       radius  = .. within the given radius
!          flag = the value of the flag to put
!     nChanList = the lenght of ChanList
!         nChan = nb. channels
!          nInt = number of time steps
! OUT: DataFlag = updated DataFlag array
!         nFlag = number of timestamp flagged
! """
!INOUT ----------------------------
  integer                             :: nFlags, nInt, nChan, nChanList
  integer, dimension(nChanList)       :: ChanList
  integer, dimension(nInt,nChan)      :: DataFlag
  real,    dimension(2,nInt)          :: AzEl
  real,    dimension(2,nChan)         :: Offsets
  real,    dimension(2)               :: AzElRef
  real                                :: radius
  integer                             :: flag
  !f2py intent(hide),depend(DataFlag) :: nInt=shape(DataFlag,0),nChan=shape(DataFlag,1)
  !f2py intent(hide),depend(chanList) :: nChanList=shape(ChanList)
  !f2py intent(in)                    :: ChanList, AzEl, Offsets, AzElRef, radius, flag
  !f2py intent(in,out)                :: DataFlag
  !f2py intent(hide,out)              :: nFlags
  !LOCAL ----------------------------
  integer                             :: iChannel, index, nToFlag
  logical, dimension(nInt)            :: flag_mask
  real, dimension(nInt)               :: dist
  real,    dimension(2,nInt)          :: iAzEl
  ! ---------------------------------
  
  nFlags = 0
  DO iChannel=1,nChanList
     ! index in the data/flag array (remember fortran indexing)
     index = ChanList(iChannel)+1

     iAzEl(1,:) = AzEl(1,:)+Offsets(1,index)-AzElRef(1)
     iAzEl(2,:) = AzEl(2,:)+Offsets(2,index)-AzElRef(2)
     
     ! find out the valid channels defined in the aperture
     flag_mask = (SQRT(SUM(iAzEl**2,DIM=1)).le.radius).and.(DataFlag(:,index).eq.0)
     nToFlag = COUNT(flag_mask)

     IF (nToFlag>1) THEN
        nFlags = nFlags+nToFlag
        where(flag_mask) DataFlag(:,index) = flag
     END IF
     
  END DO
  
END SUBROUTINE flagPosition

