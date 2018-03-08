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

!=============================================================================
subroutine arrayfitPoly_s(ChanList, Data, DataFlags, Time, Subscans, deg, poly, &
     &nChanList, nChan, nInt, nSubscans )
  ! """
  ! DES: polynomial fit to data with singular value decomposition
  !      using modified routines from Numerical Recipes.
  ! INP: 
  !      ChanList  = make the fit on this channels only
  !      Data      = the entire data array with ...
  !      DataFlags = ... corresponding flags ...
  !      Time      = ... and time, either MJD or LST or whatever
  !      Subscans  = 2d array with starting and ending index of subscans, (python indexed)
  !      deg       = polynom order 
  !
  ! OUT: poly      = the polynomial coefficient by channels and subscans
  !
  ! """
  !INOUT ----------------------------
  integer                                       :: nChan, nInt, deg
  integer, dimension(nChanList)                 :: ChanList
  real,    dimension(nInt,nChan)                :: Data
  integer, dimension(nInt,nChan)                :: DataFlags
  real,    dimension(nInt)                      :: Time
  integer, dimension(2,nSubscans)               :: Subscans
  real,    dimension(nChanList,nSubscans,deg+1) :: poly
  !f2py intent(hide),depend(ChanList)           :: nChanlist = shape(ChanList)
  !f2py intent(hide),depend(Data)               :: nInt = shape(Data,0), nChan = shape(Data,1)
  !f2py intent(hide),depend(Subscans)           :: nSubscans = shape(Subscans)
  !f2py intent(in)                              :: ChanList, Data, DataFlags, Time, Subscans, deg
  !f2py intent(out)                             :: poly
  !LOCAL ----------------------------           
  integer                                       :: iChannel, indexChan
  integer                                       :: iSubscan, startSubscan, endSubscan
  integer                                       :: nGood, subSize
  real                                          :: chisq
  real                                          :: xMax, xMin, yMax, yMin
  logical, allocatable,dimension(:)             :: flag_mask
  integer, allocatable,dimension(:)             :: iDataFlags
  real,    allocatable,dimension(:)             :: iData, iTime
  real,    allocatable,dimension(:)             :: x,y
  real,    dimension(deg+1)                     :: ipoly
  ! ---------------------------------

  DO iChannel=1, nChanList
     ! indexChan in the data/flag array (remember fortran indexing)
     indexChan = ChanList(iChannel)+1

     DO iSubscan=1, nSubscans
        startSubscan = Subscans(1,iSubscan)+1
        endSubscan   = Subscans(2,iSubscan)

        subSize = endSubscan-startSubscan+1

        ALLOCATE(flag_mask(subSize), iDataFlags(subSize), iData(subSize), iTime(subSize))
        
        iDataFlags = DataFlags(startSubscan:endSubscan,indexChan)
        iData      = Data(startSubscan:endSubscan,indexChan)
        iTime      = Time(startSubscan:endSubscan)


        ! find the valid values
        flag_mask = iDataFlags.eq.0
        nGood = COUNT(flag_mask)
        
        ! only fit if the number of good data points is greater than the degree of the polynomial
        IF (nGood>deg) THEN

           ALLOCATE(x(nGood),y(nGood))

           ! reset polynomial coefficient
           ipoly = 0
           
           x = PACK(iTime,MASK=flag_mask) 
           y = PACK(iData,MASK=flag_mask)
           
           call polyfit(x,y,deg,nGood,ipoly)

           poly(iChannel,iSubscan,:) = ipoly

           DEALLOCATE(x,y)

        ENDIF

        DEALLOCATE(flag_mask, iDataFlags, iData, iTime)

     ENDDO
  ENDDO

end subroutine arrayfitPoly_s

!=============================================================================
subroutine arrayfitPolysecant_s(ChanList, Data, DataFlags, Time, Subscans, Sec_el, deg, poly, &
     &nChanList, nChan, nInt, nSubscans )
  ! """
  ! DES: polynomial fit to data with singular value decomposition
  !      using modified routines from Numerical Recipes.
  ! INP: 
  !      ChanList  = make the fit on this channels only
  !      Data      = the entire data array with ...
  !      DataFlags = ... corresponding flags ...
  !      Time      = ... and time, either MJD or LST or whatever
  !      Subscans  = 2d array with starting and ending index of subscans, (python indexed)
  !      Sec_el    = secant of elevation
  !      deg       = polynom order 
  !
  ! OUT: poly      = the polynomial coefficient by channels and subscans
  !
  ! """
  !INOUT ----------------------------
  integer                                       :: nChan, nInt, deg
  integer, dimension(nChanList)                 :: ChanList
  real,    dimension(nInt,nChan)                :: Data
  integer, dimension(nInt,nChan)                :: DataFlags
  real,    dimension(nInt)                      :: Time
  integer, dimension(2,nSubscans)               :: Subscans
  real,    dimension(nInt)                      :: Sec_el
  real,    dimension(nChanList,nSubscans,deg+2) :: poly
  !f2py intent(hide),depend(ChanList)           :: nChanlist = shape(ChanList)
  !f2py intent(hide),depend(Data)               :: nInt = shape(Data,0), nChan = shape(Data,1)
  !f2py intent(hide),depend(Subscans)           :: nSubscans = shape(Subscans)
  !f2py intent(in)                              :: ChanList, Data, DataFlags, Time, Subscans, deg
  !f2py intent(out)                             :: poly
  !LOCAL ----------------------------           
  integer                                       :: iChannel, indexChan
  integer                                       :: iSubscan, startSubscan, endSubscan
  integer                                       :: nGood, subSize
  real                                          :: chisq
  real                                          :: xMax, xMin, yMax, yMin
  logical, allocatable,dimension(:)             :: flag_mask
  integer, allocatable,dimension(:)             :: iDataFlags
  real,    allocatable,dimension(:)             :: iData, iTime, iSecant
  real,    allocatable,dimension(:)             :: x,y,s
  real,    dimension(deg+2)                     :: ipoly
  ! ---------------------------------

  DO iChannel=1, nChanList
     ! indexChan in the data/flag array (remember fortran indexing)
     indexChan = ChanList(iChannel)+1

     DO iSubscan=1, nSubscans
        startSubscan = Subscans(1,iSubscan)+1
        endSubscan   = Subscans(2,iSubscan)

        subSize = endSubscan-startSubscan+1

        ALLOCATE(flag_mask(subSize), iDataFlags(subSize), iData(subSize), iTime(subSize), iSecant(subSize))
        
        iDataFlags = DataFlags(startSubscan:endSubscan,indexChan)
        iData      = Data(startSubscan:endSubscan,indexChan)
        iTime      = Time(startSubscan:endSubscan)
        iSecant    = Sec_el(startSubscan:endSubscan)


        ! find the valid values
        flag_mask = iDataFlags.eq.0
        nGood = COUNT(flag_mask)
        
        ! only fit if the number of good data points is greater than the degree of the polynomial 
        !   + one order for the secant
        IF (nGood>deg+1) THEN

           ALLOCATE(x(nGood),y(nGood),s(nGood))

           ! reset polynomial coefficient
           ipoly = 0
           
           x = PACK(iTime,MASK=flag_mask) 
           y = PACK(iData,MASK=flag_mask)
           s = PACK(iSecant,MASK=flag_mask)
           
           call polysecantfit(x,y,s,deg,nGood,ipoly)

           poly(iChannel,iSubscan,:) = ipoly

           DEALLOCATE(x,y,s)

        ENDIF

        DEALLOCATE(flag_mask, iDataFlags, iData, iTime, iSecant)

     ENDDO
  ENDDO

end subroutine arrayfitPolysecant_s



!=============================================================================
subroutine subtractPoly(ChanList, Data, Time, Subscans, poly, &
     &nChanList, nChan, nInt, nSubscans, nPoly )
  ! """
  ! DES: polynomial fit to data with singular value decomposition
  !      using modified routines from Numerical Recipes.
  ! INP: 
  !      ChanList  = make the fit on this channels only
  !      Data      = the entire data array with...
  !      Time      = ... time, either MJD or LST or whatever
  !      Subscans  = 2d array with starting and ending index of subscans, (python indexed)
  !      poly      = the polynomial coefficient by channels and subscans
  !
  ! OUT: poly      = the polynomial coefficient by channels and subscans
  !
  ! """
  !INOUT ----------------------------
  integer                                       :: nChanList, nChan, nInt, nSubscan, nPoly
  integer, dimension(nChanList)                 :: ChanList
  real,    dimension(nInt,nChan)                :: Data
  real,    dimension(nInt)                      :: Time
  integer, dimension(2,nSubscans)               :: Subscans
  real,    dimension(nChanList,nSubscans,nPoly) :: poly
  !f2py intent(hide),depend(ChanList)           :: nChanlist = shape(ChanList)
  !f2py intent(hide),depend(Data)               :: nInt = shape(Data,0), nChan = shape(Data,1)
  !f2py intent(hide),depend(Subscans)           :: nSubscans = shape(Subscans)
  !f2py intent(hide),depend(poly)               :: nPoly = shape(poly,2)
  !f2py intent(in)                              :: ChanList, DataFlags, Time, Subscans, poly
  !f2py intent(in,out)                          :: Data
  !LOCAL ----------------------------           
  integer                                       :: iChannel, indexChan
  real,    dimension(nSubscans,nPoly)           :: ipoly
  real,    dimension(nInt)                      :: yFit
  ! ---------------------------------

  DO iChannel=1, nChanList
     ! indexChan in the data/flag array (remember fortran indexing)
     indexChan = ChanList(iChannel)+1

     ipoly = poly(iChannel,:,:)

     ! Evaluate the polynome ... 
     CALL evalChunkedPoly(Time, Subscans, ipoly, yFit, nInt, nSubscans, nPoly)
     
     ! ... and subtract from the data
     Data(:,indexChan) = Data(:,indexChan) - yFit
     
  ENDDO

end subroutine subtractPoly


!=============================================================================
subroutine subtractPolysecant(ChanList, Data, Time, Sec, Subscans, poly, &
     &nChanList, nChan, nInt, nSubscans, nPoly )
  ! """
  ! DES: polynomial+secant subtraction given a set of parameters
  !
  ! INP: 
  !      ChanList  = make the fit on this channels only
  !      Data      = the entire data array with...
  !      Time      = ... time, either MJD or LST or whatever
  !      Sec       = secant of elevation
  !      Subscans  = 2d array with starting and ending index of subscans, (python indexed)
  !      poly      = the polynomial coefficients + the secant normalization by channels and subscans
  !
  ! OUT: poly      = the polynomial coefficient by channels and subscans
  !
  ! """
  !INOUT ----------------------------
  integer                                       :: nChanList, nChan, nInt, nSubscan, nPoly
  integer, dimension(nChanList)                 :: ChanList
  real,    dimension(nInt,nChan)                :: Data
  real,    dimension(nInt)                      :: Time, Sec
  integer, dimension(2,nSubscans)               :: Subscans
  real,    dimension(nChanList,nSubscans,nPoly) :: poly
  !f2py intent(hide),depend(ChanList)           :: nChanlist = shape(ChanList)
  !f2py intent(hide),depend(Data)               :: nInt = shape(Data,0), nChan = shape(Data,1)
  !f2py intent(hide),depend(Subscans)           :: nSubscans = shape(Subscans)
  !f2py intent(hide),depend(poly)               :: nPoly = shape(poly,2)
  !f2py intent(in)                              :: ChanList, DataFlags, Time, Subscans, poly
  !f2py intent(in,out)                          :: Data
  !LOCAL ----------------------------           
  integer                                       :: iChannel, indexChan
  real,    dimension(nSubscans,nPoly)           :: ipoly
  real,    dimension(nInt)                      :: yFit
  ! ---------------------------------

  DO iChannel=1, nChanList
     ! indexChan in the data/flag array (remember fortran indexing)
     indexChan = ChanList(iChannel)+1

     ipoly = poly(iChannel,:,:)

     ! Evaluate the polynomial ... 
     CALL evalChunkedPolysecant(Time, Sec, Subscans, ipoly, yFit, nInt, nSubscans, nPoly)
     
     ! ... and subtract from the data
     Data(:,indexChan) = Data(:,indexChan) - yFit
     
  ENDDO

end subroutine subtractPolysecant





!=============================================================================
subroutine evalChunkedPoly(x, Chunks, poly, yFit, nX, nChunks, nPoly )
  ! """
  ! DES: evaluate a polynom over several chunks
  ! INP: 
  !      x     = the absciss and ...
  !     Chunks = the start and end index of the chunks (python indexed)
  !      poly  = the polynomial coefficient for the different chunk
  !
  ! OUT:  yFit = the values of the polynomial at x
  !
  ! """
  !INOUT ----------------------------
  integer                            :: nX, nChunks, nPoly
  real,    dimension(nX)             :: x, yFit
  integer, dimension(2,nChunks)      :: Chunks
  real,    dimension(nChunks, nPoly) :: poly
  !f2py intent(hide),depend(X)       :: nX = shape(X)
  !f2py intent(hide),depend(Chunks)  :: nChunks = shape(Chunks,1)
  !f2py intent(hide),depend(poly)    :: nPoly = shape(poly,1)
  !f2py intent(in)                   :: x, Chunks, poly
  !f2py intent(hide,out)             :: yFit
  !LOCAL ----------------------------
  integer                            :: iChunk, start, end, subSize
  real                               :: xMin,xMax, yMin, yMax
  real,    dimension(nPoly)          :: ipoly
  double precision, allocatable, dimension(:)    :: xChunk, yChunk
  ! ---------------------------------

  yFix = 0
  
  DO iChunk=1, nChunks
        start = Chunks(1,iChunk)+1
        end   = Chunks(2,iChunk)

        ! allocate the proper array
        subSize = end-start+1
        ALLOCATE(xChunk(subSize),yChunk(subSize))

        xChunk = x(start:end)
        ipoly = poly(iChunk,:)

        ! Evaluate the polynom ... ( evalPoly is python indexed )
        CALL evalPoly(xChunk, ipoly, yChunk, subSize, nPoly)

        yFit(start:end) = yChunk

        DEALLOCATE(xChunk,yChunk)

     ENDDO

end subroutine evalChunkedPoly



!=============================================================================
subroutine evalChunkedPolysecant(x, sec_el, Chunks, poly, yFit, nX, nChunks, nPoly )
  ! """
  ! DES: evaluate a polynom+secant over several chunks
  ! INP: 
  !      x     = the abscissa 
  !     sec_el = secant of elevation
  !     Chunks = the start and end index of the chunks (python indexed)
  !      poly  = the polynomial coefficient for the different chunk
  !
  ! OUT:  yFit = the values of the polynomial at x
  !
  ! """
  !INOUT ----------------------------
  integer                            :: nX, nChunks, nPoly
  real,    dimension(nX)             :: x, yFit, sec_el
  integer, dimension(2,nChunks)      :: Chunks
  real,    dimension(nChunks, nPoly) :: poly
  !f2py intent(hide),depend(X)       :: nX = shape(X)
  !f2py intent(hide),depend(Chunks)  :: nChunks = shape(Chunks,1)
  !f2py intent(hide),depend(poly)    :: nPoly = shape(poly,1)
  !f2py intent(in)                   :: x, Chunks, poly
  !f2py intent(hide,out)             :: yFit
  !LOCAL ----------------------------
  integer                            :: iChunk, start, end, subSize
  real                               :: xMin,xMax, yMin, yMax
  real,    dimension(nPoly)          :: ipoly
  double precision, allocatable, dimension(:)    :: xChunk, yChunk, secChunk
  ! ---------------------------------

  yFix = 0
  
  DO iChunk=1, nChunks
        start = Chunks(1,iChunk)+1
        end   = Chunks(2,iChunk)

        ! allocate the proper array
        subSize = end-start+1
        ALLOCATE(xChunk(subSize),yChunk(subSize),secChunk(subSize))

        xChunk = x(start:end)
        secChunk = sec_el(start:end)
        ipoly = poly(iChunk,:)

        ! Evaluate the polynom ... ( evalPoly is python indexed )
        CALL evalPoly(xChunk, ipoly(1:nPoly-1), yChunk, subSize, nPoly-1) 
        yChunk = yChunk + secChunk*iPoly(nPoly)

        yFit(start:end) = yChunk 

        DEALLOCATE(xChunk,yChunk,secChunk)

     ENDDO

end subroutine evalChunkedPolysecant



!=============================================================================
subroutine evalPoly(x, poly, y, nX, nPoly )
  ! """
  ! DES: evaluate a polynom 
  ! INP: 
  !      x     = the abcsiss
  !      poly  = the polynomial coefficient
  !
  ! OUT:  y    = the polynomial values at x
  !
  ! """
  !INOUT ----------------------------
  integer                            :: nX
  double precision, dimension(nX)    :: x, y
  real, dimension(nPoly)             :: poly
  !f2py intent(hide),depend(X)       :: nX = shape(X)
  !f2py intent(hide),depend(poly)    :: nPoly = shape(poly)
  !f2py intent(in)                   :: x, poly
  !f2py intent(hide,out)             :: y
  !LOCAL ----------------------------
  integer                            :: i
  ! ---------------------------------

!  Cramer's rule:
  y = 0
  DO i=nPoly,1,-1
     y = y*x+poly(i)
  ENDDO
     
end subroutine evalPoly

