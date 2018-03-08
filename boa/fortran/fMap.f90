!module fMap

!contains

subroutine ksmooth(image,nx,ny,kernel,nk)
! """
! NAM: ksmooth (subroutine)
! DES: Smooth an image with a given kernel 
!      the kernel should have a squared support with odd dimension
! INP: image : input rank-2 array('f') with bounds (nx,ny); nx>1,ny>1
!      kernel  input rank-2 array('f') with bounds (nk,nk); nk have to be odd and >1
! OUT: image : rank-2 array('f') with bounds (nx,ny). 
! """ 
!-----------------------------------
!INOUT
 integer                              :: nx,ny, nk
 real, dimension(nx,ny)   :: image
 real, dimension(nk,nk)   :: kernel
 !f2py intent(hide),depend(image)     :: nx=shape(image,0),ny=shape(image,1)
 !f2py intent(hide),depend(kernel)    :: nk=shape(kernel,0)
 !f2py intent(in,out)                 :: image
 !f2py intent(in)                     :: kernel
!-----------------------------------
!LOCAL
 real, allocatable,dimension(:,:) :: WorkImage,subimage
 logical, allocatable,dimension(:,:)          :: mask
 integer                                      :: i,j,ng
!-----------------------------------

!---- from gsmooth

 ng = (nk-1)/2           ! range of kernel in whole pixel

 ! create an image with size of original image + 2 * the size of the kernel....

 ALLOCATE(WorkImage(nx+2*ng,ny+2*ng))

 ! ... fill it first with NaN (SQRT(-1) should provide NaN) ...
 WorkImage = -1
 WorkImage = SQRT(WorkImage)
 !  ... and put the original image in the middle
 WorkImage(ng+1:nx+ng,ng+1:ny+ng) = image


! all work is done on a subimage of the size of the kernel
 ALLOCATE(subimage(nk,nk),mask(nk,nk))


 DO i=ng+1,ng+nx       ! loop through all pixels in image
    DO j=ng+1,ng+ny
       ! Extract the corresponding subimage
       subimage = WorkImage(i-ng:i+ng,j-ng:j+ng)
       
       ! test for the presence of NaN (IEEE x.eq.y should fail for NaN even if x=y=NaN)
       ! for ifc one need to use the -mp1 option for that to work
       mask = (subimage.eq.subimage).and.(kernel.eq.kernel)
       
       ! If there is something in the subimage and the kernel AND the
       ! original pixel is not a NaN then....
       IF (any(mask).and.(image(i-ng,j-ng).eq.image(i-ng,j-ng))) THEN
          image(i-ng,j-ng) = sum( kernel*subimage, MASK=mask)/SUM(kernel, MASK=mask)
       ENDIF
    ENDDO
 ENDDO
 
 DEALLOCATE(WorkImage,subimage,mask)

end subroutine ksmooth


SUBROUTINE wcs2pix(X,Y,AXIS1,AXIS2,I,J, nInt)
! """
! DES: compute the pixel coordinates from real world coordinates
! INP:
!        X/Y  = real world coordinates
!     AXIS1/2 = the description of the 2 WCS AXIS [ NAXIS?, CRPIX?, CDELT?, CRVAL?]
! OUT:   i/j  = the corresponding pixels coordinates
!
! WAR: result is 0 indexed, remember to add +1 if used in fortran
! 
!     should be changed to use libwcs at some point
! """
!INOUT ----------------------------
  integer                      :: nInt
  real, dimension(nInt)        :: X,Y
  real, dimension(nInt)        :: I,J
  real, dimension(5)           :: AXIS1,AXIS2
  !f2py intent(hide),depend(X) :: nInt=shape(X)
  !f2py intent(in)             :: X,Y,AXIS1,AXIS2
  !f2py intent(out)            :: I,J

  ! CRVAL? = AXIS?(4)
  ! CDELT? = AXIS?(3) + unit in AXIS?(5)
  ! CRPIX? = AXIS?(2)

  ! result is 0 indexed remember to add +1 if used in fortran
  I = (X-AXIS1(4))*cos(Y*3.1415926535/180.*AXIS2(5))/AXIS1(3)+AXIS1(2)
  J = (Y-AXIS2(4))/AXIS2(3)+AXIS2(2)

END SUBROUTINE wcs2pix

SUBROUTINE wcs2phy(I,J,AXIS1,AXIS2,X,Y, nInt)
! """
! DES: compute the pixel coordinates from real world coordinates
! INP:
!        i/j  = pixel coordinates
!     AXIS1/2 = the description of the 2 WCS AXIS [ NAXIS?, CRPIX?, CDELT?, CRVAL?]
! OUT:   X/Y  = corresponding real world coordinates
!
! WAR: result is 0 indexed, remember to add +1 if used in fortran
! """
!INOUT ----------------------------
  integer                      :: nInt
  real, dimension(nInt)        :: X,Y
  real, dimension(nInt)        :: I,J
  real, dimension(5)           :: AXIS1,AXIS2
  !f2py intent(hide),depend(X) :: nInt=shape(X)
  !f2py intent(in)             :: I,J,AXIS1,AXIS2
  !f2py intent(out)            :: X,Y

  Y = AXIS2(4) + (J-AXIS2(2))*AXIS2(3)
  X = AXIS1(4) + (I-AXIS1(2))*AXIS1(3)/cos(Y*3.1415926535/180.*AXIS2(5))

END SUBROUTINE wcs2phy


SUBROUTINE HorizontalProjection(ChanList, Data, DataWeights, DataFlags, flag, &
     &AzEl, Offsets, AXIS1, AXIS2, map, weight, coverage, &
     &nChanList, nChan,nInt, dimX, dimY)

! """
! DES: Project the signal in Horizontal Coordinates according to a list of Channel
! INP:
!              ChanList = flag only this channels indexes
!                  Data = entire data array ...
!           DataWeights = ... with corresponding weights ..
!              DataFlag = ... and flags
!                 flag  = flag to select the data
!                 AzEl  = Reference Position of the array with...
!               Offsets = ... channel Offsets for the Used Channel only
!             nChanList = the lenght of ChanList
!                 nChan = nb. channels
!                  nInt = number of time steps
!            dimX dimY  = the dimension of the maps
!   map/weight/coverage = a initial version of the maps.
!               AXIS1/2 = the WCS description of the map axis
! OUT:
!   map/weight/coverage = the resulting signal, weight and coverage (weight=1) map
!         
! """
!INOUT -----------------------------
  integer                         :: nChanList, nChan, nInt, dimX, dimY
  integer, dimension(nChanList)   :: ChanList
  real,    dimension(nInt,nChan)  :: Data, DataWeights
  integer, dimension(nInt,nChan)  :: DataFlags
  integer                         :: flag
  real,    dimension(2,nInt)      :: AzEl
  real,    dimension(2,nChan)     :: Offsets
  real,    dimension(5)           :: AXIS1,AXIS2
  real,    dimension(dimX,dimY)   :: map, weight, coverage
  !f2py intent(hide),depend(ChanList) :: nChanlist = shape(ChanList)
  !f2py intent(hide),depend(Data) :: nInt = shape(Data,0), nChan = shape(Data,1)
  !f2py inten(hide),depend(map)   :: dimX = shape(map,0), dimY = shape(map,1)
  !f2py intent(in)                :: ChanList, Data, DataFlags, DataWeights, AzEl, Offsets, AXIS, AXIS2
  !f2py intent(in,out)            :: map, weight, coverage
!-----------------------------------
!LOCAL
  integer                             :: iChannel, index, nGood
  logical, dimension(nInt)            :: flag_mask
  real,    allocatable,dimension(:)   :: iAz, iEl, iData, iWeights
!  integer, allocatable,dimension(:)   :: I, J
  real, allocatable,dimension(:)      :: I, J
  real                                :: fracI,fracJ
  integer                             :: intI,intJ
!-----------------------------------

  DO iChannel=1,nChanList
     ! index in the data/flag array (remember fortran indexing)
     index = ChanList(iChannel)+1

     ! find the value we want to flag i.e. the flagIn within time_mask
     flag_mask = DataFlags(:,index).eq.flag
     nGood = COUNT(flag_mask)

     IF (nGood>1) THEN
        ALLOCATE(iAz(nGood))
        ALLOCATE(iEl(nGood))
        ALLOCATE(iData(nGood))
        ALLOCATE(iWeights(nGood))
        ALLOCATE(I(nGood))
        ALLOCATE(J(nGood))

        iAz      = PACK(AzEl(1,:),            mask=flag_mask)+Offsets(1,index)
        iEl      = PACK(AzEl(2,:),            mask=flag_mask)+Offsets(2,index)
        iData    = PACK(Data(:,index),        mask=flag_mask)
        iWeights = PACK(DataWeights(:,index), mask=flag_mask)
        
        ! Compute the corresponding pixel index
        CALL wcs2pix(iAz, iEl, AXIS1, AXIS2, I, J, nGood)

        ! wcspix if 0 indexed so...
        I = I+1
        J = J+1

        DO iGood = 1, nGood
           intI = floor(I(iGood)+0.5)
           intJ = floor(J(iGood)+0.5)
           IF ((1.le.intI).and.(intI.le.dimX).and.(1.le.intJ).and.(intJ.le.dimY)) THEN             
               map(intI,intJ)      = map(intI,intJ)      + iData(iGood)*iWeights(iGood)
               weight(intI,intJ)   = weight(intI,intJ)   + iWeights(iGood)
               coverage(intI,intJ) = coverage(intI,intJ) + 1

!            IF ((1.le.I(iGood)).and.(I(iGood).le.dimX-1).and.(1.le.J(iGood)).and.(J(iGood).le.dimY-1)) THEN
!               ! spread signal among the four neighbouring pixels
!                intI = floor(I(iGood))
!                intJ = floor(J(iGood))
!                fracI = I(iGood) - intI
!                fracJ = J(iGood) - intJ
!                map(intI,intJ)          = map(intI,intJ)          + iData(iGood)*iWeights(iGood)*(1.-fracI)*(1.-fracJ)
!                weight(intI,intJ)       = weight(intI,intJ)       + iWeights(iGood)*(1.-fracI)*(1.-fracJ)
!                coverage(intI,intJ)     = coverage(intI,intJ)     + (1.-fracI)*(1.-fracJ)
!                map(intI+1,intJ)        = map(intI+1,intJ)        + iData(iGood)*iWeights(iGood)*fracI*(1.-fracJ)
!                weight(intI+1,intJ)     = weight(intI+1,intJ)     + iWeights(iGood)*fracI*(1.-fracJ)
!                coverage(intI+1,intJ)   = coverage(intI+1,intJ)   + fracI*(1.-fracJ)
!                map(intI,intJ+1)        = map(intI,intJ+1)        + iData(iGood)*iWeights(iGood)*(1.-fracI)*fracJ
!                weight(intI,intJ+1)     = weight(intI,intJ+1)     + iWeights(iGood)*(1.-fracI)*fracJ
!                coverage(intI,intJ+1)   = coverage(intI,intJ+1)   + (1.-fracI)*fracJ
!                map(intI+1,intJ+1)      = map(intI+1,intJ+1)      + iData(iGood)*iWeights(iGood)*fracI*fracJ
!                weight(intI+1,intJ+1)   = weight(intI+1,intJ+1)   + iWeights(iGood)*fracI*fracJ
!                coverage(intI+1,intJ+1) = coverage(intI+1,intJ+1) + fracI*fracJ

           ENDIF

        ENDDO
        
        DEALLOCATE(iAz)
        DEALLOCATE(iEl)
        DEALLOCATE(iData)
        DEALLOCATE(iWeights)
        DEALLOCATE(I)
        DEALLOCATE(J)

     END IF
     
  END DO

! Allow division by 0, i.e having NaN where weight = 0  
  map = map/weight

END SUBROUTINE HorizontalProjection


SUBROUTINE EquatorialProjection(ChanList, Data, DataWeights, DataFlags, flag, &
     &RaDec, OffsetsAzEl, rotAngles, refChOffsets, &
     &AXIS1, AXIS2, neighbour, map, weight, coverage, &
     &nChanList, nChan,nInt, dimX, dimY)

! """
! DES: Project the signal in Equatorial Coordinates according to a list of Channel
! INP:
!              ChanList = flag only this channels indexes
!                  Data = entire data array ...
!           DataWeights = ... with corresponding weights ..
!              DataFlag = ... and flags
!                 flag  = flag to select the data
!                 RaDec = Reference Position of the array (EQUATORIAL!) with...
!           OffsetsAzEl = ... channel Offsets (HORIZONTAL!) for the Used Channels only
!             rotAngles = rotation angle at each time stamp
!             nChanList = the lenght of ChanList
!                 nChan = nb. channels
!                  nInt = number of time steps
!            dimX dimY  = the dimension of the maps
!   map/weight/coverage = a initial version of the maps.
!               AXIS1/2 = the WCS description of the map axis
!             neighbour = do we divide the flux to 4 neighbouring pixels?
! OUT:
!   map/weight/coverage = the resulting signal, weight and coverage (weight=1) map
!         
! """
!INOUT -----------------------------
  integer                         :: nChanList, nChan, nInt, dimX, dimY
  integer, dimension(nChanList)   :: ChanList
  real,    dimension(nInt,nChan)  :: Data, DataWeights
  integer, dimension(nInt,nChan)  :: DataFlags
  integer                         :: flag
  real,    dimension(2,nInt)      :: RaDec
  real,    dimension(2,1)         :: refChOffsets
  real,    dimension(nInt)        :: rotAngles 
  real,    dimension(2,nChan)     :: OffsetsAzEl
  real,    dimension(5)           :: AXIS1,AXIS2
  real,    dimension(dimX,dimY)   :: map, weight, coverage
  !f2py intent(hide),depend(ChanList) :: nChanlist = shape(ChanList)
  !f2py intent(hide),depend(Data) :: nInt = shape(Data,0), nChan = shape(Data,1)
  !f2py inten(hide),depend(map)   :: dimX = shape(map,0), dimY = shape(map,1)
  !f2py intent(in)                :: ChanList, Data, DataFlags, DataWeights, RaDec, OffsetsAzEl, AXIS, AXIS2, rotAngles, refChOffsets
  !f2py intent(in,out)            :: map, weight, coverage
!-----------------------------------
!LOCAL
  integer                             :: iChannel, index, nGood
  real,    dimension(2,nInt)          :: ChRaDec
  real,    dimension(2,1)             :: OffsetsRaDec, OffsetsUsedRaDec
  logical, dimension(nInt)            :: flag_mask
  real,    allocatable,dimension(:)   :: iRa, iDec, iData, iWeights
!  integer, allocatable,dimension(:)   :: I, J
  real, allocatable,dimension(:)      :: I, J
  real                                :: fracI,fracJ
  integer                             :: intI,intJ
!-----------------------------------

  DO iChannel=1,nChanList

     ! index in the data/flag array (remember fortran indexing)
     index = ChanList(iChannel)+1
     
     ! find the value we want to flag i.e. the flagIn within time_mask
     flag_mask = DataFlags(:,index).eq.flag
     nGood = COUNT(flag_mask)

     IF (nGood>1) THEN

        ALLOCATE(iRa(nGood))
        ALLOCATE(iDec(nGood))
        ALLOCATE(iData(nGood))
        ALLOCATE(iWeights(nGood))
        ALLOCATE(I(nGood))
        ALLOCATE(J(nGood))

        OffsetsAzEl(1,index) = OffsetsAzEl(1,index) - refChOffsets(1,1)
        OffsetsAzEl(2,index) = OffsetsAzEl(2,index) - refChOffsets(2,1)
        DO iInteg=1, nInt 
           !cosdec = cos(RaDec(2,iInteg)*3.14/180.)
           call RotateOffset(OffsetsAzEl(:,index),OffsetsRaDec,rotAngles(iInteg))
           !OffsetsRaDec(:,1)=OffsetsRaDec(:,1)-refChOffsets(:,1)
           !OffsetsRaDec(1,1)=-1.*OffsetsRaDec(1,1)
           chRaDec(2,iInteg)=RaDec(2,iInteg)+OffsetsRaDec(2,1)/3600.
           OffsetsRaDec(1,1) = OffsetsRaDec(1,1)/cos(chRaDec(2,iInteg)*3.1415926535/180.)
           chRaDec(1,iInteg)=RaDec(1,iInteg)+OffsetsRaDec(1,1)/3600.
        ENDDO
        
        iRa      = PACK(chRaDec(1,:),         mask=flag_mask)
        iDec     = PACK(chRaDec(2,:),         mask=flag_mask)
        iData    = PACK(Data(:,index),        mask=flag_mask)
        iWeights = PACK(DataWeights(:,index), mask=flag_mask)
        
        ! Compute the corresponding pixel index
        CALL wcs2pix(iRa, iDec, AXIS1, AXIS2, I, J, nGood)

        ! wcspix if 0 indexed so...
        I = I+1
        J = J+1

        DO iGood = 1, nGood
           IF (neighbour == 0) THEN
              intI = floor(I(iGood)+0.5)
              intJ = floor(J(iGood)+0.5)
              IF ((1.le.intI).and.(intI.le.dimX).and.(1.le.intJ).and.(intJ.le.dimY)) THEN
                 map(intI,intJ)      = map(intI,intJ)      + iData(iGood)*iWeights(iGood)
                 weight(intI,intJ)   = weight(intI,intJ)   + iWeights(iGood)
                 coverage(intI,intJ) = coverage(intI,intJ) + 1
              ENDIF
           ELSE
              intI = floor(I(iGood))
              intJ = floor(J(iGood))
              IF ((1.le.I(iGood)).and.(I(iGood).le.dimX-1).and.(1.le.J(iGood)).and.(J(iGood).le.dimY-1)) THEN
                 ! spread signal among the four neighbouring pixels
                 fracI = I(iGood) - intI
                 fracJ = J(iGood) - intJ
                 map(intI,intJ)          = map(intI,intJ)          + iData(iGood)*iWeights(iGood)*(1.-fracI)*(1.-fracJ)
                 weight(intI,intJ)       = weight(intI,intJ)       + iWeights(iGood)*(1.-fracI)*(1.-fracJ)
                 coverage(intI,intJ)     = coverage(intI,intJ)     + 1
                 map(intI+1,intJ)        = map(intI+1,intJ)        + iData(iGood)*iWeights(iGood)*fracI*(1.-fracJ)
                 weight(intI+1,intJ)     = weight(intI+1,intJ)     + iWeights(iGood)*fracI*(1.-fracJ)
                 coverage(intI+1,intJ)   = coverage(intI+1,intJ)   + 1
                 map(intI,intJ+1)        = map(intI,intJ+1)        + iData(iGood)*iWeights(iGood)*(1.-fracI)*fracJ
                 weight(intI,intJ+1)     = weight(intI,intJ+1)     + iWeights(iGood)*(1.-fracI)*fracJ
                 coverage(intI,intJ+1)   = coverage(intI,intJ+1)   + 1
                 map(intI+1,intJ+1)      = map(intI+1,intJ+1)      + iData(iGood)*iWeights(iGood)*fracI*fracJ
                 weight(intI+1,intJ+1)   = weight(intI+1,intJ+1)   + iWeights(iGood)*fracI*fracJ
                 coverage(intI+1,intJ+1) = coverage(intI+1,intJ+1) + 1
              ENDIF
           ENDIF
        ENDDO
        
        DEALLOCATE(iRa)
        DEALLOCATE(iDec)
        DEALLOCATE(iData)
        DEALLOCATE(iWeights)
        DEALLOCATE(I)
        DEALLOCATE(J)

     END IF
  END DO

! Allow division by 0, i.e having NaN where weight = 0
  map = map/weight

END SUBROUTINE EquatorialProjection


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Routines related with source model (computing, flagging, adding)
!
SUBROUTINE ComputeModel(ChanList, Data, Model,  &
     &RaDec, OffsetsAzEl, rotAngles, refChOffsets, &
     &AXIS1, AXIS2, DataArray, &
     &nChanList, nChan,nInt, dimX, dimY)
! """
! DES: Computes a data array (same size as data: nchannels x ndumps)
!      according to a source model provided as a map
! INP:
!              ChanList = do computations only for these channels
!                  Data = entire data array
!                 Model = a 2D map containing the source model
!                 RaDec = Reference Position of the array (EQUATORIAL!) with...
!           OffsetsAzEl = ... channel Offsets (HORIZONTAL!) for the Used Channels only
!             rotAngles = rotation angle at each time stamp
!          refChOffsets = offsets (HORIZONTAL!) for the reference channel
!               AXIS1/2 = the WCS description of the map axis
! Derived:
!             nChanList = the length of ChanList
!                 nChan = total number of channels
!                  nInt = number of dumps
!            dimX, dimY = dimensions of model map
! OUT:
!             DataArray = the resulting array of data
!         
! """
!INOUT -----------------------------
  integer                         :: nChanList, nChan, nInt, dimX, dimY
  integer, dimension(nChanList)   :: ChanList
  real,    dimension(nInt,nChan)  :: Data, DataArray
  real,    dimension(2,nInt)      :: RaDec
  real,    dimension(2,1)         :: refChOffsets
  real,    dimension(nInt)        :: rotAngles 
  real,    dimension(2,nChan)     :: OffsetsAzEl
  real,    dimension(5)           :: AXIS1,AXIS2
  real,    dimension(dimX,dimY)   :: Model
  !f2py intent(hide),depend(ChanList) :: nChanlist = shape(ChanList)
  !f2py intent(hide),depend(Data) :: nInt = shape(Data,0), nChan = shape(Data,1)
  !f2py inten(hide),depend(Model) :: dimX = shape(Model,0), dimY = shape(Model,1)
  !f2py intent(in)                :: ChanList, Data, Model, RaDec, OffsetsAzEl
  !f2py intent(in)                :: AXIS1, AXIS2, rotAngles, refChOffsets
  !f2py intent(out)               :: DataArray
!-----------------------------------
!LOCAL
  integer                             :: iChannel, index, nGood
  real,    dimension(2,nInt)          :: ChRaDec
  real,    dimension(2,1)             :: OffsetsRaDec, OffsetsUsedRaDec
  real,    allocatable,dimension(:)   :: iRa, iDec
  real,    allocatable,dimension(:)   :: I, J
  integer                             :: intI,intJ
  real                                :: fracI,fracJ
!-----------------------------------

  ! We will loop through the time serie for each channel
  ALLOCATE(iRa(nInt))
  ALLOCATE(iDec(nInt))
  ALLOCATE(I(nInt))
  ALLOCATE(J(nInt))
  
  DO iChannel=1,nChanList
     ! index in the data array (remember fortran indexing)
     index = ChanList(iChannel)+1
     
     OffsetsAzEl(1,index) = OffsetsAzEl(1,index) - refChOffsets(1,1)
     OffsetsAzEl(2,index) = OffsetsAzEl(2,index) - refChOffsets(2,1)
     DO iInteg=1, nInt 
        call RotateOffset(OffsetsAzEl(:,index),OffsetsRaDec,rotAngles(iInteg))
        chRaDec(2,iInteg)=RaDec(2,iInteg)+OffsetsRaDec(2,1)/3600.
        OffsetsRaDec(1,1) = OffsetsRaDec(1,1)/cos(chRaDec(2,iInteg)*3.1415926535/180.)
        chRaDec(1,iInteg)=RaDec(1,iInteg)+OffsetsRaDec(1,1)/3600.
     ENDDO
        
     iRa      = chRaDec(1,:)
     iDec     = chRaDec(2,:)
        
     ! Compute the corresponding pixel coordinates
     CALL wcs2pix(iRa, iDec, AXIS1, AXIS2, I, J, nInt)

     ! wcs2pix is 0 indexed so...
     I = I+1
     J = J+1

     DO iGood = 1, nInt
        intI = floor(I(iGood)+0.5)
        intJ = floor(J(iGood)+0.5)
        IF ((1.le.intI).and.(intI.le.dimX).and.(1.le.intJ).and.(intJ.le.dimY)) THEN
           ! get signal at closest pixel
           fracI = I(iGood) - intI
           fracJ = J(iGood) - intJ
           IF (Model(intI,intJ).eq.Model(intI,intJ)) THEN
              ! this is true only if Model(intI,intJ) is not NaN
              DataArray(iGood,index) = Model(intI,intJ)
           ELSE
              DataArray(iGood,index) = 0.
           ENDIF
        ELSE
           DataArray(iGood,index) = 0.
        ENDIF
     ENDDO
  END DO

  DEALLOCATE(iRa)
  DEALLOCATE(iDec)
  DEALLOCATE(I)
  DEALLOCATE(J)
END SUBROUTINE ComputeModel

SUBROUTINE FlagSource(ChanList, Data, Model, threshold, flag, &
     &RaDec, OffsetsAzEl, rotAngles, refChOffsets, &
     &AXIS1, AXIS2, FlagArray, &
     &nChanList, nChan,nInt, dimX, dimY)

! """
! DES: Computes an array of flags (same size as data: nchannels x ndumps)
!      according to a source model provided as a map
! INP:
!              ChanList = do computations only for these channels
!                  Data = entire data array
!                 Model = a 2D map containing the source model
!             threshold = the value (in Model pixels) above which flag is set
!                  flag = the flag value to set
!                 RaDec = Reference Position of the array (EQUATORIAL!) with...
!           OffsetsAzEl = ... channel Offsets (HORIZONTAL!) for the Used Channels only
!             rotAngles = rotation angle at each time stamp
!          refChOffsets = offsets (HORIZONTAL!) for the reference channel
!               AXIS1/2 = the WCS description of the map axis
! Derived:
!             nChanList = the length of ChanList
!                 nChan = total number of channels
!                  nInt = number of dumps
!            dimX, dimY = dimensions of model map
! OUT:
!             FlagArray = the resulting array of flags
!         
! """
!INOUT -----------------------------
  integer                         :: nChanList, nChan, nInt, dimX, dimY
  real                            :: threshold
  integer, dimension(nChanList)   :: ChanList
  real,    dimension(nInt,nChan)  :: Data
  real,    dimension(nInt,nChan)  :: FlagArray
  integer                         :: flag
  real,    dimension(2,nInt)      :: RaDec
  real,    dimension(2,1)         :: refChOffsets
  real,    dimension(nInt)        :: rotAngles 
  real,    dimension(2,nChan)     :: OffsetsAzEl
  real,    dimension(5)           :: AXIS1,AXIS2
  real,    dimension(dimX,dimY)   :: Model
  !f2py intent(hide),depend(ChanList) :: nChanlist = shape(ChanList)
  !f2py intent(hide),depend(Data) :: nInt = shape(Data,0), nChan = shape(Data,1)
  !f2py inten(hide),depend(Model) :: dimX = shape(Model,0), dimY = shape(Model,1)
  !f2py intent(in)                :: ChanList, Data, Model, RaDec, OffsetsAzEl, AXIS1, AXIS2
  !f2py intent(in)                :: rotAngles, refChOffsets, threshold
  !f2py intent(in,out)            :: FlagArray
!-----------------------------------
!LOCAL
  logical, dimension(nint)        :: mask
  integer                         :: mask_count
  integer, dimension(:), allocatable :: aFlag
!-----------------------------------

  ! compute model
  call ComputeModel(ChanList, Data, Model, RaDec, OffsetsAzEl, rotAngles, refChOffsets, &
       &AXIS1, AXIS2, FlagArray, nChanList, nChan, nInt, dimX, dimY)

  ! loop on channels to fill the Flag array
  DO iChannel=1,nChanList
     ! index in the data array (remember fortran indexing)
     index = ChanList(iChannel)+1
     mask  = FlagArray(:,index).ge.threshold
     mask_count = COUNT(mask)
     IF (mask_count>0) THEN
        allocate(aFlag(mask_count))
        aFlag = flag
        FlagArray(:,index) = UNPACK(aFlag,MASK=mask,FIELD=0)
        deallocate(aFlag)
     ELSE
        FlagArray(:,index) = 0
     ENDIF
  ENDDO
END SUBROUTINE FlagSource

SUBROUTINE AddSource(ChanList, Data, Model, &
     &RaDec, OffsetsAzEl, rotAngles, refChOffsets, &
     &AXIS1, AXIS2, factor, &
     &nChanList, nChan,nInt, dimX, dimY)

! """
! DES: add factor x simulated data to the data, where simulated
!      data have the same size as data (nchannels x ndumps) and
!      are computed according to source model
! INP:
!              ChanList = do computations only for these channels
!                  Data = entire data array
!                 Model = a 2D map containing the source model
!                 RaDec = Reference Position of the array (EQUATORIAL!) with...
!           OffsetsAzEl = ... channel Offsets (HORIZONTAL!) for the Used Channels only
!             rotAngles = rotation angle at each time stamp
!          refChOffsets = offsets (HORIZONTAL!) for the reference channel
!               AXIS1/2 = the WCS description of the map axis
!                factor = applied to the sim'ed data (e.g. -1. to *subtract* source)
! Derived:
!             nChanList = the length of ChanList
!                 nChan = total number of channels
!                  nInt = number of dumps
!            dimX, dimY = dimensions of model map
! OUT:
!               the input Data array + factor x Sim. data is returned
!         
! """
!INOUT -----------------------------
  integer                         :: nChanList, nChan, nInt, dimX, dimY
  real                            :: factor
  integer, dimension(nChanList)   :: ChanList
  real,    dimension(nInt,nChan)  :: Data
  real,    dimension(2,nInt)      :: RaDec
  real,    dimension(2,1)         :: refChOffsets
  real,    dimension(nInt)        :: rotAngles 
  real,    dimension(2,nChan)     :: OffsetsAzEl
  real,    dimension(5)           :: AXIS1,AXIS2
  real,    dimension(dimX,dimY)   :: Model
  !f2py intent(hide),depend(ChanList) :: nChanlist = shape(ChanList)
  !f2py intent(hide),depend(Data) :: nInt = shape(Data,0), nChan = shape(Data,1)
  !f2py inten(hide),depend(Model) :: dimX = shape(Model,0), dimY = shape(Model,1)
  !f2py intent(in)                :: ChanList, Model, RaDec, OffsetsAzEl, AXIS1, AXIS2
  !f2py intent(in)                :: rotAngles, refChOffsets
  !f2py intent(in,out)            :: Data
!-----------------------------------
!LOCAL
  real,    dimension(nInt,nChan)  :: SimData
!-----------------------------------
  ! compute model
  call ComputeModel(ChanList, Data, Model, RaDec, OffsetsAzEl, rotAngles, refChOffsets, &
       &AXIS1, AXIS2, SimData, nChanList, nChan, nInt, dimX, dimY)
  Data = Data + factor * SimData

END SUBROUTINE AddSource


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE RotateOffset(azel_off,radec_off,angle)

  !INOUT ----------------------------
  real                         :: angle
  real, dimension(2,1)         :: azel_off, radec_off
  !f2py intent(in)             :: azel_off, angle
  !f2py intent(out)            :: radec_off  
  !LOCAl-----------------------------
  real                         :: ang
  real, dimension(2,2)         :: rotMatrix

  ang = angle*3.14159265359/180.

  rotMatrix(1,1)=-1.*cos(ang)
  rotMatrix(1,2)=sin(ang)
  rotMatrix(2,1)=sin(ang)
  rotMatrix(2,2)=cos(ang)

  radec_off=matmul(rotMatrix,azel_off)

END SUBROUTINE RotateOffset

!----------------------------------------
subroutine mapsum(image1,weight1,coverage1,image2,weight2,coverage2,nx,ny,imageresult,coverageresult,weightresult)
!INOUT
integer                              :: nx,ny
real, dimension(nx,ny)   :: image1,image2,coverage1,coverage2,weight1,weight2
real, dimension(nx,ny)   :: imageresult,coverageresult,weightresult
!f2py intent(hide),depend(image1)     :: nx=shape(image,0),ny=shape(image,1)
!f2py intent(in)                      :: image1,image2,coverage1,coverage2,weight1,weight2
!f2py intent(out)                     :: imageresult,coverageresult,weightresult
!-----------------------------------------

imageresult = -1
imageresult = SQRT(imageresult)
coverageresult = 0
weightresult = 0

WHERE ((image1.eq.image1).and.(image2.eq.image2)) 
imageresult=image1*weight1+image2*weight2
weightresult = weight1+weight2
coverageresult = coverage1+coverage2
END WHERE

WHERE ((image1.eq.image1).and.(image2.ne.image2)) 
imageresult=image1*weight1
weightresult = weight1
coverageresult = coverage1
END WHERE

WHERE ((image1.ne.image1).and.(image2.eq.image2)) 
imageresult=image2*weight2
weightresult = weight2
coverageresult = coverage2
END WHERE

imageresult= imageresult/weightresult

END SUBROUTINE mapsum

