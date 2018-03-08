! """
! DES: module collecting f90 subroutines for sky noise subtraction:
! """ 
!=============================================================================
SUBROUTINE correlatedNoise(CNoise, data, dataFlag, chanList, boloWeight, FFCF, clip, nInt, nChan, nChanlist)
! """
! DES: compute correlated noise
! INP: 
!      data     : entire data array
!      dataFlag : corresponding flag array
!      chanList : indexes of channel where we want to compute the cn (0 indexed)
!      Weight   : weight matrix 
!      FFCF     : flat field correction factors
!      clip     : sigma clipping factor (5 or so)
!      nInt     : nb. integ.
!      nChan    : nb. channels
!  nChanList    : nb. of channel to compute
!
! OUT: CNoise   : correlated noise
!
! """
! 20070421FB fixed some errors 
!
!-----------------------------------
!INOUT
  integer                             :: nInt, nChan, nChanList
  real,    dimension(nInt,nChan)      :: data, CNoise
  integer, dimension(nInt,nChan)      :: dataFlag
  integer, dimension(nChanList)       :: chanList 
  real,    dimension(nChan,nChan)     :: boloWeight,FFCF
  real                                :: clip
  !f2py intent(hide),depend(data)     :: nInt=shape(data,0),nChan=shape(data,1)
  !f2py inten(hide), depend(chanList) :: nChanList = shape(chanList)
  !f2py intent(in)                    :: data,dataFlag,chanList,boloWeight,FFCF,clip
  !f2py intent(out)                   :: CNoise
!-----------------------------------
  integer                             :: iChan, chanIndex, count1, count2
  integer, allocatable, dimension(:)  :: nGood
  real,    allocatable, dimension(:)  :: mean, rms
  real,    allocatable,dimension(:,:) :: res,subWeight
  logical, allocatable,dimension(:,:) :: flag_mask, flag_mask_chan
!-----------------------------------

ALLOCATE(flag_mask(nInt,nChan))

! first flag all the channels false, i.e. bad
flag_mask = .false.
! ... then unflag the channels in ChanList ...
DO iChan=1, nChanList
   ! index in the data/flag array (remember fortran indexing)
   flag_mask(:,chanList(iChan) + 1) = .true.   ! is good data
END DO

! ... and finally get the unflag data only
flag_mask = flag_mask .and. (dataFlag.eq.0)

! Number of good channel vs time
ALLOCATE(nGood(nInt))
nGood= MAX(1,COUNT(flag_mask,DIM=2))

! compute mean and rms along channel dimension
! ie mean and rms of all the bolo vs time
! Variance: Numerically-stable "two-pass" formula, which offers less
!           round-off error. Page 613, Numerical Recipes in C.

ALLOCATE(mean(nInt),rms(nInt))
ALLOCATE(res(nInt,nChan))

mean = SUM(data,MASK=flag_mask,DIM=2)/nGood
res  = data-SPREAD(mean,DIM=2,NCOPIES=nChan)
!rms  = clip*SQRT(SUM(res*res,MASK=flag_mask,DIM=2)-SUM(res,MASK=flag_mask,DIM=2)**2/nGood) / (nGood-1) ! wrong
rms  = clip*SQRT((SUM(res*res,MASK=flag_mask,DIM=2)-(SUM(res,MASK=flag_mask,DIM=2)**2/nGood))/real(nGood-1))

DEALLOCATE(res,nGood)

! add sigma clipping to mask --- I am not sure this is correct yet. Seems to clip too many points.
! thus i keep output. maybe someone else can look at this in a quiet minute/hour (FB)
count1 = COUNT(flag_mask)
flag_mask = flag_mask .and. (ABS(data-SPREAD(mean,DIM=2,NCOPIES=nChan)).le.SPREAD(rms,DIM=2,NCOPIES=nChan))
count2 = COUNT(flag_mask)
write(*,'("clip ",i10," of ",i10,$)') count1-count2,count1

DEALLOCATE(mean,rms)
ALLOCATE(flag_mask_chan(nInt,nChan))
ALLOCATE(subWeight(nInt,nChan))

DO iChan=1, nChanList
!DO iChan=1, 1

   chanIndex = chanList(iChan)+1

   flag_mask_chan = flag_mask
   flag_mask_chan(:,chanIndex)=.false.

   ! FFCF are FF correction factors with respect to channel in first index
   subWeight = SPREAD(boloWeight(chanIndex,:)*FFCF(chanIndex,:),DIM=1,NCOPIES=nInt)

   write(*,'(i4,$)') chanIndex
!   print*,flag_mask_chan(100,:)
!   print*,' ' 
!   print*,subWeight(100,:)
!   print*,' ' 
!   print*,data(100,:)

   ! multiply weights into data and sum:
   CNoise(:,chanIndex) = SUM(data * subWeight,MASK=flag_mask_chan,DIM=2)

!   print*,' ' 
!   print*,CNoise(1:120,chanIndex)

   ! normalize by the weight - not by flatfield!
   ! CNoise(:,chanIndex)=CNoise(:,chanIndex)/SUM(subWeight,MASK=flag_mask_chan,DIM=2)
   subWeight = SPREAD(boloWeight(chanIndex,:),DIM=1,NCOPIES=nInt)
   CNoise(:,chanIndex)=CNoise(:,chanIndex)/SUM(subWeight,MASK=flag_mask_chan,DIM=2)

ENDDO
print*,' '
 
DEALLOCATE(subWeight,flag_mask,flag_mask_chan)  

RETURN
END SUBROUTINE correlatedNoise

!=============================================================================
SUBROUTINE cmatrix(data,flag,chanList,ccmatrix,nInt,nChan,nChanList)
! """
! DES: compute correlation matrix of an array over a sub sample of chanList
! INP: data      : entire data array
!      flag      : entire flag array
!      chanList  : the indexes used to compute the correlation
!      nInt      : nb. integ.
!      nChan     : nb. channels
!      nChanList : nb. indexes
! OUT: ccmatrix (nChan,nChan) = correlation matrix
! """
!-----------------------------------
!INOUT
  integer                             :: nInt,nChan
  real,    dimension(nInt,nChan)      :: data
  integer, dimension(nInt,nChan)      :: flag
  real,    dimension(nChan,nChan)     :: ccmatrix
  integer, dimension(nChanList)       :: chanList
  !f2py intent(hide),depend(data)     :: nInt=shape(data,0),nChan=shape(data,1)
  !f2py intent(hide),depend(chanList) :: nChanList = shape(chanList)
  !f2py intent(in)                    :: data,flag,chanList
  !f2py intent(out)                   :: ccmatrix
!-----------------------------------
!LOCAL
  integer                             :: iChan1, iChan2, indexChan1, indexChan2
  integer                             :: nGood, imax
  real                                :: sum1, sum2
  logical,dimension(nInt)             :: flag_mask
  real,allocatable,dimension(:)       :: arr1,arr2
!-----------------------------------

! initialize the matrix with NaN
ccmatrix = -1
ccmatrix = sqrt(ccmatrix)

DO iChan1=1,nChanList
   indexChan1 = chanList(iChan1)+1
   DO iChan2=iChan1,nChanList
      indexChan2 = chanList(iChan2)+1

      ! test for non flagged data in both chan1 and chan2
      flag_mask  = (flag(:,indexChan1).eq.0).and.(flag(:,indexChan2).eq.0)
      nGood = COUNT(flag_mask)

      ! we need at least 2 numbers to compute a correlation coefficient
      if(nGood>1) then
         
         ! extract common non flagged data from chan1 and chan2
         ALLOCATE(arr1(nGood),arr2(nGood))
         arr1 = PACK(data(:,indexChan1),MASK=flag_mask)
         arr2 = PACK(data(:,indexChan2),MASK=flag_mask)

         ! more accurate formula :
         ! http://mathworld.wolfram.com/CorrelationCoefficient.html
         sum1 = sum(arr1)
         sum2 = sum(arr2)
         cc = (nGood*sum(arr1*arr2)-sum1*sum2) /  &
              &sqrt((nGood*sum(arr1*arr1)-sum1**2)*(nGood*sum(arr2*arr2)-sum2**2))
         
         ccmatrix(indexChan1,indexChan2) = cc
         ccmatrix(indexChan2,indexChan1) = cc
         DEALLOCATE(arr1,arr2)
      ENDIF
   END DO
END DO

END SUBROUTINE cmatrix

!=============================================================================
SUBROUTINE wmatrix(ccmatrix, chanSep, chanList, minCorr, a, b, core, beta, wwmatrix, nChan, nChanList)
! """
! DES: compute the weight matrix, a non linear rescaling of the correlation matrix
! INP: ccmatrix  : the correlation matrix
!      chanSep   : the channel separation in arcsec
!      minCorr   : cut of value for the correlation
!      a,b       : coefficient for the rescaling ( wmatrix_nm = ( cmatrix_nm - a * min_m( cmatrix_nm ) )**b )
!      core,beta : additional coefficient for rescaling ( weight_nm  = weight_nm * 1.0 / ( 1 + ( dist_nm / core(in arcmin) )**beta ) )
!      chanList  : the indexes used to compute the correlation
!      nChan     : nb. channels
!      nChanList : nb. indexes
! OUT: wwmatrix (nChan,nChan) = weight matrix
! """
!-----------------------------------
!INOUT
  integer                             :: nChan, nChanList
  real,    dimension(nChan,nChan)     :: ccmatrix, wwmatrix, chanSep
  integer, dimension(nChanList)       :: chanList
  real                                :: minCorr,a,b, core, beta
  !f2py intent(hide),depend(ccmatrix) :: nChan=shape(ccmatrix,0)
  !f2py intent(hide),depend(chanList) :: nChanList = shape(chanList)
  !f2py intent(in)                    :: ccmatrix, chanSep, chanList, minCorr, a, b
  !f2py intent(out)                   :: wwmatrix
!-----------------------------------
!LOCAL
  integer                             :: iChan1, iChan2, indexChan1, indexChan2
  integer                             :: mask_count
  real, dimension(:), allocatable     :: aMinCorr
  logical, dimension(nChan,nChan)     :: flag_mask
  logical, dimension(nChan)           :: min_flag_mask
!-----------------------------------
! first flag all the channels...
flag_mask = .false.
! ... then unflag the channels in ChanList ...

DO iChan1=1, nChanList
   indexChan1 = chanList(iChan1)+1
   DO iChan2=iChan1,nChanList
      indexChan2 = chanList(iChan2)+1
      ! index in the data/flag array (remember fortran indexing)
      flag_mask(indexChan1,indexChan2) = .true.
      flag_mask(indexChan2,indexChan1) = .true.
   END DO
END DO

wwmatrix(:,1)  = MINVAL(ccmatrix,DIM=2,MASK=flag_mask)  ! will return largest possible number is MASK=.false.
min_flag_mask  = wwmatrix(:,1).lt.minCorr               ! check where minimum is < minCorr
mask_count     = COUNT(min_flag_mask)

IF (mask_count>1) THEN                                  ! limit minimum to minCorr
   allocate(aMinCorr(mask_count))
   aMinCorr = minCorr
   wwmatrix(:,1) = UNPACK(aMinCorr,MASK=min_flag_mask,FIELD=wwmatrix(:,1))
   deallocate(aMinCorr)
ENDIF

wwmatrix = (ccmatrix-SPREAD(a*wwmatrix(:,1),DIM=2,NCOPIES=nChan))**b   ! spread the weights

IF (core > 0) THEN 
   wwmatrix = wwmatrix * 1.0 / ( 1. +( ChanSep/60. / core )**beta )       ! taper with beta profile
ENDIF
IF (core < 0) THEN 
   wwmatrix = wwmatrix * (1.0 - (1.0 / ( 1. +( ChanSep/60. / (-core) )**beta )))  ! 1 - beta profile
ENDIF

END SUBROUTINE

!=============================================================================
subroutine correlationFit(correlateFrom, correlateTo, flag, chanList, minMaxSlope, slope, intercep, FFCF, nInt,nChan,nChanList)
! """
! DES: order-1 (linear) fit to data vs. data2
!      Produces fit coefficients and fit array matching the data (including
!      the flagged!).
! INP: correlateFrom  : data or skynoise array
!      correlateTo    : data
!      flag           : corresponding flag array
!    chanList         : indexes of channel where we want to compute the fit (0 indexed)
!      nInt           : nb. integ.
!      nChan          : nb. channels
!  nChanList          : nb. of channel to compute
!
! OUT: slope          : slopes of linear fit   (nChanList x nChanList)
!     intercep        : intercep of linear fit (nChanList x nChanList)
!      FFCF           : Flat field correction factor
! """

!INOUT ----------------------------
  integer                                  :: nInt, nChan, nChanList
  real,    dimension(nInt,nChan)           :: correlateFrom, correlateTo
  integer, dimension(nInt,nChan)           :: flag
  real,    dimension(nChanList,nChanList)  :: slope, intercep
  integer, dimension(nChanList)            :: chanList 
  real,    dimension(nChan,nChan)          :: FFCF
  real,    dimension(2)                    :: minMaxSlope
  !f2py intent(hide),depend(correlateFrom) :: nInt=shape(correlateFrom,0),nChan=shape(correlateFrom,1)
  !f2py intent(hide),depend(chanList)      :: nChanList=shape(chanList)
  !f2py intent(in)                         :: data, flag, chanList, minMaxSlope
  !f2py intent(out,hide)                   :: slope, intercep, FFCF
!LOCAL ---------------------------- 
  integer                                  :: iChannel, index, indexRef, nGood
  real                                     :: chisq, iSlope, iIntercep
  logical, dimension(nInt)                 :: flag_mask
  real,allocatable,dimension(:)            :: x,y
! ---------------------------------

! Initialize the Flat Field
FFCF = 0

! loop through the requested channels
DO iRefChannel=1, nChanList
   ! index in the data/flag array (remember fortran indexing)
   indexRef = chanList(iRefChannel) + 1

   DO iChannel=iRefChannel, nChanList
      ! index in the data/flag array (remember fortran indexing)
      index = chanList(iChannel) + 1
      
      ! find the non flagged value
      flag_mask = (flag(:,indexRef).eq.0).and.(flag(:,index).eq.0)
      nGood = COUNT(flag_mask)
      
      if(nGood>1) then
     
         ! Extract the non flagged data
         ALLOCATE( x(nGood), y(nGood) )

         x = PACK(correlateFrom(:,indexRef) ,MASK=flag_mask) 
         y = PACK(  correlateTo(:,index),    MASK=flag_mask)
         
         call fit(x,y,iIntercep,iSlope,nGood)

         intercep(iRefChannel, iChannel) = iIntercep
         if(iChannel.ne.iRefChannel) intercep(iChannel, iRefChannel) = -1.0*iIntercep/iSlope
!         if(iChannel.ne.iRefChannel) intercep(iChannel, iRefChannel) = iIntercep

         slope(iRefChannel, iChannel)    = iSlope
         if(iChannel.ne.iRefChannel) slope(iChannel, iRefChannel)    = 1./iSlope
!         if(iChannel.ne.iRefChannel) slope(iChannel, iRefChannel)    = iSlope

         FFCF(indexRef, index)           = iSlope
         if(index.ne.indexRef) FFCF(index, indexRef)           = 1./iSlope

         DEALLOCATE(x,y)
         
      END IF
   END DO
END DO

! Saturate the FFCF and inverse

WHERE( FFCF.ge.maxval(minMaxSlope) ) FFCF = maxval(minMaxSlope)
WHERE( FFCF.le.minval(minMaxSlope) ) FFCF = minval(minMaxSlope)

! Array filled with NaN for non computed channels... 
! Diagnostic tools for further routines which should not use them

FFCF = 1./FFCF

! FFCF is NOT normalized at this point...

END SUBROUTINE correlationFit


SUBROUTINE printcipalComponentAnalysis(data, nComp, dataBack, nx,ny)
!
! Do not use yet
!

integer                   :: nx, ny, nComp
real, dimension(nx,ny)    :: data, dataBack

!f2py intent(hide),depend(data)  :: nx=shape(data,0),ny=shape(data,1)
!f2py intent(in)                 :: data,nComp
!f2py intent(out)                :: dataBack

real, dimension(ny)       :: mean
real, dimension(nx,ny)    :: adjustedData
real, dimension(nx,nComp) :: finalData

real, dimension(nComp,ny)           :: desiredFeature
double precision, dimension(ny,ny)  :: eigenVectors
double precision, dimension(ny)     :: eigenValues
double precision, dimension(3*ny-1) :: work
integer                             :: info

mean = SUM(data, DIM=1)/nX
adjustedData = data-SPREAD(mean,DIM=1,NCOPIES=nX)

! fill eigenvectors matrix with covariance matrix
eigenVectors = MATMUL(transpose(adjustedData),adjustedData)/(nY+1)

CALL DSYEV('V','U',ny,eigenVectors,ny,eigenvalues,work,3*ny-1,info)

print *,info
desiredFeature = eigenVectors(1:nComp,:)

finalData = MATMUL(adjustedData,transpose(desiredFeature))
dataBack  = MATMUL(finalData,desiredFeature)+SPREAD(mean,DIM=1,NCOPIES=nX)

END SUBROUTINE printcipalComponentAnalysis
