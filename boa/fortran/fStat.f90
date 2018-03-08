! Copyright (C) 2002-2007
! Max-Planck-Institut fuer Radioastronomie Bonn
! Argelander Institut fuer Astronomie
! Astronomisches Institut der Ruhr-Universitaet Bochum
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
!$$$$

!=============================================================================
subroutine minmax(array,nx,extrema)
! """
! NAM: minmax
! DES: return maximum and minimum of a 1D array
! """ 
!INOUT ----------------------------
  real, dimension(nx)              :: array
  integer                          :: nx
  real, dimension(2)               :: extrema
  !f2py intent(hide),depend(array) :: nx = shape(array)
  !f2py intent(in)                 :: array
  !f2py intent(out)                :: extrema
!-----------------------------------

  extrema(1) = minval(array)
  extrema(2) = maxval(array)

END SUBROUTINE minmax


!=============================================================================
SUBROUTINE arrayStat_s(data,flag,index_s,mean_s,med_s,sdev_s,mdev_s,nx,ny,ns)
! """
! NAM: arrayStat  (sub)
! DES: compute mean/sdev and median/mdev of a given array sliced by index_s
!      with a given mask (flag=0), along the first dimension 
! INP: data        = entire data array
!      flag        = entire flag array
!      index_s     = the indexes of the subscans
! OUT: mean = 2D float array with the (ny,nz) dimension
!      med  = 2D float array
!      sdev = 2D float array
!      mdev = 2D float array
! """
!-----------------------------------
!INOUT
  integer :: nx,ny,ns
  real,    dimension(nx,ny)          :: data
  integer, dimension(nx,ny)          :: flag
  integer, dimension(2,ns)           :: index_s
  real,    dimension(ny,ns)          :: mean_s, med_s, sdev_s, mdev_s
  !f2py intent(in)                   :: data, flag, index_s
  !f2py intent(hide),depend(flag)    :: nx=shape(flag,0), ny=shape(flag,1)
  !f2py intent(hide),depend(index_s) :: ns=shape(index_s,0)
  !f2py intent(out)                  :: mean_s, med_s, sdev_s, mdev_s
!------------------------------------
!LOCAL
  integer                            :: index, nItem, lo, hi
  real,allocatable,dimension(:,:)    :: tData 
  integer,allocatable,dimension(:,:) :: tFlag
  real, dimension(ny)                :: mean_i, med_i, sdev_i, mdev_i
!------------------------------------

DO index=1,ns        ! loop through the subscans

   lo = index_s(1,index)+1
   hi = index_s(2,index)

   nItem = hi-lo+1

   ALLOCATE(tData(nItem,ny))
   ALLOCATE(tFlag(nItem,ny))

   tData = data(lo:hi,:)
   tFlag = flag(lo:hi,:)

   call arrayStat(tData,tFlag,mean_i,med_i,sdev_i,mdev_i, nItem,ny)
   mean_s(:,index) = mean_i
   med_s(:,index)  = med_i
   sdev_s(:,index) = sdev_i
   mdev_s(:,index) = mdev_i

   DEALLOCATE(tData,tFlag)

END DO

END SUBROUTINE arrayStat_s

!=============================================================================
SUBROUTINE arrayStat(data,flag,mean,med,sdev,mdev,nx,ny)
! """
! NAM: arrayStat  (sub)
! DES: compute mean/sdev and median/median dev of a given 2D array (data) 
!      with a given mask (flag=0), along the first dimension 
! INP: data        = entire data array
!      flag        = entire flag array
! OUT: mean = float array with the second dimension of data
!      med  = float array
!      sdev = float array
!      mdev = float array
! """
!-----------------------------------
!INOUT
  integer :: nx,ny
  real,    dimension(nx,ny)       :: data
  integer, dimension(nx,ny)       :: flag
  real,    dimension(ny)          :: mean,med,sdev,mdev
  !f2py intent(hide),depend(flag) :: nx=shape(flag,0), ny=shape(flag,1)
  !f2py intent(out)               :: mean, med, sdev, mdev
  !f2py intent(in)                :: data,flag
!-----------------------------------
!LOCAL
  integer                         :: index,mask_count
  real                            :: mean_i, med_i, sdev_i, mdev_i
  logical,dimension(nx)           :: flag_mask
  real,allocatable,dimension(:)   :: array
!-----------------------------------

DO index=1,ny        ! loop through channels

   ! take only the 0 flagged data
   flag_mask  = flag(:,index).eq.0
   mask_count = COUNT(flag_mask)
   

   if(mask_count>0) then
      ALLOCATE(array(mask_count))
      array = PACK(data(:,index),MASK=flag_mask)
      CALL f_stat(array,mean_i,sdev_i,med_i,mdev_i,mask_count)
      DEALLOCATE(array)

      mean(index) = mean_i
      sdev(index) = sdev_i
      med(index)  = med_i
      mdev(index) = mdev_i
   else
      mean(index) = -1
      sdev(index) = -1
      med(index)  = -1
      mdev(index) = -1
   end if
      

END DO

END SUBROUTINE arrayStat

!=============================================================================
SUBROUTINE arrayMedian(data,flag,med,nx,ny)
! """
! NAM: arrayMedian  (sub)
! DES: compute median of a given 2D array (data) with a given mask (flag=0),
!      along the first dimension 
! INP: data        = entire data array
!      flag        = entire flag array
! OUT: med = float array with the second dimension of data
! """
!-----------------------------------
!INOUT
  integer :: nx,ny
  real,    dimension(nx,ny)       :: data
  integer, dimension(nx,ny)       :: flag
  real,    dimension(ny)          :: med
  !f2py intent(hide),depend(flag) :: nx=shape(flag,0), ny=shape(flag,1)
  !f2py intent(out)               :: med
  !f2py intent(in)                :: data,flag
!-----------------------------------
!LOCAL
  integer                         :: index,mask_count
  real                            :: med_i
  logical,dimension(nx)           :: flag_mask
  real,allocatable,dimension(:)   :: array
!-----------------------------------

DO index=1,ny        ! loop through channels

   ! take only the 0 flagged data
   flag_mask  = flag(:,index).eq.0
   mask_count = COUNT(flag_mask)
   
   if(mask_count>0) then
      ALLOCATE(array(mask_count))
      array = PACK(data(:,index),MASK=flag_mask)
      CALL f_median(array,med_i,mask_count)
      DEALLOCATE(array)

      med(index)  = med_i
   else
      med(index)  = 0
   end if

END DO

END SUBROUTINE arrayMedian

!=============================================================================
SUBROUTINE arrayMean(data,flag,mean,nx,ny)
! """
! NAM: arrayMean  (sub)
! DES: compute mean of a given 2D array (data) with a given mask (flag=0),
!      along the first dimension 
! INP: data        = entire data array
!      flag        = entire flag array
! OUT: mean = float array with the second dimension of data
! """
!-----------------------------------
!INOUT
  integer :: nx,ny
  real,    dimension(nx,ny)       :: data
  integer, dimension(nx,ny)       :: flag
  real,    dimension(ny)          :: mean
  !f2py intent(hide),depend(flag) :: nx=shape(flag,0), ny=shape(flag,1)
  !f2py intent(out)               :: mean
  !f2py intent(in)                :: data,flag
!-----------------------------------
!LOCAL
  integer                         :: index,mask_count
  real                            :: mean_i
  logical,dimension(nx)           :: flag_mask
  real,allocatable,dimension(:)   :: array
!-----------------------------------

DO index=1,ny        ! loop through channels

   ! take only the 0 flagged data
   flag_mask  = flag(:,index).eq.0
   mask_count = COUNT(flag_mask)
   
   if(mask_count>0) then
      ALLOCATE(array(mask_count))
      array = PACK(data(:,index),MASK=flag_mask)
      CALL f_mean(array,mean_i,mask_count)
      DEALLOCATE(array)

      mean(index)  = mean_i
   else
      mean(index)  = -1
   end if

END DO

END SUBROUTINE arrayMean

!=============================================================================
SUBROUTINE f_stat(data,mean,sdev,med,mdev,i)
! """
! NAM: f_stat  (sub)
! DES: compute mean median, sdev and median dev of a given array
! INP: data         = the data array
! OUT: mean = f scalar
!      med  = f scalar
!      sdev  = f scalar
!      mdef  = f scalar
! """
!-----------------------------------
!INOUT
  integer                          :: i
  real,    dimension(i)            :: data
  real                             :: mean,sdev,med,mdev
  !f2py intent(hide),depend(data)  :: i=shape(data)
  !f2py intent(out)                :: mean,med,sdev,mdev
  !f2py intent(in)                 :: data
!-----------------------------------
!LOCAL
  real,     dimension(i)         :: res,sortedData
!-----------------------------------

  ! need at least two point
  if (i <= 1) then
     mean = -1
     med  = -1
     sdev = -1
     mdef = -1
     return
  end if

  CALL f_mean(data,mean,i)
  CALL f_median(data,med,i)
  CALL f_rms(data,mean,sdev,i)
  CALL f_rms(data,med,mdev,i)

END SUBROUTINE f_stat

!=============================================================================
SUBROUTINE f_mean(data,mean,i)
! """
! NAM: mean (sub)
! DES: compute mean of a given array
! INP: data = the data array
! OUT: mean = f scalar
! """
!-----------------------------------
!INOUT
  integer                          :: i                     
  real,    dimension(i)            :: data
  real                             :: mean
  !f2py intent(hide),depend(data)  :: i=shape(data)
  !f2py intent(in)                 :: data
  !f2py intent(out)                :: mean
!-----------------------------------

  mean = SUM(data(:)) / i

END SUBROUTINE f_mean

!=============================================================================
SUBROUTINE f_rms(data,mean,sdev,i)
! """
! NAM: rms  (sub)
! DES: compute mean and sdev of a given array
! INP: data = the data array
!    : mean = the mean of the data scalar
! OUT: sdev = f scalar
! """
!-----------------------------------
!INOUT
  integer :: i                     
  real,    dimension(i)            :: data
  real                             :: mean,sdev
  !f2py intent(hide),depend(data)  :: i=shape(data)
  !f2py intent(in)                 :: data, mean
  !f2py intent(out)                :: sdev
!-----------------------------------
!LOCAL
  real,     dimension(i)         :: res
!-----------------------------------
  res(:) = data(:) - mean ! these ARE the residuals
  
  ! Variance: Numerically-stable "two-pass" formula, which offers less
  !           round-off error. Page 613, Numerical Recipes in C.
  !sdev = sqrt( ( sum( res(:)*res(:) ) - sum( res(:) )**2 / i) / (i-1) )
  sdev = sqrt( sum( res(:)*res(:) ) / i) 

END SUBROUTINE f_rms

!=============================================================================
SUBROUTINE f_median(data,med,i)
! """
! NAM: median  (sub)
! DES: compute median of a given array
! INP: data         = the data array
! OUT: med = f scalar
! """
!-----------------------------------
!INOUT
  integer                          :: i
  real,    dimension(i)            :: data
  real                             :: med
  !f2py intent(hide),depend(data)  :: i=shape(data)
  !f2py intent(in)                 :: data
  !f2py intent(out)                :: med
!-----------------------------------
!LOCAL
  real,     dimension(i)         :: sortedData
!-----------------------------------
  ! sort the data with the shell method
  call shell(data,sortedData,i)

  if (mod(i,2).eq.0) then
     med = (sortedData(i/2)+sortedData(i/2+1))/2
  else
     med = sortedData((i+1)/2) 
  end if

END SUBROUTINE f_median

!=============================================================================
SUBROUTINE shell(datain,dataout,n)
! """
! NAM: shell  (sub)
! DES: sort an input array with the shell method
!      from N.R. 8.2
! INP: datain   = the input array
! OUT: dataout  = the sorted array
! """
!-----------------------------------
!INOUT
  integer :: n                     ! these must appear in the call parameter list but are hidden
  real,    dimension(n)            :: datain,dataout
  !f2py intent(hide),depend(datain)  :: n=shape(datain)
  !f2py intent(out)                :: dataout
  !f2py intent(in)                 :: datain
!-----------------------------------
!LOCAL
  integer                        :: i,j,inc
  real                           :: v
!-----------------------------------

  dataout = datain

  ! Determine the starting increment
  inc=1
  do while (inc <= n)
     inc=3*inc+1
  end do
  
  do while (inc > 1)
     inc=inc/3
     
     do i=inc+1 ,n
        v=dataout(i)
        j=i
        do
           if (dataout(j-inc) <= v) exit
           dataout(j)=dataout(j-inc)
           j=j-inc
           if (j <= inc) exit
        end do
        dataout(j)=v
     end do
  end do
  
END SUBROUTINE shell

!=============================================================================
SUBROUTINE quickSort(datain,dataout,n)
! """
! NAM: shell  (sub)
! DES: sort an input array with the shell method
!      N.R. 8.1
! INP: datain   = the input array
! OUT: dataout  = the sorted array
! """
!-----------------------------------
!INOUT
  real,    dimension(n)            :: datain,dataout
  integer :: n                     ! these must appear in the call parameter list but are hidden
  !f2py intent(hide),depend(datain)  :: n=shape(datain)
  !f2py intent(out)                :: dataout
  !f2py intent(in)                 :: datain
!-----------------------------------
!LOCAL
  integer                        :: i, j, k, l, r, s, stackl(15), stackr(15), ww
  real                           :: w, x
!-----------------------------------

! modified from http://www.it.uu.se/edu/course/homepage/algpar1/ht02/Quicksort_examples/quicksort.f90

  dataout = datain
  
  s = 1
  stackl(1) = 1
  stackr(1) = n
  
  !     KEEP TAKING THE TOP REQUEST FROM THE STACK UNTIL S = 0.

  DO WHILE (s.NE.0)
     
     l = stackl(s)
     r = stackr(s)
     s = s - 1
     
     !     KEEP SPLITTING DATAOUT(L), ... , DATAOUT(R) UNTIL L >= R.

     DO WHILE (l.LT.r) 
        i = l
        j = r
        k = (l+r) / 2
        x = dataout(k)
        
        !     REPEAT UNTIL I > J.
        
        DO
           DO
              IF (dataout(i).LT.x) THEN                ! Search from lower end
                 i = i + 1
                 CYCLE
              ELSE
                 EXIT
              END IF
           END DO
           
           DO
              IF (x.LT.dataout(j)) THEN                ! Search from upper end
                 j = j - 1
                 CYCLE
              ELSE
                 EXIT
              END IF
           END DO
           
           IF (i.LE.j) THEN                     ! Swap positions i & j
              w = dataout(i)
              dataout(i) = dataout(j)
              dataout(j) = w
              i = i + 1
              j = j - 1
              IF (i.GT.j) EXIT
           ELSE
              EXIT
           END IF
        END DO
        
        IF (j-l.GE.r-i) THEN
           IF (l.LT.j) THEN
              s = s + 1
              stackl(s) = l
              stackr(s) = j
           END IF
           l = i
        ELSE
           IF (i.LT.r) THEN
              s = s + 1
              stackl(s) = i
              stackr(s) = r
           END IF
           r = j
        END IF
        
     END DO
  END DO
  RETURN
  
END SUBROUTINE quickSort

!=============================================================================
SUBROUTINE compress(data,n,datacomp,sdev,i)
! """
! NAM: compress  (sub)
! DES: compress a 1D array over n data
! INP: data = the data array
!    :    n = the number of data point to compress
! OUT: datacomp = the resulting array
!          sdev = the corresponding sdev
! """
!-----------------------------------
!INOUT
  integer :: i                     
  real,    dimension(i)            :: data
  real,    dimension(i/n)          :: datacomp, sdev
  !f2py intent(hide),depend(data)  :: i=shape(data)
  !f2py intent(in)                 :: data, n
  !f2py intent(out)                :: datacomp,sdev
  
  integer :: j
  real    :: iMean, iSdev
  

DO j=0,i/n

   CALL f_mean(data(j*n+1:(j+1)*n+1),iMean,n)
   datacomp(j+1) = iMean
   CALL f_rms(data(j*n+1:(j+1)*n+1),iMean,iSdev,n)
   sdev(j+1)       = iSdev
     
ENDDO

END SUBROUTINE compress


!=============================================================================
SUBROUTINE compressWithFlag(data,flag,n,datacomp,flagcomp,sdev,i)
! """
! NAM: compress  (sub)
! DES: compress a 1D array with flag over n data
! INP: data = the data array
!      flag = the flag array
!    :    n = the number of data point to compress
! OUT: datacomp = resulting data array
!      flagcomp = resulting flag array
!          sdev = corresponding sdev
! """
!-----------------------------------
!INOUT
  integer                          :: i                     
  real,    dimension(i)            :: data
  integer, dimension(i)            :: flag
  real,    dimension(i/n)          :: datacomp, sdev
  integer, dimension(i/n)          :: flagcomp
  !f2py intent(hide),depend(data)  :: i=shape(data)
  !f2py intent(in)                 :: data, flag, n
  !f2py intent(out)                :: datacomp,flagcomp,sdev
  
  integer                           :: j
  real                              :: iMean, iSdev
  real,    allocatable,dimension(:) :: tData
  logical, dimension(n)             :: flag_mask

DO j=0,i/n

   ! take only the 0 flagged data
   flag_mask  = flag(j*n+1:(j+1)*n+1).eq.0
   mask_count = COUNT(flag_mask)
 
   ! we need at least 2 unflagged points
   if(mask_count>1) then
      ALLOCATE(tData(mask_count))

      tData = PACK(data(j*n+1:(j+1)*n+1),MASK=flag_mask)
      CALL f_mean(tData, iMean, mask_count)
      CALL f_rms(tData,iMean,iSdev,mask_count)

      DEALLOCATE(tData)
      
      datacomp(j+1) = iMean
      sdev(j+1)     = iSdev
      flagcomp(j+1) = 0
   else
      datacomp(j+1) = 0
      sdev(j+1)     = 0
      flagcomp(j+1) = -1
   end if
     
ENDDO

END SUBROUTINE compressWithFlag

!=============================================================================
SUBROUTINE arrayCompress(data,flag,n,datacomp,flagcomp,sdev,nx,ny)
! """
! NAM: arrayCompress  (sub)
! DES: compress a 2D array with flag over n data in the first dimension
! INP: data = entire data array
!      flag = entire flag array
!    :    n = the number of data point to compress
! OUT: datacomp = resulting data array
!      flagcomp = resulting flag array
!          sdev = the corresponding sdev
! """
!-----------------------------------
!INOUT
  integer                          :: nx,ny
  real,    dimension(nx,ny)        :: data
  integer, dimension(nx,ny)        :: flag
  real,    dimension(nx/n,ny)      :: datacomp, sdev
  integer, dimension(nx/n,ny)      :: flagcomp
  !f2py intent(hide),depend(flag)  :: nx=shape(flag,0), ny=shape(flag,1)
  !f2py intent(in)                 :: data,flag, n
  !f2py intent(out)                :: datacomp,flagcomp,sdev
  
  integer :: index
  real,    dimension(n)            :: idatacomp, isdev
  integer, dimension(n)            :: iflagcomp
  

DO index=1,ny        ! loop through channels

   CALL compressWithFlag(data(:,index),flag(:,index),n,idatacomp,iflagcomp,isdev,nx)
   datacomp(:,index) = idatacomp
   flagcomp(:,index) = iflagcomp
   sdev(:,index)     = isdev

ENDDO

END SUBROUTINE arrayCompress

!=============================================================================
SUBROUTINE histogram(datain,mini,step,nBin,distrib,n)
! """
! NAM: histogram  (sub)
! DES: return distribution function
! INP: datain   = the input array
!      mini     = lower value
!      step     = bin size
!      nBin     = number of bins
! OUT: distrib  = the histogram values
! """
!-----------------------------------
!INOUT
  integer :: n,nBin
  real    :: mini,step
  real,    dimension(n)            :: datain
  real,    dimension(nBin)         :: distrib
  !f2py intent(hide),depend(datain)  :: n=shape(datain)
  !f2py intent(out)                :: distrib
  !f2py intent(in)                 :: datain,mini,step,nBin
!-----------------------------------
!LOCAL
  integer                     :: i
  real                        :: val1,val2
  logical, dimension(n)       :: mask
!-----------------------------------

  do i=1,nBin
     val1 = mini + (i-1)*step
     val2 = mini + i*step
     mask  = (datain.ge.val1).and.(datain.lt.val2)
     distrib(i) = COUNT(mask)
  end do
  
END SUBROUTINE histogram

!=============================================================================
SUBROUTINE clipping(datain,kappa,dataout,nOut,n)
! """
! NAM: clipping  (sub)
! DES: perform kappa-sigma clipping
! INP: datain   = the input array
!      kappa    = number of sigma for clipping
! OUT: dataout  = the clipped array
! """
!-----------------------------------
!INOUT
  integer                          :: n,nOut
  real                             :: kappa
  real, dimension(n)               :: datain,dataout
  !f2py intent(hide),depend(datain):: n=shape(datain)
  !f2py intent(out)                :: dataout,nOut
  !f2py intent(in)                 :: datain,kappa
!-----------------------------------
!LOCAL
  real                        :: mean,rms
  logical, dimension(n)       :: mask
!-----------------------------------
  CALL f_mean(datain,mean,n)
  CALL f_rms(datain,mean,rms,n)
  mask  = (datain.ge.mean - kappa*rms).and.(datain.le.mean + kappa*rms)
  nOut = COUNT(mask)
  dataout = PACK(datain,MASK=mask)

END SUBROUTINE clipping

!=============================================================================
SUBROUTINE meanDistribution(datain,nbX,nbY,nbTotal,cell,allMean,nbOut,nx,ny)
! """
! NAM: meanDistribution (sub)
! DES: return distribution of means (by cell) in a 2D array
! INP: datain   = the input array
!      kappa    = number of sigma for clipping
! OUT: allMean  = a 1D array of mean values
! """
!-----------------------------------
!INOUT
  integer                            :: nbX,nbY,nbTotal,cell,nx,ny
  real, dimension(nx,ny)             :: datain
  real, dimension(nbTotal)           :: allMean
  !f2py intent(hide),depend(datain)  :: nx=shape(datain,0), ny=shape(datain,1)
  !f2py intent(in)                   :: datain,nbX,nbY,nbTotal,cell
  !f2py intent(out)                  :: allMean,nbOut
!-----------------------------------
!LOCAL
  real                          :: mean
  integer                       :: i,j,k,nOut,nbOut
  real, dimension(cell*cell)    :: subimage,subclean
  logical, dimension(cell*cell) :: mask
  logical, dimension(nbTotal)   :: allmask
!-----------------------------------
  do i=1,nbX
     do j=1,nbY
        do k=0,cell-1
           subimage(k*cell+1:k*cell+cell) = datain(i+1:i+cell,j+k)
        end do
        mask = (subimage.eq.subimage)  ! true if not NaN
        nOut = COUNT(mask)
        if (nOut>0) then
           subclean = PACK(subimage,MASK=mask)
           CALL f_mean(subclean,mean,nOut)
           allMean(j+nbY*(i-1)) = mean
           allmask(j+nbY*(i-1)) = .true.
        else
           allmask(j+nbY*(i-1)) = .false.
        end if
     end do
  end do
  allMean = PACK(allMean,MASK=allmask)
  nbOut = COUNT(allmask)

END SUBROUTINE meanDistribution

!=============================================================================
SUBROUTINE slidingRms(datain,flags,nbInteg,rms,nInt)
! """
! NAM: slidingRms
! DES: return RMS computed over a sliding window
! INP: datain   = the input array
!      flags    = associated flags
!      nbInteg  = size of the sliding window
! OUT: rms      = array of rms, same size as datain
! """
!-----------------------------------
!INOUT
  integer                            :: nInt,nbInteg
  real, dimension(nInt)              :: datain,rms
  integer, dimension(nInt)           :: flags
  !f2py intent(hide),depend(datain)  :: nInt=shape(datain)
  !f2py intent(in)                   :: datain,flags,nbInteg
  !f2py intent(out)                  :: rms
!-----------------------------------
!LOCAL
  integer                          :: half,i,mask_count,first,prev,next
  real,allocatable,dimension(:)    :: tmpData
  logical,allocatable,dimension(:) :: flag_mask
  logical, dimension(nInt)         :: zero_mask
  real                             :: sdev
!-----------------------------------
  half = nbInteg/2
  ! first datapoints until half
  DO i=1,half
     ALLOCATE(flag_mask(i+half))
     flag_mask = flags(1:i+half).eq.0
     mask_count = COUNT(flag_mask)
     ! need at least half/2 points to compute an rms
     IF(mask_count>half/2) then
        ALLOCATE(tmpData(mask_count))
        tmpData = PACK(datain(1:i+half),MASK=flag_mask)
        CALL f_rms(tmpData,0.,sdev,mask_count)
        rms(i) = sdev
        DEALLOCATE(tmpData)
     ELSE
        rms(i) = 0.
     END IF
     DEALLOCATE(flag_mask)
  END DO

  ! datapoints between half and end-half
  ALLOCATE(flag_mask(2*half+1))
  DO i=half+1,nInt-half
     flag_mask = flags(i-half:i+half).eq.0
     mask_count = COUNT(flag_mask)
     IF(mask_count>half/2) then
        ALLOCATE(tmpData(mask_count))
        tmpData = PACK(datain(i-half:i+half),MASK=flag_mask)
        CALL f_rms(tmpData,0.,sdev,mask_count)
        rms(i) = sdev
        DEALLOCATE(tmpData)
     ELSE
        rms(i) = 0.
     END IF
  END DO
  DEALLOCATE(flag_mask)
     
  ! last half datapoints
  DO i=nInt-half+1,nInt
     ALLOCATE(flag_mask(nInt-i+half+1))
     flag_mask = flags(i-half:nInt).eq.0
     mask_count = COUNT(flag_mask)
     ! need at least 3 points to compute an rms
     IF(mask_count>half/2) then
        ALLOCATE(tmpData(mask_count))
        tmpData = PACK(datain(i-half:nInt),MASK=flag_mask)
        CALL f_rms(tmpData,0.,sdev,mask_count)
        rms(i) = sdev
        DEALLOCATE(tmpData)
     ELSE
        rms(i) = 0.
     END IF
     DEALLOCATE(flag_mask)
  END DO

  ! now replace zeros with inter- or extrapolated values
  zero_mask = rms.eq.0.
  mask_count = COUNT(zero_mask)
  IF(mask_count>0) THEN
     i = 1
     IF (zero_mask(1)) THEN
        ! zeros at the start: replace with first non-zero value
        first = 2
        DO WHILE (zero_mask(first))
           first = first+1
        END DO
        rms(1:first-1) = rms(first)
        i = first
     END IF
     DO WHILE (i <= nInt)
        prev = i-1
        DO WHILE (.NOT.(zero_mask(i)).AND.(i <= nInt))
           prev = i
           i = i+1
        END DO
        IF (i <= nInt) THEN
           ! rms(i) is 0
           next = i
           DO WHILE (zero_mask(next).AND.(next < nint))
              next = next+1
           END DO
           IF (next.eq.nInt) THEN
              ! zeros at end: replace with last non-zero
              rms(prev+1:nInt) = rms(prev)
           ELSE
              ! some gap in the middle of the time serie
              DO i=prev+1,next-1
                 rms(i) = rms(prev) + (rms(next)-rms(prev))*(i-prev)/(next-prev)
              END DO
           END IF
           i = next
        END IF
        i = i+1
     END DO
  END IF


END SUBROUTINE slidingRms
