SUBROUTINE wavetransform1d(np,data,ns,wtrans)
! """
! NAM: waveTransform1D (subroutine)
! DES: Computes the wavelet transform of an image using the 'a trou' algorithm and a B3-spline wavelet
! INP: data   : input rank-1 array('f') to be decomposed
!      ns     : number of wavelet scales
! OUT: wtrans : rank-2 array('f') with shape (len(data),ns). 
! """ 
!-----------------------------------
!INOUT
  integer*8                                    :: np
  real*4, dimension(np)                        :: data
  integer*4                                    :: ns
  real*4, dimension(np,ns)                     :: wtrans
  integer*4,parameter                          :: filtsize=5
  real*4, dimension(filtsize)                  :: filter = (/ 1./16. , 1./4.,  3./8.,  1./4., 1./16. /)
  !f2py integer*8, intent(hide),depend(data)   :: np=len(data)
  !f2py real*4, intent(in)                     :: data
  !f2py integer*4, intent(in)                  :: ns
  !f2py real*4, intent(out)                    :: wtrans
!-----------------------------------
!LOCAL
  integer*8, dimension(filtsize)               :: x_filt
  integer*4                                    :: s
  integer*8                                    :: i, j, x
  real*4, dimension(np)                        :: smoothed
  real*4, dimension(np)                        :: c
!-----------------------------------
 
  ! Initialize pixel filter
  DO i=1,filtsize
    x_filt(i) = i - FLOOR( filtsize/2. ) -1
  END DO
  
  s = FLOOR( LOG(real(np)) / LOG(2.) )
  IF (ns > s) THEN
   	WRITE(*,*) "Warning: Maximum number of scales for input is ",s
   	WRITE(*,*) "Warning: Returning empty array ..."	
   	RETURN
  END IF
  
  DO i = 1,np
  	c(i) = data(i)
  END DO
  DO s = 0,(ns-2)
	DO i = 1,np
		smoothed(i) = 0.
		DO j = 1,filtsize
		   x = x_filt(j) + i
		   IF (x < 1) THEN
		     x = 2 - x
		   END IF
		   IF (x > np) THEN
		     x = 2*np -x
		   END IF
		  smoothed(i) =  smoothed(i) + filter(j)*c(x)
		END DO
	        wtrans(i,s+1) = c(i) - smoothed(i)
	END DO
	DO i = 1,np
	  	c(i) = smoothed(i)
	ENd DO
	DO i = 1,filtsize
		x_filt(i)=x_filt(i)*2
	END DO
  END DO

  wtrans(:,ns) = smoothed

END SUBROUTINE wavetransform1d

SUBROUTINE wavetransform2d(nx,ny,data,ns,wtrans)
! """
! NAM: waveTransform2D (subroutine)
! DES: Computes the wavelet transform of an image using the 'a trou' algorithm and a B3-spline wavelet
! INP: data   : input rank-2array('f') to be decomposed
!      ns     : number of wavelet scales
! OUT: wtrans : rank-3 array('f') with shape (len(data),ns). 
! """ 
!-----------------------------------
!INOUT
  integer*8                                    :: nx
  integer*8                                    :: ny
  real*4, dimension(nx,ny)                     :: data
  integer*4                                    :: ns
  real*4, dimension(nx,ny,ns)                  :: wtrans
  integer*4,parameter                          :: filtsize=5
  real*4, dimension(filtsize)                  :: filter = (/ 1./16. , 1./4.,  3./8.,  1./4., 1./16. /)
  real*4, dimension(filtsize,filtsize)         :: filter2d
  !f2py integer*8, intent(hide),depend(data)   :: nx=shape(data,0),ny=shape(data,1)
  !f2py real*4, intent(in)                     :: data
  !f2py integer*4, intent(in)                  :: ns
  !f2py real*4, intent(out)                    :: wtrans
!-----------------------------------
!LOCAL
  integer*8, dimension(filtsize)      	       :: pix_filt
  integer*4                                    :: s
  integer*8                                    :: i, j, k, l, y, x
  real*4, dimension(nx,ny)                     :: smoothed
  real*4, dimension(nx,ny)                     :: c
!-----------------------------------
 
  ! Initialize pixel filters
  DO i=1,filtsize
  	pix_filt(i) = i - FLOOR( filtsize/2. ) -1
  	DO j=1,filtsize
    		filter2d(i,j) = filter(i)*filter(j)
  	END DO
  END DO
  
  s = FLOOR( LOG(real(min(nx,ny))) / LOG(2.) )
  IF (ns > s) THEN
   	WRITE(*,*) "Warning: Maximum number of scales for input is ",s
   	WRITE(*,*) "Warning: Returning empty array ..."	
   	RETURN
  END IF
  
  DO i = 1,nx
  	DO j = 1,ny
  		c(i,j) = data(i,j)
  	END DO
  END DO

  DO s = 0,(ns-2)

	DO i = 1,nx
		DO j = 1,ny
			smoothed(i,j) = 0.
			DO k = 1,filtsize
				DO l = 1,filtsize
		       			x = pix_filt(k) + i
		       			y = pix_filt(l) + j
		       			IF (x < 1) THEN
		          			x = 2 - x
		      		 	END IF
		       			IF (x > nx) THEN
		          			x = 2*nx -x
		      			END IF
		       			IF (y < 1) THEN
		          			y = 2 - y
		      		 	END IF
		       			IF (y > ny) THEN
		          			y = 2*ny -y
		      			END IF
		       			smoothed(i,j) =  smoothed(i,j) + filter2d(k,l)*c(x,y)
		 		END DO
		 	END DO
	          	wtrans(i,j,s+1) = c(i,j) - smoothed(i,j)
	     	END DO
	END DO

	DO i = 1,nx
		DO j = 1,ny
	  		c(i,j) = smoothed(i,j)
		END DO
	END DO

	DO i = 1,filtsize
		pix_filt(i) = pix_filt(i)*2
	END DO

  END DO

  wtrans(:,:,ns) = smoothed

END SUBROUTINE wavetransform2d
