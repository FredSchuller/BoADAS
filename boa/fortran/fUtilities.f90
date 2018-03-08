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

subroutine safeExp(x,y,n)
! """
! NAM: safeExp (subroutine)
! DES: Attempt to create a fast exponential for python
! """ 
!INOUT ----------------------------
  integer                      :: n
  real, dimension(n)           :: x,y
  !f2py intent(hide),depend(x) :: n = shape(x,0)
  !f2py intent(in)             :: x
  !f2py intent(out)            :: y
!-----------------------------------
!LOCAL
  integer                      :: i
!-----------------------------------

  do i = 1, n
     y(i) = exp(x(i))
  enddo

end subroutine safeExp

subroutine modelBaseEllipticalGaussian(p,n_p,position,n_position,returned_result)
! """
! NAM: modelBase2dgauss (subroutine)
! DES: Compute a model of a 2D gaussian + gradient 
!      (see corresponding python routine)
! """ 
!INOUT ----------------------------
  integer*8                             :: n_p, n_position
  double precision, dimension(n_p)                :: p
  double precision, dimension(n_position)         :: returned_result
  double precision, dimension(2,n_position)       :: position
  !f2py intent(hide),depend(position) :: n_position = shape(position,2)
  !f2py intent(hide),depend(p)        :: n_p = shape(p,0)
  !f2py intent(in)                    :: p, position
  !f2py intent(out)                   :: returned_result
!-----------------------------------
!LOCAL
  integer*8                           :: i
  double precision                    :: pi
  double precision                    :: xp, yp, U
  double precision                    :: fwhm2sigma, sigma_x, sigma_y
!-----------------------------------
  

  fwhm2sigma = 1./(2*sqrt(2.*log(2.)))
  sigma_x = p(7)*fwhm2sigma
  sigma_y = p(8)*fwhm2sigma

  pi = 3.1415926535897931
  
  do i=1,n_position
     xp = (position(1,i)-p(5))*cos(p(9))-(position(2,i)-p(6))*sin(p(9))
     yp = (position(1,i)-p(5))*sin(p(9))+(position(2,i)-p(6))*cos(p(9))
     U = (xp/sigma_x)**2+(yp/sigma_y)**2
     returned_result(i) = p(1)+p(2)*position(1,i)+p(3)*position(2,i)+p(4)*exp(-U/2)
  enddo
    
end subroutine modelBaseEllipticalGaussian

!=============================================================================
subroutine compress(data, flag, flag_value, dataOut, nData, nOut)
! """
! NAM: compress (subroutine)
! DES: compress array based on mask
! INP: data_input (f): 1-D array with data to be returned where flag = flag_value
!      flag       (i): 1-D array, must be same size as data_input
!      flag_value (i): select data where flag=flag_value
! OUT: tuple (x,nData): 1-D array of same size as data_input where first n elements are
!                   those of data_input where flag = flag_value  
! USE: 
!      Unfortunately it is not possible to pass from fortran to python an output
!      array the size of which is computed in fortran. Therefore we must truncate
!      the returned array at the size of the selected true mask elements.
!
!      Example:
!            input_array = array(range(5),'f')
!            flag_array  = array([0,1,0,1,0],'i')
!            compress_array,nmax = f90.f1.compress(input_array,flag_array,1)
!            compress_array = compress_array[0:nmax]
!
! WRN: On input array bigger than 2100000 elements, this function
!      crash because of the f2py interface
! """ 
!-----------------------------------
!INOUT
  integer                          :: nData, nOut
  real,    dimension(nData)        :: data, dataOut
  integer, dimension(nData)        :: flag
  integer, intent(in)              :: flag_value
  !f2py intent(hide),depend(data)  :: nData=shape(data)
  !f2py intent(in)                 :: data
  !f2py intent(hide,out)           :: nOut, dataOut
  
!-----------------------------------
!LOCAL
  logical, dimension(nData)        :: flag_mask   
  integer                          :: mask_true
!-----------------------------------

  flag_mask  = flag.eq.flag_value     ! make a mask
  mask_true  = COUNT(flag_mask)       ! counts the number of TRUE in mask

  if(mask_true>0) then 
     dataOut = PACK(data,MASK=flag_mask)
     nOut    = mask_true
  else
     dataOut = data
     nOut = 0
  end if

end subroutine compress

!=============================================================================
subroutine icompress(data, flag, flag_value, dataOut, nData, nOut)
! """
! NAM: compress (subroutine)
! DES: compress array based on mask
! INP: data_input (i): 1-D array with data to be returned where flag = flag_value
!      flag       (i): 1-D array, must be same size as data_input
!      flag_value (i): select data where flag=flag_value
! OUT: tuple (x,nData): 1-D array of same size as data_input where first n elements are
!                   those of data_input where flag = flag_value  
! USE: 
!      Unfortunately it is not possible to pass from fortran to python an output
!      array the size of which is computed in fortran. Therefore we must truncate
!      the returned array at the size of the selected true mask elements.
!
!      Example:
!            input_array = array(range(5),'f')
!            flag_array  = array([0,1,0,1,0],'i')
!            compress_array,nmax = f90.f1.compress(input_array,flag_array,1)
!            compress_array = compress_array[0:nmax]
! """ 
!-----------------------------------
!INOUT
  integer                          :: nData, nOut
  integer, dimension(nData)        :: data, dataOut
  integer, dimension(nData)        :: flag
  integer, intent(in)              :: flag_value
  !f2py intent(hide),depend(data)  :: nData=shape(data)
  !f2py intent(in)                 :: data
  !f2py intent(hide,out)           :: nOut, dataOut
  
!-----------------------------------
!LOCAL
  logical, dimension(nData)        :: flag_mask   
  integer                          :: mask_true
!-----------------------------------

  flag_mask  = flag.eq.flag_value     ! make a mask
  mask_true  = COUNT(flag_mask)       ! counts the number of TRUE in mask

  if(mask_true>0) then 
     dataOut = PACK(data,MASK=flag_mask)
     nOut    = mask_true
  else
     dataOut = data
     nOut = 0
  end if

end subroutine icompress

!=============================================================================
subroutine ncompress(data, flag, flag_value, dataOut, nData, nOut)
! """
! NAM: compress (subroutine)
! DES: compress array based on mask
! INP: data_input (f): 1-D array with data to be returned where flag != flag_value
!      flag       (i): 1-D array, must be same size as data_input
!      flag_value (i): select data where flag=flag_value
! OUT: tuple (x,nData): 1-D array of same size as data_input where first n elements are
!                   those of data_input where flag = flag_value  
! USE: 
!      Unfortunately it is not possible to pass from fortran to python an output
!      array the size of which is computed in fortran. Therefore we must truncate
!      the returned array at the size of the selected true mask elements.
!
!      Example:
!            input_array = array(range(5),'f')
!            flag_array  = array([0,1,0,1,0],'i')
!            compress_array,nmax = f90.f1.compress(input_array,flag_array,1)
!            compress_array = compress_array[0:nmax]
! """ 
!-----------------------------------
!INOUT
  integer                          :: nData, nOut
  real,    dimension(nData)        :: data, dataOut
  integer, dimension(nData)        :: flag
  integer, intent(in)              :: flag_value
  !f2py intent(hide),depend(data)  :: nData=shape(data)
  !f2py intent(in)                 :: data
  !f2py intent(hide,out)           :: nOut, dataOut
  
!-----------------------------------
!LOCAL
  logical, dimension(nData)        :: flag_mask   
  integer                          :: mask_true
!-----------------------------------

  flag_mask  = flag.ne.flag_value     ! make a mask
  mask_true  = COUNT(flag_mask)       ! counts the number of TRUE in mask

  if(mask_true>0) then 
     dataOut = PACK(data,MASK=flag_mask)
     nOut    = mask_true
  else
     dataOut = data
     nOut = 0
  end if

end subroutine ncompress




!=============================================================================
subroutine replaceNaN(array,nx,ny)
! """
! DES: replace the NaN in a 2D array by a value below 
!      the minimum value of the array
! """ 
!INOUT ----------------------------
  real, dimension(nx,ny)           :: array
  integer                          :: nx,ny
  !f2py intent(hide),depend(array) :: nx = shape(array,0), ny = shape(array,1)
  !f2py intent(in,out)             :: array
!-----------------------------------
!LOCAL
  logical, dimension(nx,ny)        :: mask
  real                             :: min, max
!-----------------------------------

  ! test for the presence of NaN (IEEE x.eq.y should fail for NaN even if x=y=NaN)
  ! for ifc one need to use the -mp1 option for that to work

  mask = array.eq.array

  min = MINVAL(array, mask=mask)
  max = MAXVAL(array, mask=mask)
  
  ! Usualy 255 colors 
  WHERE(array.ne.array) array = min -(max-min)/255

end subroutine replaceNaN

!=============================================================================
subroutine maskNaN(array,mask,nx)
! """
! DES: compute a mask where input 1D array is NaN
! """ 
!INOUT ----------------------------
  real, dimension(nx)              :: array
  integer*4, dimension(nx)         :: mask
  integer                          :: nx
  !f2py intent(hide),depend(array) :: nx = shape(array,0)
  !f2py intent(in)                 :: array
  !f2py intent(out)                :: mask
!-----------------------------------
  mask = array.ne.array
end subroutine maskNaN

!=============================================================================
subroutine MreplaceNaN(array,nx,ny,nz)
! """
! DES: replace the NaN in a 3D array by a value below 
!      the minimum value of the array
! """ 
!INOUT ----------------------------
  real, dimension(nx,ny,nz)        :: array
  integer                          :: nx,ny,nz
  !f2py intent(hide),depend(array) :: nx = shape(array,0), ny = shape(array,1), nz = shape(array,2)
  !f2py intent(in,out)             :: array
!-----------------------------------
!LOCAL
  logical, dimension(nx,ny,nz)     :: mask
  real                             :: min, max
!-----------------------------------

  ! test for the presence of NaN (IEEE x.eq.y should fail for NaN even if x=y=NaN)
  ! for ifc one need to use the -mp1 option for that to work

  mask = array.eq.array

  min = MINVAL(array, mask=mask)
  max = MAXVAL(array, mask=mask)
  
  ! Usualy 255 colors 
  WHERE(array.ne.array) array = min -(max-min)/255

end subroutine MreplaceNaN


subroutine matrixmultiply(m1,m2,prod,n1,nx,n2)
! """
! DES: matrix multiplication 
!     
! """ 
!INOUT ----------------------------
  integer                   :: n1, nx, n2
  double precision, dimension(n1,nx)    :: m1
  double precision, dimension(nx,n2)    :: m2
  double precision, dimension(n1,n2)    :: prod
  !f2py intent(hide),depend(m1)    :: n1=shape(m1,0),nx=shape(m1,1)
  !f2py intent(hide),depend(m2)    :: n2=shape(m2,1)
  !f2py intent(in)                 :: m1,m2
  !f2py intent(out)                :: prod
!-----------------------------------
!LOCAL
 
!-----------------------------------

  prod=matmul(m1,m2)

end subroutine matrixmultiply
