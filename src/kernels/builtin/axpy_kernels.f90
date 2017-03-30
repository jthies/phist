/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
!> \file axpy_kernels.f90
!! Fast parallel BLAS-axpy style subroutines for different blocksizes for mvec_module,
!! cases which may emplyoy non-temporal stores (NT) are delegated to
!! axpy_kernels_nt.c
!! \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
!!

subroutine dset_1(nrows, y, val)
  implicit none
  integer, intent(in) :: nrows
  real(kind=8), intent(out) :: y(nrows)
  real(kind=8), intent(in) :: val
  integer :: i
!dir$ assume_aligned y:64

!$omp parallel do schedule(static)
  do i = 1, nrows, 1
    y(i) = val
  end do
end subroutine dset_1

subroutine dset_general(nvec, nrows, y, ldy, val)
  implicit none
  integer, intent(in) :: nvec, nrows, ldy
  real(kind=8), intent(out) :: y(ldy,*)
  real(kind=8), intent(in) :: val
  integer :: i
!dir$ assume_aligned y:8

!$omp parallel do schedule(static)
  do i = 1, nrows, 1
    y(1:nvec,i) = val
  end do
end subroutine dset_general

subroutine dcopy_1(nrows, x, y)
  implicit none
  integer, intent(in) :: nrows
  real(kind=8), intent(in) :: x(nrows)
  real(kind=8), intent(out) :: y(nrows)
  integer :: i
!dir$ assume_aligned x:64, y:64

!$omp parallel do schedule(static)
  do i = 1, nrows, 1
    y(i) = x(i)
  end do
end subroutine dcopy_1

subroutine dcopy_general(nvec, nrows, x, ldx, y, ldy)
  implicit none
  integer, intent(in) :: nvec, nrows, ldx, ldy
  real(kind=8), intent(in) :: x(ldx,*)
  real(kind=8), intent(out) :: y(ldy,*)
  integer :: i
!dir$ assume_aligned x:8, y:8

!$omp parallel do schedule(static)
  do i = 1, nrows, 1
    y(1:nvec,i) = x(1:nvec,i)
  end do
end subroutine dcopy_general


subroutine dscal_1(nrows, alpha, x)
  implicit none
  integer, intent(in) :: nrows
  real(kind=8), intent(in) :: alpha(1)
  real(kind=8), intent(inout) :: x(1,nrows)
  integer :: i
!dir$ assume_aligned x:64

!$omp parallel do schedule(static)
  do i = 1, nrows, 1
    x(:,i) = alpha(:)*x(:,i)
  end do
end subroutine dscal_1

subroutine dscal_2(nrows, alpha, x)
  implicit none
  integer, intent(in) :: nrows
  real(kind=8), intent(in) :: alpha(2)
  real(kind=8), intent(inout) :: x(2,nrows)
  integer :: i
!dir$ assume_aligned x:64

!$omp parallel do schedule(static)
  do i = 1, nrows, 1
    x(:,i) = alpha(:)*x(:,i)
  end do
end subroutine dscal_2

subroutine dscal_4(nrows, alpha, x)
  implicit none
  integer, intent(in) :: nrows
  real(kind=8), intent(in) :: alpha(4)
  real(kind=8), intent(inout) :: x(4,nrows)
  integer :: i
!dir$ assume_aligned x:64

!$omp parallel do schedule(static)
  do i = 1, nrows, 1
    x(:,i) = alpha(:)*x(:,i)
  end do
end subroutine dscal_4

subroutine dscal_8(nrows, alpha, x)
  implicit none
  integer, intent(in) :: nrows
  real(kind=8), intent(in) :: alpha(8)
  real(kind=8), intent(inout) :: x(8,nrows)
  integer :: i
!dir$ assume_aligned x:64

!$omp parallel do schedule(static)
  do i = 1, nrows, 1
    x(:,i) = alpha(:)*x(:,i)
  end do
end subroutine dscal_8


subroutine dscal_strided_1(nrows, alpha, x, ldx)
  implicit none
  integer, intent(in) :: nrows, ldx
  real(kind=8), intent(in) :: alpha(1)
  real(kind=8), intent(inout) :: x(ldx,*)
  integer :: i
!dir$ assume_aligned x:8

!$omp parallel do schedule(static)
  do i = 1, nrows, 1
    x(1,i) = alpha(1)*x(1,i)
  end do
end subroutine dscal_strided_1

subroutine dscal_strided_2(nrows, alpha, x, ldx)
  implicit none
  integer, intent(in) :: nrows, ldx
  real(kind=8), intent(in) :: alpha(2)
  real(kind=8), intent(inout) :: x(ldx,*)
  integer :: i
!dir$ assume_aligned x:8

!$omp parallel do schedule(static)
  do i = 1, nrows, 1
    x(1:2,i) = alpha(:)*x(1:2,i)
  end do
end subroutine dscal_strided_2

subroutine dscal_strided_4(nrows, alpha, x, ldx)
  implicit none
  integer, intent(in) :: nrows, ldx
  real(kind=8), intent(in) :: alpha(4)
  real(kind=8), intent(inout) :: x(ldx,*)
  integer :: i
!dir$ assume_aligned x:8

!$omp parallel do schedule(static)
  do i = 1, nrows, 1
    x(1:4,i) = alpha(:)*x(1:4,i)
  end do
end subroutine dscal_strided_4

subroutine dscal_strided_8(nrows, alpha, x, ldx)
  implicit none
  integer, intent(in) :: nrows, ldx
  real(kind=8), intent(in) :: alpha(8)
  real(kind=8), intent(inout) :: x(ldx,*)
  integer :: i
!dir$ assume_aligned x:8

!$omp parallel do schedule(static)
  do i = 1, nrows, 1
    x(1:8,i) = alpha(:)*x(1:8,i)
  end do
end subroutine dscal_strided_8

subroutine dscal_general(nrows, nvec, alpha, x, ldx)
  implicit none
  integer, intent(in) :: nrows, nvec, ldx
  real(kind=8), intent(in) :: alpha(nvec)
  real(kind=8), intent(inout) :: x(ldx,*)
  integer :: i
!dir$ assume_aligned x:8

!$omp parallel do schedule(static)
  do i = 1, nrows, 1
    x(1:nvec,i) = alpha(:)*x(1:nvec,i)
  end do
end subroutine dscal_general


subroutine daxpby_1(nrows, alpha, x, beta, y)
  implicit none
  integer, intent(in) :: nrows
  real(kind=8), intent(in) :: alpha
  real(kind=8), intent(in) :: beta
  real(kind=8), intent(in) :: x(nrows)
  real(kind=8), intent(inout) :: y(nrows)
  integer :: i
!dir$ assume_aligned x:64, y:64

  if( beta .eq. 0 ) then
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      y(i) = alpha*x(i)
    end do
  else
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      y(i) = alpha*x(i) + beta*y(i)
    end do
  end if
end subroutine daxpby_1


subroutine daxpy_NT_2(nrows, alpha, x, y)
  implicit none
  integer, intent(in) :: nrows
  real(kind=8), intent(in) :: alpha(2)
  real(kind=8), intent(in) :: x(2,nrows)
  real(kind=8), intent(inout) :: y(2,nrows)

  interface
    subroutine daxpy_NT_2_c(nrows,alpha,x,y) bind(C)
      use, intrinsic :: iso_c_binding
      integer(C_INT), value :: nrows
      real(C_DOUBLE), intent(in) :: alpha(*)
      real(C_DOUBLE), intent(in) :: x(*)
      real(C_DOUBLE), intent(out) :: y(*)
    end subroutine daxpy_NT_2_c
  end interface


  call daxpy_NT_2_c(nrows,alpha,x,y)

end subroutine daxpy_NT_2

subroutine daxpy_NT_4(nrows, alpha, x, y)
  implicit none
  integer, intent(in) :: nrows
  real(kind=8), intent(in) :: alpha(4)
  real(kind=8), intent(in) :: x(4,nrows)
  real(kind=8), intent(inout) :: y(4,nrows)

  interface
    subroutine daxpy_NT_4_c(nrows,alpha,x,y) bind(C)
      use, intrinsic :: iso_c_binding
      integer(C_INT), value :: nrows
      real(C_DOUBLE), intent(in) :: alpha(*)
      real(C_DOUBLE), intent(in) :: x(*)
      real(C_DOUBLE), intent(out) :: y(*)
    end subroutine daxpy_NT_4_c
  end interface

  call daxpy_NT_4_c(nrows,alpha,x,y)

end subroutine daxpy_NT_4

subroutine daxpy_NT_8(nrows, alpha, x, y)
  implicit none
  integer, intent(in) :: nrows
  real(kind=8), intent(in) :: alpha(8)
  real(kind=8), intent(in) :: x(8,nrows)
  real(kind=8), intent(inout) :: y(8,nrows)

  interface
    subroutine daxpy_NT_8_c(nrows,alpha,x,y) bind(C)
      use, intrinsic :: iso_c_binding
      integer(C_INT), value :: nrows
      real(C_DOUBLE), intent(in) :: alpha(*)
      real(C_DOUBLE), intent(in) :: x(*)
      real(C_DOUBLE), intent(out) :: y(*)
    end subroutine daxpy_NT_8_c
  end interface


  call daxpy_NT_8_c(nrows,alpha,x,y)

end subroutine daxpy_NT_8


subroutine daxpy_NT_strided_2(nrows, alpha, x, ldx, y, ldy)
  implicit none
  integer, intent(in) :: nrows, ldx, ldy
  real(kind=8), intent(in) :: alpha(2)
  real(kind=8), intent(in) :: x(2,nrows)
  real(kind=8), intent(inout) :: y(2,nrows)

  interface
    subroutine daxpy_NT_strided_2_c(nrows,alpha,x,ldx,y,ldy) bind(C)
      use, intrinsic :: iso_c_binding
      integer(C_INT), value :: nrows, ldx, ldy
      real(C_DOUBLE), intent(in) :: alpha(*)
      real(C_DOUBLE), intent(in) :: x(*)
      real(C_DOUBLE), intent(out) :: y(*)
    end subroutine daxpy_NT_strided_2_c
  end interface


  call daxpy_NT_strided_2_c(nrows,alpha,x,ldx,y,ldy)

end subroutine daxpy_NT_strided_2

subroutine daxpy_NT_strided_4(nrows, alpha, x, ldx, y, ldy)
  implicit none
  integer, intent(in) :: nrows, ldx, ldy
  real(kind=8), intent(in) :: alpha(4)
  real(kind=8), intent(in) :: x(4,nrows)
  real(kind=8), intent(inout) :: y(4,nrows)

  interface
    subroutine daxpy_NT_strided_4_c(nrows,alpha,x,ldx,y,ldy) bind(C)
      use, intrinsic :: iso_c_binding
      integer(C_INT), value :: nrows, ldx, ldy
      real(C_DOUBLE), intent(in) :: alpha(*)
      real(C_DOUBLE), intent(in) :: x(*)
      real(C_DOUBLE), intent(out) :: y(*)
    end subroutine daxpy_NT_strided_4_c
  end interface

  call daxpy_NT_strided_4_c(nrows,alpha,x,ldx,y,ldy)

end subroutine daxpy_NT_strided_4

subroutine daxpy_NT_strided_8(nrows, alpha, x, ldx, y, ldy)
  implicit none
  integer, intent(in) :: nrows, ldx, ldy
  real(kind=8), intent(in) :: alpha(8)
  real(kind=8), intent(in) :: x(8,nrows)
  real(kind=8), intent(inout) :: y(8,nrows)

  interface
    subroutine daxpy_NT_strided_8_c(nrows,alpha,x,ldx,y,ldy) bind(C)
      use, intrinsic :: iso_c_binding
      integer(C_INT), value :: nrows, ldx, ldy
      real(C_DOUBLE), intent(in) :: alpha(*)
      real(C_DOUBLE), intent(in) :: x(*)
      real(C_DOUBLE), intent(out) :: y(*)
    end subroutine daxpy_NT_strided_8_c
  end interface


  call daxpy_NT_strided_8_c(nrows,alpha,x,ldx,y,ldy)

end subroutine daxpy_NT_strided_8


subroutine daxpby_2(nrows, alpha, x, beta, y)
  implicit none
  integer, intent(in) :: nrows
  real(kind=8), intent(in) :: alpha(2)
  real(kind=8), intent(in) :: beta
  real(kind=8), intent(in) :: x(2,nrows)
  real(kind=8), intent(inout) :: y(2,nrows)
  integer :: i
!dir$ assume_aligned x:64, y:64

  if( beta .eq. 0 ) then
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      y(:,i) = alpha(:)*x(:,i)
    end do
  else
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      y(:,i) = alpha(:)*x(:,i) + beta*y(:,i)
    end do
  end if
end subroutine daxpby_2

subroutine daxpby_4(nrows, alpha, x, beta, y)
  implicit none
  integer, intent(in) :: nrows
  real(kind=8), intent(in) :: alpha(4)
  real(kind=8), intent(in) :: beta
  real(kind=8), intent(in) :: x(4,nrows)
  real(kind=8), intent(inout) :: y(4,nrows)
  integer :: i
!dir$ assume_aligned x:64, y:64

  if( beta .eq. 0 ) then
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      y(:,i) = alpha(:)*x(:,i)
    end do
  else
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      y(:,i) = alpha(:)*x(:,i) + beta*y(:,i)
    end do
  end if
end subroutine daxpby_4

subroutine daxpby_8(nrows, alpha, x, beta, y)
  implicit none
  integer, intent(in) :: nrows
  real(kind=8), intent(in) :: alpha(8)
  real(kind=8), intent(in) :: beta
  real(kind=8), intent(in) :: x(8,nrows)
  real(kind=8), intent(inout) :: y(8,nrows)
  integer :: i
!dir$ assume_aligned x:64, y:64

  if( beta .eq. 0 ) then
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      y(:,i) = alpha(:)*x(:,i) + beta*y(:,i)
    end do
  else
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      y(:,i) = alpha(:)*x(:,i) + beta*y(:,i)
    end do
  end if
end subroutine daxpby_8


subroutine daxpby_strided_1(nrows, alpha, x, ldx, beta, y, ldy)
  implicit none
  integer, intent(in) :: nrows, ldx, ldy
  real(kind=8), intent(in) :: alpha
  real(kind=8), intent(in) :: beta
  real(kind=8), intent(in) :: x(ldx,*)
  real(kind=8), intent(inout) :: y(ldy,*)
  integer :: i
!dir$ assume_aligned x:8, y:8

  if( beta .eq. 0 ) then
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      y(1,i) = alpha*x(1,i)
    end do
  else
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      y(1,i) = alpha*x(1,i) + beta*y(1,i)
    end do
  end if
end subroutine daxpby_strided_1

subroutine daxpby_strided_2(nrows, alpha, x, ldx, beta, y, ldy)
  implicit none
  integer, intent(in) :: nrows, ldx, ldy
  real(kind=8), intent(in) :: alpha(2)
  real(kind=8), intent(in) :: beta
  real(kind=8), intent(in) :: x(ldx,*)
  real(kind=8), intent(inout) :: y(ldy,*)
  integer :: i
!dir$ assume_aligned x:8, y:8

  if( beta .eq. 0 ) then
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      y(1:2,i) = alpha*x(1:2,i)
    end do
  else
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      y(1:2,i) = alpha*x(1:2,i) + beta*y(1:2,i)
    end do
  end if
end subroutine daxpby_strided_2

subroutine daxpby_strided_4(nrows, alpha, x, ldx, beta, y, ldy)
  implicit none
  integer, intent(in) :: nrows, ldx, ldy
  real(kind=8), intent(in) :: alpha(4)
  real(kind=8), intent(in) :: beta
  real(kind=8), intent(in) :: x(ldx,*)
  real(kind=8), intent(inout) :: y(ldy,*)
  integer :: i
!dir$ assume_aligned x:8, y:8

  if( beta .eq. 0 ) then
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      y(1:4,i) = alpha*x(1:4,i)
    end do
  else
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      y(1:4,i) = alpha*x(1:4,i) + beta*y(1:4,i)
    end do
  end if
end subroutine daxpby_strided_4

subroutine daxpby_strided_8(nrows, alpha, x, ldx, beta, y, ldy)
  implicit none
  integer, intent(in) :: nrows, ldx, ldy
  real(kind=8), intent(in) :: alpha(8)
  real(kind=8), intent(in) :: beta
  real(kind=8), intent(in) :: x(ldx,*)
  real(kind=8), intent(inout) :: y(ldy,*)
  integer :: i
!dir$ assume_aligned x:8, y:8

  if( beta .eq. 0 ) then
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      y(1:8,i) = alpha*x(1:8,i)
    end do
  else
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      y(1:8,i) = alpha*x(1:8,i) + beta*y(1:8,i)
    end do
  end if
end subroutine daxpby_strided_8

subroutine daxpby_generic(nrows, nvec, alpha, x, ldx, beta, y, ldy)
  implicit none
  integer, intent(in) :: nrows, nvec, ldx, ldy
  real(kind=8), intent(in) :: alpha(nvec)
  real(kind=8), intent(in) :: beta
  real(kind=8), intent(in) :: x(ldx,*)
  real(kind=8), intent(inout) :: y(ldy,*)
  integer :: i
!dir$ assume_aligned x:8, y:8

  if( beta .eq. 0 ) then
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      y(1:nvec,i) = alpha*x(1:nvec,i)
    end do
  else
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      y(1:nvec,i) = alpha*x(1:nvec,i) + beta*y(1:nvec,i)
    end do
  end if
end subroutine daxpby_generic

