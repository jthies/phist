/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
!> \file gemm_kernels_sC.f90
!! Fast parallel BLAS-gemm like subroutines for different blocksizes for mvecT_times_mvec in mvec_module
!! \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
!!

#include "phist_config.h"

subroutine dgemm_sC_1_1(nrows,v,w,M)
  implicit none
  integer, intent(in) :: nrows
  real(kind=8), intent(in) :: v(1,nrows)
  real(kind=8), intent(in) :: w(1,nrows)
  real(kind=8), intent(out) :: M(1,1)
  integer :: i, j
!dir$ assume_aligned v:64, w:64, M:64

  M = 0.
!$omp parallel do reduction(+:M) schedule(static)
  do i = 1, nrows, 1
    do j = 1, 1, 1
      M(:,j) = M(:,j) + v(:,i)*w(j,i)
    end do
  end do

end subroutine dgemm_sC_1_1

subroutine dgemm_sC_2_2(nrows,v,w,M)
  implicit none
  integer, intent(in) :: nrows
  real(kind=8), intent(in) :: v(2,nrows)
  real(kind=8), intent(in) :: w(2,nrows)
  real(kind=8), intent(out) :: M(2,2)
  integer :: i, j
!dir$ assume_aligned v:64, w:64, M:64

  M = 0.
!$omp parallel do reduction(+:M) schedule(static)
  do i = 1, nrows, 1
    do j = 1, 2, 1
      M(:,j) = M(:,j) + v(:,i)*w(j,i)
    end do
  end do

end subroutine dgemm_sC_2_2

subroutine dgemm_sC_4_4(nrows,v,w,M)
  implicit none
  integer, intent(in) :: nrows
  real(kind=8), intent(in) :: v(4,nrows)
  real(kind=8), intent(in) :: w(4,nrows)
  real(kind=8), intent(out) :: M(4,4)
  integer :: i, j
!dir$ assume_aligned v:64, w:64, M:64

  M = 0.
!$omp parallel do reduction(+:M) schedule(static)
  do i = 1, nrows, 1
    do j = 1, 4, 1
      M(:,j) = M(:,j) + v(:,i)*w(j,i)
    end do
  end do

end subroutine dgemm_sC_4_4

subroutine dgemm_sC_8_8(nrows,v,w,M)
  implicit none
  integer, intent(in) :: nrows
  real(kind=8), intent(in) :: v(8,nrows)
  real(kind=8), intent(in) :: w(8,nrows)
  real(kind=8), intent(out) :: M(8,8)
  integer :: i, j
!dir$ assume_aligned v:64, w:64, M:64

  M = 0.
!$omp parallel do reduction(+:M) schedule(static)
  do i = 1, nrows, 1
    do j = 1, 8, 1
      M(:,j) = M(:,j) + v(:,i)*w(j,i)
    end do
  end do

end subroutine dgemm_sC_8_8


subroutine dgemm_sC_1(nrows,nvecw,v,w,M)
  implicit none
  integer, intent(in) :: nrows, nvecw
  real(kind=8), intent(in) :: v(1,nrows)
  real(kind=8), intent(in) :: w(nvecw,nrows)
  real(kind=8), intent(out) :: M(1,nvecw)
  integer :: i, j
!dir$ assume_aligned v:64, w:64, M:64

  M = 0.
!$omp parallel do reduction(+:M) schedule(static)
  do i = 1, nrows, 1
    do j = 1, nvecw, 1
      M(:,j) = M(:,j) + v(:,i)*w(j,i)
    end do
  end do

end subroutine dgemm_sC_1

subroutine dgemm_sC_2(nrows,nvecw,v,w,M)
  implicit none
  integer, intent(in) :: nrows, nvecw
  real(kind=8), intent(in) :: v(2,nrows)
  real(kind=8), intent(in) :: w(nvecw,nrows)
  real(kind=8), intent(out) :: M(2,nvecw)
  integer :: i, j
!dir$ assume_aligned v:64, w:64, M:64

  M = 0.
!$omp parallel do reduction(+:M) schedule(static)
  do i = 1, nrows, 1
    do j = 1, nvecw, 1
      M(:,j) = M(:,j) + v(:,i)*w(j,i)
    end do
  end do

end subroutine dgemm_sC_2

subroutine dgemm_sC_4(nrows,nvecw,v,w,M)
  implicit none
  integer, intent(in) :: nrows, nvecw
  real(kind=8), intent(in) :: v(4,nrows)
  real(kind=8), intent(in) :: w(nvecw,nrows)
  real(kind=8), intent(out) :: M(4,nvecw)
  integer :: i, j
!dir$ assume_aligned v:64, w:64, M:64

  M = 0.
!$omp parallel do reduction(+:M) schedule(static)
  do i = 1, nrows, 1
    do j = 1, nvecw, 1
      M(:,j) = M(:,j) + v(:,i)*w(j,i)
    end do
  end do

end subroutine dgemm_sC_4

subroutine dgemm_sC_8(nrows,nvecw,v,w,M)
  implicit none
  integer, intent(in) :: nrows, nvecw
  real(kind=8), intent(in) :: v(8,nrows)
  real(kind=8), intent(in) :: w(nvecw,nrows)
  real(kind=8), intent(out) :: M(8,nvecw)
  integer :: i, j
!dir$ assume_aligned v:64, w:64, M:64

  M = 0.
!$omp parallel do reduction(+:M) schedule(static)
  do i = 1, nrows, 1
    do j = 1, nvecw, 1
      M(:,j) = M(:,j) + v(:,i)*w(j,i)
    end do
  end do

end subroutine dgemm_sC_8


subroutine dgemm_sC_strided_1(nrows,nvecw,v,w,ldw,M)
  implicit none
  integer, intent(in) :: nrows, nvecw, ldw
  real(kind=8), intent(in) :: v(1,nrows)
  real(kind=8), intent(in) :: w(ldw,nrows)
  real(kind=8), intent(out) :: M(1,nvecw)
  integer :: i, j
!dir$ assume_aligned v:64, w:8, M:64

  M = 0.
!$omp parallel do reduction(+:M) schedule(static)
  do i = 1, nrows, 1
    do j = 1, nvecw, 1
      M(:,j) = M(:,j) + v(:,i)*w(j,i)
    end do
  end do

end subroutine dgemm_sC_strided_1

subroutine dgemm_sC_strided_2(nrows,nvecw,v,w,ldw,M)
  implicit none
  integer, intent(in) :: nrows, nvecw, ldw
  real(kind=8), intent(in) :: v(2,nrows)
  real(kind=8), intent(in) :: w(ldw,nrows)
  real(kind=8), intent(out) :: M(2,nvecw)
  integer :: i, j
!dir$ assume_aligned v:64, w:8, M:64

  M = 0.
!$omp parallel do reduction(+:M) schedule(static)
  do i = 1, nrows, 1
    do j = 1, nvecw, 1
      M(:,j) = M(:,j) + v(:,i)*w(j,i)
    end do
  end do

end subroutine dgemm_sC_strided_2

subroutine dgemm_sC_strided_4(nrows,nvecw,v,w,ldw,M)
  implicit none
  integer, intent(in) :: nrows, nvecw, ldw
  real(kind=8), intent(in) :: v(4,nrows)
  real(kind=8), intent(in) :: w(ldw,nrows)
  real(kind=8), intent(out) :: M(4,nvecw)
  integer :: i, j
!dir$ assume_aligned v:64, w:8, M:64

  M = 0.
!$omp parallel do reduction(+:M) schedule(static)
  do i = 1, nrows, 1
    do j = 1, nvecw, 1
      M(:,j) = M(:,j) + v(:,i)*w(j,i)
    end do
  end do

end subroutine dgemm_sC_strided_4

subroutine dgemm_sC_strided_8(nrows,nvecw,v,w,ldw,M)
  implicit none
  integer, intent(in) :: nrows, nvecw, ldw
  real(kind=8), intent(in) :: v(8,nrows)
  real(kind=8), intent(in) :: w(ldw,nrows)
  real(kind=8), intent(out) :: M(8,nvecw)
  integer :: i, j
!dir$ assume_aligned v:64, w:8, M:64

  M = 0.
!$omp parallel do reduction(+:M) schedule(static)
  do i = 1, nrows, 1
    do j = 1, nvecw, 1
      M(:,j) = M(:,j) + v(:,i)*w(j,i)
    end do
  end do

end subroutine dgemm_sC_strided_8


subroutine dgemm_sC_generic(nrows,nvecv,nvecw,v,ldv,w,ldw,M)
  implicit none
  integer, intent(in) :: nrows, nvecv, nvecw, ldv, ldw
  real(kind=8), intent(in) :: v(ldv,nrows)
  real(kind=8), intent(in) :: w(ldw,nrows)
  real(kind=8), intent(out) :: M(nvecv,nvecw)
  integer :: i, j
!dir$ assume_aligned v:8, w:8, M:64

  M = 0.
!$omp parallel do reduction(+:M) schedule(static)
  do i = 1, nrows, 1
    do j = 1, nvecw, 1
      M(:,j) = M(:,j) + v(1:nvecv,i)*w(j,i)
    end do
  end do

end subroutine dgemm_sC_generic


