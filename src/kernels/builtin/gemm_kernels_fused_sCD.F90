/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/

!> \file gemm_fused_kernels_sCD.f90
!! Fast parallel fused BLAS-gemm subroutines for different blocksizes for mvecT_times_mvec_times_sdMat_inplace in mvec_module
!! \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
!!

#include "phist_config_fortran.h"

subroutine dgemm_fused_sCD_1_self(nrows,w,N,M)
  implicit none
  integer, intent(in) :: nrows
  real(kind=8), intent(inout) :: w(1,nrows)
  real(kind=8), intent(in) :: N(1,1)
  real(kind=8), intent(out) :: M(1,1)
  real(kind=8) :: work(1)
  integer :: i, j
!dir$ assume_aligned w:64, M:64, N:64

  M = 0.
!$omp parallel do reduction(+:M) private(work) schedule(static)
  do i = 1, nrows, 1
    work = w(:,i)
    w(1,i) = sum(N(1,1)*work)
    M = M + w(1,i)*w(1,i)
  end do

end subroutine dgemm_fused_sCD_1_self

subroutine dgemm_fused_sCD_2_self(nrows,w,N,M)
  implicit none
  integer, intent(in) :: nrows
  real(kind=8), intent(inout) :: w(2,nrows)
  real(kind=8), intent(in) :: N(2,2)
  real(kind=8), intent(out) :: M(2,2)
  real(kind=8) :: work(2)
  integer :: i, j
!dir$ assume_aligned w:64, M:64, N:64

  M = 0.
!$omp parallel do reduction(+:M) private(work) schedule(static)
  do i = 1, nrows, 1
    work = w(:,i)
    do j = 1, 2, 1
      w(j,i) = sum(N(:,j)*work(:))
      M(:j,j) = M(:j,j) + w(:j,i)*w(j,i)
    end do
  end do

  M(2,1) = M(1,2)

end subroutine dgemm_fused_sCD_2_self

subroutine dgemm_fused_sCD_4_self(nrows,w,N,M)
  implicit none
  integer, intent(in) :: nrows
  real(kind=8), intent(inout) :: w(4,nrows)
  real(kind=8), intent(in) :: N(4,4)
  real(kind=8), intent(out) :: M(4,4)
  real(kind=8) :: work(4)
  integer :: i, j
!dir$ assume_aligned w:64, M:64, N:64

  M = 0.
!$omp parallel do reduction(+:M) private(work) schedule(static)
  do i = 1, nrows, 1
    work = w(:,i)
    do j = 1, 4, 1
      w(j,i) = sum(N(:,j)*work(:))
      M(:j,j) = M(:j,j) + w(:j,i)*w(j,i)
    end do
  end do

  do i = 1, 4, 1
    do j = i+1, 4, 1
      M(j,i) = M(i,j)
    end do
  end do

end subroutine dgemm_fused_sCD_4_self

subroutine dgemm_fused_sCD_8_self(nrows,w,N,M)
  implicit none
  integer, intent(in) :: nrows
  real(kind=8), intent(inout) :: w(8,nrows)
  real(kind=8), intent(in) :: N(8,8)
  real(kind=8), intent(out) :: M(8,8)
  real(kind=8) :: work(8)
  integer :: i, j
!dir$ assume_aligned w:64, M:64, N:64

  M = 0.
!$omp parallel do reduction(+:M) private(work) schedule(static)
  do i = 1, nrows, 1
    work = w(:,i)
    do j = 1, 8, 1
      w(j,i) = sum(N(:,j)*work(:))
      M(:,j) = M(:,j) + w(:,i)*w(j,i)
    end do
  end do

end subroutine dgemm_fused_sCD_8_self



subroutine dgemm_fused_sCD_1_k(nrows,nvecw,v,w,N,M)
  implicit none
  integer, intent(in) :: nrows, nvecw
  real(kind=8), intent(in) :: v(1,nrows)
  real(kind=8), intent(inout) :: w(nvecw,nrows)
  real(kind=8), intent(in) :: N(nvecw,nvecw)
  real(kind=8), intent(out) :: M(1,nvecw)
  real(kind=8) :: work(nvecw)
  integer :: i, j
!dir$ assume_aligned v:64, w:64, M:64, N:64

  M = 0.
!$omp parallel do reduction(+:M) private(work) schedule(static)
  do i = 1, nrows, 1
    work = w(:,i)
    do j = 1, nvecw, 1
      w(j,i) = sum(N(:,j)*work(:))
      M(:,j) = M(:,j) + v(:,i)*w(j,i)
    end do
  end do

end subroutine dgemm_fused_sCD_1_k

subroutine dgemm_fused_sCD_2_k(nrows,nvecw,v,w,N,M)
  implicit none
  integer, intent(in) :: nrows, nvecw
  real(kind=8), intent(in) :: v(2,nrows)
  real(kind=8), intent(inout) :: w(nvecw,nrows)
  real(kind=8), intent(in) :: N(nvecw,nvecw)
  real(kind=8), intent(out) :: M(2,nvecw)
  real(kind=8) :: work(nvecw)
  integer :: i, j
!dir$ assume_aligned v:64, w:64, M:64, N:64

  M = 0.
!$omp parallel do reduction(+:M) private(work) schedule(static)
  do i = 1, nrows, 1
    work = w(:,i)
    do j = 1, nvecw, 1
      w(j,i) = sum(N(:,j)*work(:))
      M(:,j) = M(:,j) + v(:,i)*w(j,i)
    end do
  end do

end subroutine dgemm_fused_sCD_2_k

subroutine dgemm_fused_sCD_4_k(nrows,nvecw,v,w,N,M)
  implicit none
  integer, intent(in) :: nrows, nvecw
  real(kind=8), intent(in) :: v(4,nrows)
  real(kind=8), intent(inout) :: w(nvecw,nrows)
  real(kind=8), intent(in) :: N(nvecw,nvecw)
  real(kind=8), intent(out) :: M(4,nvecw)
  real(kind=8) :: work(nvecw)
  integer :: i, j
!dir$ assume_aligned v:64, w:64, M:64, N:64

  M = 0.
!$omp parallel do reduction(+:M) private(work) schedule(static)
  do i = 1, nrows, 1
    work = w(:,i)
    do j = 1, nvecw, 1
      w(j,i) = sum(N(:,j)*work(:))
      M(:,j) = M(:,j) + v(:,i)*w(j,i)
    end do
  end do

end subroutine dgemm_fused_sCD_4_k

subroutine dgemm_fused_sCD_8_k(nrows,nvecw,v,w,N,M)
  implicit none
  integer, intent(in) :: nrows, nvecw
  real(kind=8), intent(in) :: v(8,nrows)
  real(kind=8), intent(inout) :: w(nvecw,nrows)
  real(kind=8), intent(in) :: N(nvecw,nvecw)
  real(kind=8), intent(out) :: M(8,nvecw)
  real(kind=8) :: work(nvecw)
  integer :: i, j
!dir$ assume_aligned v:64, w:64, M:64, N:64

  M = 0.
!$omp parallel do reduction(+:M) private(work) schedule(static)
  do i = 1, nrows, 1
    work = w(:,i)
    do j = 1, nvecw, 1
      w(j,i) = sum(N(:,j)*work(:))
      M(:,j) = M(:,j) + v(:,i)*w(j,i)
    end do
  end do

end subroutine dgemm_fused_sCD_8_k


subroutine dgemm_fused_sCD_1_strided_k(nrows,nvecw,v,w,ldw,N,M)
  implicit none
  integer, intent(in) :: nrows, nvecw, ldw
  real(kind=8), intent(in) :: v(1,nrows)
  real(kind=8), intent(inout) :: w(ldw,nrows)
  real(kind=8), intent(in) :: N(nvecw,nvecw)
  real(kind=8), intent(out) :: M(1,nvecw)
  real(kind=8) :: work(nvecw)
  integer :: i, j
!dir$ assume_aligned v:64, w:8, M:64, N:64

  M = 0.
!$omp parallel do reduction(+:M) private(work) schedule(static)
  do i = 1, nrows, 1
    work = w(1:nvecw,i)
    do j = 1, nvecw, 1
      w(j,i) = sum(N(:,j)*work(:))
      M(:,j) = M(:,j) + v(:,i)*w(j,i)
    end do
  end do

end subroutine dgemm_fused_sCD_1_strided_k

subroutine dgemm_fused_sCD_2_strided_k(nrows,nvecw,v,w,ldw,N,M)
  implicit none
  integer, intent(in) :: nrows, nvecw, ldw
  real(kind=8), intent(in) :: v(2,nrows)
  real(kind=8), intent(inout) :: w(ldw,nrows)
  real(kind=8), intent(in) :: N(nvecw,nvecw)
  real(kind=8), intent(out) :: M(2,nvecw)
  real(kind=8) :: work(nvecw)
  integer :: i, j
!dir$ assume_aligned v:64, w:8, M:64, N:64

  M = 0.
!$omp parallel do reduction(+:M) private(work) schedule(static)
  do i = 1, nrows, 1
    work = w(1:nvecw,i)
    do j = 1, nvecw, 1
      w(j,i) = sum(N(:,j)*work(:))
      M(:,j) = M(:,j) + v(:,i)*w(j,i)
    end do
  end do

end subroutine dgemm_fused_sCD_2_strided_k

subroutine dgemm_fused_sCD_4_strided_k(nrows,nvecw,v,w,ldw,N,M)
  implicit none
  integer, intent(in) :: nrows, nvecw, ldw
  real(kind=8), intent(in) :: v(4,nrows)
  real(kind=8), intent(inout) :: w(ldw,nrows)
  real(kind=8), intent(in) :: N(nvecw,nvecw)
  real(kind=8), intent(out) :: M(4,nvecw)
  real(kind=8) :: work(nvecw)
  integer :: i, j
!dir$ assume_aligned v:64, w:8, M:64, N:64

  M = 0.
!$omp parallel do reduction(+:M) private(work) schedule(static)
  do i = 1, nrows, 1
    work = w(1:nvecw,i)
    do j = 1, nvecw, 1
      w(j,i) = sum(N(:,j)*work(:))
      M(:,j) = M(:,j) + v(:,i)*w(j,i)
    end do
  end do

end subroutine dgemm_fused_sCD_4_strided_k

subroutine dgemm_fused_sCD_8_strided_k(nrows,nvecw,v,w,ldw,N,M)
  implicit none
  integer, intent(in) :: nrows, nvecw, ldw
  real(kind=8), intent(in) :: v(8,nrows)
  real(kind=8), intent(inout) :: w(ldw,nrows)
  real(kind=8), intent(in) :: N(nvecw,nvecw)
  real(kind=8), intent(out) :: M(8,nvecw)
  real(kind=8) :: work(nvecw)
  integer :: i, j
!dir$ assume_aligned v:64, w:8, M:64, N:64

  M = 0.
!$omp parallel do reduction(+:M) private(work) schedule(static)
  do i = 1, nrows, 1
    work = w(1:nvecw,i)
    do j = 1, nvecw, 1
      w(j,i) = sum(N(:,j)*work(:))
      M(:,j) = M(:,j) + v(:,i)*w(j,i)
    end do
  end do

end subroutine dgemm_fused_sCD_8_strided_k


subroutine dgemm_fused_sCD_k_1(nrows,nvecv,v,w,N,M)
  implicit none
  integer, intent(in) :: nrows, nvecv
  real(kind=8), intent(in) :: v(nvecv,nrows)
  real(kind=8), intent(inout) :: w(1,nrows)
  real(kind=8), intent(in) :: N(1,1)
  real(kind=8), intent(out) :: M(nvecv,1)
  real(kind=8) :: work(1)
  integer :: i, j
!dir$ assume_aligned v:64, w:64, M:64, N:64

  M = 0.
!$omp parallel do reduction(+:M) private(work) schedule(static)
  do i = 1, nrows, 1
    work = w(:,i)
    do j = 1, 1, 1
      w(j,i) = sum(N(:,j)*work(:))
      M(:,j) = M(:,j) + v(:,i)*w(j,i)
    end do
  end do

end subroutine dgemm_fused_sCD_k_1


subroutine dgemm_fused_sCD_k_2(nrows,nvecv,v,w,N,M)
  implicit none
  integer, intent(in) :: nrows, nvecv
  real(kind=8), intent(in) :: v(nvecv,nrows)
  real(kind=8), intent(inout) :: w(2,nrows)
  real(kind=8), intent(in) :: N(2,2)
  real(kind=8), intent(out) :: M(nvecv,2)
  real(kind=8) :: work(2)
  integer :: i, j
!dir$ assume_aligned v:64, w:64, M:64, N:64

  M = 0.
!$omp parallel do reduction(+:M) private(work) schedule(static)
  do i = 1, nrows, 1
    work = w(:,i)
    do j = 1, 2, 1
      w(j,i) = sum(N(:,j)*work(:))
      M(:,j) = M(:,j) + v(:,i)*w(j,i)
    end do
  end do

end subroutine dgemm_fused_sCD_k_2


subroutine dgemm_fused_sCD_k_4(nrows,nvecv,v,w,N,M)
  implicit none
  integer, intent(in) :: nrows, nvecv
  real(kind=8), intent(in) :: v(nvecv,nrows)
  real(kind=8), intent(inout) :: w(4,nrows)
  real(kind=8), intent(in) :: N(4,4)
  real(kind=8), intent(out) :: M(nvecv,4)
  real(kind=8) :: work(4)
  integer :: i, j
!dir$ assume_aligned v:64, w:64, M:64, N:64

  M = 0.
!$omp parallel do reduction(+:M) private(work) schedule(static)
  do i = 1, nrows, 1
    work = w(:,i)
    do j = 1, 4, 1
      w(j,i) = sum(N(:,j)*work(:))
      M(:,j) = M(:,j) + v(:,i)*w(j,i)
    end do
  end do

end subroutine dgemm_fused_sCD_k_4


subroutine dgemm_fused_sCD_k_8(nrows,nvecv,v,w,N,M)
  implicit none
  integer, intent(in) :: nrows, nvecv
  real(kind=8), intent(in) :: v(nvecv,nrows)
  real(kind=8), intent(inout) :: w(8,nrows)
  real(kind=8), intent(in) :: N(8,8)
  real(kind=8), intent(out) :: M(nvecv,8)
  real(kind=8) :: work(8)
  integer :: i, j
!dir$ assume_aligned v:64, w:64, M:64, N:64

  M = 0.
!$omp parallel do reduction(+:M) private(work) schedule(static)
  do i = 1, nrows, 1
    work = w(:,i)
    do j = 1, 8, 1
      w(j,i) = sum(N(:,j)*work(:))
      M(:,j) = M(:,j) + v(:,i)*w(j,i)
    end do
  end do

end subroutine dgemm_fused_sCD_k_8


subroutine dgemm_fused_sCD_strided_k_1(nrows,nvecv,v,ldv,w,N,M)
  implicit none
  integer, intent(in) :: nrows, nvecv, ldv
  real(kind=8), intent(in) :: v(ldv,nrows)
  real(kind=8), intent(inout) :: w(1,nrows)
  real(kind=8), intent(in) :: N(1,1)
  real(kind=8), intent(out) :: M(nvecv,1)
  real(kind=8) :: work(1)
  integer :: i, j
!dir$ assume_aligned v:8, w:64, M:64, N:64

  M = 0.
!$omp parallel do reduction(+:M) private(work) schedule(static)
  do i = 1, nrows, 1
    work = w(:,i)
    do j = 1, 1, 1
      w(j,i) = sum(N(:,j)*work(:))
      M(:,j) = M(:,j) + v(1:nvecv,i)*w(j,i)
    end do
  end do

end subroutine dgemm_fused_sCD_strided_k_1


subroutine dgemm_fused_sCD_strided_k_2(nrows,nvecv,v,ldv,w,N,M)
  implicit none
  integer, intent(in) :: nrows, nvecv, ldv
  real(kind=8), intent(in) :: v(ldv,nrows)
  real(kind=8), intent(inout) :: w(2,nrows)
  real(kind=8), intent(in) :: N(2,2)
  real(kind=8), intent(out) :: M(nvecv,2)
  real(kind=8) :: work(2)
  integer :: i, j
!dir$ assume_aligned v:8, w:64, M:64, N:64

  M = 0.
!$omp parallel do reduction(+:M) private(work) schedule(static)
  do i = 1, nrows, 1
    work = w(:,i)
    do j = 1, 2, 1
      w(j,i) = sum(N(:,j)*work(:))
      M(:,j) = M(:,j) + v(1:nvecv,i)*w(j,i)
    end do
  end do

end subroutine dgemm_fused_sCD_strided_k_2


subroutine dgemm_fused_sCD_strided_k_4(nrows,nvecv,v,ldv,w,N,M)
  implicit none
  integer, intent(in) :: nrows, nvecv, ldv
  real(kind=8), intent(in) :: v(ldv,nrows)
  real(kind=8), intent(inout) :: w(4,nrows)
  real(kind=8), intent(in) :: N(4,4)
  real(kind=8), intent(out) :: M(nvecv,4)
  real(kind=8) :: work(4)
  integer :: i, j
!dir$ assume_aligned v:8, w:64, M:64, N:64

  M = 0.
!$omp parallel do reduction(+:M) private(work) schedule(static)
  do i = 1, nrows, 1
    work = w(:,i)
    do j = 1, 4, 1
      w(j,i) = sum(N(:,j)*work(:))
      M(:,j) = M(:,j) + v(1:nvecv,i)*w(j,i)
    end do
  end do

end subroutine dgemm_fused_sCD_strided_k_4


subroutine dgemm_fused_sCD_strided_k_8(nrows,nvecv,v,ldv,w,N,M)
  implicit none
  integer, intent(in) :: nrows, nvecv, ldv
  real(kind=8), intent(in) :: v(ldv,nrows)
  real(kind=8), intent(inout) :: w(8,nrows)
  real(kind=8), intent(in) :: N(8,8)
  real(kind=8), intent(out) :: M(nvecv,8)
  real(kind=8) :: work(8)
  integer :: i, j
!dir$ assume_aligned v:8, w:64, M:64, N:64

  M = 0.
!$omp parallel do reduction(+:M) private(work) schedule(static)
  do i = 1, nrows, 1
    work = w(:,i)
    do j = 1, 8, 1
      w(j,i) = sum(N(:,j)*work(:))
      M(:,j) = M(:,j) + v(1:nvecv,i)*w(j,i)
    end do
  end do

end subroutine dgemm_fused_sCD_strided_k_8



subroutine dgemm_fused_sCD_generic(nrows,nvecv,nvecw,v,ldv,w,ldw,N,M)
  implicit none
  integer, intent(in) :: nrows, nvecv, nvecw, ldv, ldw
  real(kind=8), intent(in) :: v(ldv,nrows)
  real(kind=8), intent(inout) :: w(ldw,nrows)
  real(kind=8), intent(in) :: N(nvecw,nvecw)
  real(kind=8), intent(out) :: M(nvecv,nvecw)
  real(kind=8) :: work(nvecw)
  integer :: i, j
!dir$ assume_aligned v:8, w:8, M:64, N:64

  M = 0.
!$omp parallel do reduction(+:M) private(work) schedule(static)
  do i = 1, nrows, 1
    work = w(1:nvecw,i)
    do j = 1, nvecw, 1
      w(j,i) = sum(N(:,j)*work(:))
      M(:,j) = M(:,j) + v(1:nvecv,i)*w(j,i)
    end do
  end do

end subroutine dgemm_fused_sCD_generic


subroutine dgemm_fused_sCD_generic_self(nrows,nvecw,w,ldw,N,M)
  implicit none
  integer, intent(in) :: nrows, nvecw, ldw
  real(kind=8), intent(inout) :: w(ldw,nrows)
  real(kind=8), intent(in) :: N(nvecw,nvecw)
  real(kind=8), intent(out) :: M(nvecw,nvecw)
  real(kind=8) :: work(nvecw)
  integer :: i, j
!dir$ assume_aligned w:8, M:64, N:64

  M = 0.
!$omp parallel do reduction(+:M) private(work) schedule(static)
  do i = 1, nrows, 1
    work = w(1:nvecw,i)
    do j = 1, nvecw, 1
      w(j,i) = sum(N(:,j)*work(:))
      M(:j,j) = M(:j,j) + w(:j,i)*w(j,i)
    end do
  end do

  do i = 1, nvecw, 1
    do j = i+1, nvecw, 1
      M(j,i) = M(i,j)
    end do
  end do


end subroutine dgemm_fused_sCD_generic_self


