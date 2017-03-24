/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
!> \file gemm_kernels_sB_add_sD.f90
!! Fast parallel BLAS-gemm like subroutines augmented with the product of the result for different blocksizes for mvec_times_sdmat_add_mvec_times_sdmat in mvec_module
!! \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
!!

subroutine dgemm_sB_add_sD_1_k(nrows,k, A, B, C, D)
  implicit none
  integer, intent(in) :: nrows, k
  real(kind=8), intent(in) :: A(k,nrows)
  real(kind=8), intent(in) :: B(1,k)
  real(kind=8), intent(inout) :: C(1,nrows)
  real(kind=8), intent(in) :: D(1,1)
  integer :: i, j
!dir$ assume_aligned A:64, B:64, C:64, D:64

!$omp parallel do schedule(static)
  do i = 1, nrows, 1
    C(1,i) = D(1,1)*C(1,i)
    do j = 1, k, 1
      C(1,i) = C(1,i) + B(1,j)*A(j,i)
    end do
  end do

end subroutine dgemm_sB_add_sD_1_k

subroutine dgemm_sB_add_sD_2_k(nrows,k, A, B, C, D)
  implicit none
  integer, intent(in) :: nrows, k
  real(kind=8), intent(in) :: A(k,nrows)
  real(kind=8), intent(in) :: B(2,k)
  real(kind=8), intent(inout) :: C(2,nrows)
  real(kind=8), intent(in) :: D(2,2)
  integer :: i, j
  real(kind=8) :: work(2)
!dir$ assume_aligned A:64, B:64, C:64, D:64

!$omp parallel do private(work) schedule(static)
  do i = 1, nrows, 1
    work = 0.
    do j = 1, k, 1
      work = work + B(:,j)*A(j,i)
    end do
    do j = 1, 2, 1
      work = work + D(:,j)*C(j,i)
    end do
    C(:,i) = work
  end do

end subroutine dgemm_sB_add_sD_2_k


subroutine dgemm_sB_add_sD_4_k(nrows,k, A, B, C, D)
  implicit none
  integer, intent(in) :: nrows, k
  real(kind=8), intent(in) :: A(k,nrows)
  real(kind=8), intent(in) :: B(4,k)
  real(kind=8), intent(inout) :: C(4,nrows)
  real(kind=8), intent(in) :: D(4,4)
  integer :: i, j
  real(kind=8) :: work(4)
!dir$ assume_aligned A:64, B:64, C:64, D:64

!$omp parallel do private(work) schedule(static)
  do i = 1, nrows, 1
    work = 0.
    do j = 1, k, 1
      work = work + B(:,j)*A(j,i)
    end do
    do j = 1, 4, 1
      work = work + D(:,j)*C(j,i)
    end do
    C(:,i) = work
  end do

end subroutine dgemm_sB_add_sD_4_k

subroutine dgemm_sB_add_sD_8_k(nrows,k, A, B, C, D)
  implicit none
  integer, intent(in) :: nrows, k
  real(kind=8), intent(in) :: A(k,nrows)
  real(kind=8), intent(in) :: B(8,k)
  real(kind=8), intent(inout) :: C(8,nrows)
  real(kind=8), intent(in) :: D(8,8)
  integer :: i, j
  real(kind=8) :: work(8)
!dir$ assume_aligned A:64, B:64, C:64, D:64

!$omp parallel do private(work) schedule(static)
  do i = 1, nrows, 1
    work = 0.
    do j = 1, k, 1
      work = work + B(:,j)*A(j,i)
    end do
    do j = 1, 8, 1
      work = work + D(:,j)*C(j,i)
    end do
    C(:,i) = work
  end do

end subroutine dgemm_sB_add_sD_8_k



subroutine dgemm_sB_add_sD_1_strided_k(nrows,k, A, lda, B, C, D)
  implicit none
  integer, intent(in) :: nrows, k, lda
  real(kind=8), intent(in) :: A(lda,*)
  real(kind=8), intent(in) :: B(1,k)
  real(kind=8), intent(inout) :: C(1,nrows)
  real(kind=8), intent(in) :: D(1,1)
  integer :: i, j
!dir$ assume_aligned A:8, B:64, C:64, D:64

!$omp parallel do schedule(static)
  do i = 1, nrows, 1
    C(1,i) = D(1,1)*C(1,i)
    do j = 1, k, 1
      C(1,i) = C(1,i) + B(1,j)*A(j,i)
    end do
  end do

end subroutine dgemm_sB_add_sD_1_strided_k

subroutine dgemm_sB_add_sD_2_strided_k(nrows,k, A, lda, B, C, D)
  implicit none
  integer, intent(in) :: nrows, k, lda
  real(kind=8), intent(in) :: A(lda,*)
  real(kind=8), intent(in) :: B(2,k)
  real(kind=8), intent(inout) :: C(2,nrows)
  real(kind=8), intent(in) :: D(2,2)
  integer :: i, j
  real(kind=8) :: work(2)
!dir$ assume_aligned A:64, B:64, C:64, D:64

!$omp parallel do private(work) schedule(static)
  do i = 1, nrows, 1
    work = 0.
    do j = 1, k, 1
      work = work + B(:,j)*A(j,i)
    end do
    do j = 1, 2, 1
      work = work + D(:,j)*C(j,i)
    end do
    C(:,i) = work
  end do

end subroutine dgemm_sB_add_sD_2_strided_k

subroutine dgemm_sB_add_sD_4_strided_k(nrows,k, A, lda, B, C, D)
  implicit none
  integer, intent(in) :: nrows, k, lda
  real(kind=8), intent(in) :: A(lda,*)
  real(kind=8), intent(in) :: B(4,k)
  real(kind=8), intent(inout) :: C(4,nrows)
  real(kind=8), intent(in) :: D(4,4)
  integer :: i, j
  real(kind=8) :: work(4)
!dir$ assume_aligned A:64, B:64, C:64, D:64

!$omp parallel do private(work) schedule(static)
  do i = 1, nrows, 1
    work = 0.
    do j = 1, k, 1
      work = work + B(:,j)*A(j,i)
    end do
    do j = 1, 4, 1
      work = work + D(:,j)*C(j,i)
    end do
    C(:,i) = work
  end do

end subroutine dgemm_sB_add_sD_4_strided_k

subroutine dgemm_sB_add_sD_8_strided_k(nrows,k, A, lda, B, C, D)
  implicit none
  integer, intent(in) :: nrows, k, lda
  real(kind=8), intent(in) :: A(lda,*)
  real(kind=8), intent(in) :: B(8,k)
  real(kind=8), intent(inout) :: C(8,nrows)
  real(kind=8), intent(in) :: D(8,8)
  integer :: i, j
  real(kind=8) :: work(8)
!dir$ assume_aligned A:64, B:64, C:64, D:64

!$omp parallel do private(work) schedule(static)
  do i = 1, nrows, 1
    work = 0.
    do j = 1, k, 1
      work = work + B(:,j)*A(j,i)
    end do
    do j = 1, 8, 1
      work = work + D(:,j)*C(j,i)
    end do
    C(:,i) = work
  end do

end subroutine dgemm_sB_add_sD_8_strided_k



subroutine dgemm_sB_add_sD_generic(m,n,k,A,lda,B,C,ldc, D)
  implicit none
  integer, intent(in) :: m,n,k,lda,ldc
  real(kind=8), intent(in) :: A(lda,*)
  real(kind=8), intent(in) :: B(k,n)
  real(kind=8), intent(inout) :: C(ldc,*)
  real(kind=8), intent(in) :: D(n,n)
  integer :: i, j
  real(kind=8) :: work(n)
!dir$ assume_aligned A:8, B:64, C:8

!$omp parallel do private(work) schedule(static)
  do i = 1, m, 1
    work = 0.
    do j = 1, n, 1
      work(j) = work(j) + sum(B(:,j)*A(1:k,i))
    end do
    do j = 1, n, 1
      work = work + D(:,j)*C(j,i)
    end do
    C(1:n,i) = work
  end do

end subroutine dgemm_sB_add_sD_generic

