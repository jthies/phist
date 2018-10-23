/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/

!> \file dot_kernels.f90
!! Fast parallel BLAS-dot like routines for different blocksizes for mvec_module
!! \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>"
!!

subroutine ddot_strided_1(nrows, v, ldv, w, ldw, vdot)
  implicit none
  integer, intent(in)       :: nrows, ldv, ldw
  real(kind=8), intent(in)  :: v(ldv,*), w(ldw,*)
  real(kind=8), intent(out) :: vdot
  integer :: i
!dir$ assume_aligned v:8, w:8

  vdot = 0.
!$omp parallel do reduction(+:vdot) schedule(static)
  do i = 1, nrows, 1
    vdot = vdot + v(1,i)*w(1,i)
  end do
end subroutine ddot_strided_1


subroutine ddot_strided_2(nrows, v, ldv, w, ldw, vdot)
  implicit none
  integer, intent(in)       :: nrows, ldv, ldw
  real(kind=8), intent(in)  :: v(ldv,*), w(ldw,*)
  real(kind=8), intent(out) :: vdot(2)
  integer :: i
!dir$ assume_aligned v:8, w:8

  vdot = 0.
!$omp parallel do reduction(+:vdot) schedule(static)
  do i = 1, nrows, 1
    vdot = vdot + v(1:2,i)*w(1:2,i)
  end do
end subroutine ddot_strided_2


subroutine ddot_strided_4(nrows, v, ldv, w, ldw, vdot)
  implicit none
  integer, intent(in)       :: nrows, ldv, ldw
  real(kind=8), intent(in)  :: v(ldv,*), w(ldw,*)
  real(kind=8), intent(out) :: vdot(4)
  integer :: i
!dir$ assume_aligned v:8, w:8

  vdot = 0.
!$omp parallel do reduction(+:vdot) schedule(static)
  do i = 1, nrows, 1
    vdot = vdot + v(1:4,i)*w(1:4,i)
  end do
end subroutine ddot_strided_4


subroutine ddot_strided_8(nrows, v, ldv, w, ldw, vdot)
  implicit none
  integer, intent(in)       :: nrows, ldv, ldw
  real(kind=8), intent(in)  :: v(ldv,*), w(ldw,*)
  real(kind=8), intent(out) :: vdot(8)
  integer :: i
!dir$ assume_aligned v:8, w:8

  vdot = 0.
!$omp parallel do reduction(+:vdot) schedule(static)
  do i = 1, nrows, 1
    vdot = vdot + v(1:8,i)*w(1:8,i)
  end do
end subroutine ddot_strided_8


subroutine ddot_1(nrows, v, w, vdot)
  implicit none
  integer, intent(in)       :: nrows
  real(kind=8), intent(in)  :: v(nrows), w(nrows)
  real(kind=8), intent(out) :: vdot
  integer :: i
!dir$ assume_aligned v:64, w:64

  vdot = 0.
!$omp parallel do reduction(+:vdot) schedule(static)
  do i = 1, nrows, 1
    vdot = vdot + v(i)*w(i)
  end do
end subroutine ddot_1


subroutine ddot_2(nrows, v, w, vdot)
  implicit none
  integer, intent(in)       :: nrows
  real(kind=8), intent(in)  :: v(2,nrows), w(2,nrows)
  real(kind=8), intent(out) :: vdot(2)
  integer :: i
!dir$ assume_aligned v:64, w:64

  vdot = 0.
!$omp parallel do reduction(+:vdot) schedule(static)
  do i = 1, nrows, 1
    vdot = vdot + v(:,i)*w(:,i)
  end do
end subroutine ddot_2


subroutine ddot_4(nrows, v, w, vdot)
  implicit none
  integer, intent(in)       :: nrows
  real(kind=8), intent(in)  :: v(4,nrows), w(4,nrows)
  real(kind=8), intent(out) :: vdot(4)
  integer :: i
!dir$ assume_aligned v:64, w:64

  vdot = 0.
!$omp parallel do reduction(+:vdot) schedule(static)
  do i = 1, nrows, 1
    vdot = vdot + v(:,i)*w(:,i)
  end do
end subroutine ddot_4


subroutine ddot_8(nrows, v, w, vdot)
  implicit none
  integer, intent(in)       :: nrows
  real(kind=8), intent(in)  :: v(8,nrows), w(8,nrows)
  real(kind=8), intent(out) :: vdot(8)
  integer :: i
!dir$ assume_aligned v:64, w:64

  vdot = 0.
!$omp parallel do reduction(+:vdot) schedule(static)
  do i = 1, nrows, 1
    vdot = vdot + v(:,i)*w(:,i)
  end do
end subroutine ddot_8


subroutine ddot_general(nrows, nvec, v, ldv, w, ldw, vdot)
  implicit none
  integer, intent(in)       :: nrows, nvec, ldv, ldw
  real(kind=8), intent(in)  :: v(ldv,*), w(ldw,*)
  real(kind=8), intent(out) :: vdot(nvec)
  integer :: i
!dir$ assume_aligned v:8, w:8

  vdot = 0.
!$omp parallel do reduction(+:vdot) schedule(static)
  do i = 1, nrows, 1
    vdot(1:nvec) = vdot(1:nvec) + v(1:nvec,i)*w(1:nvec,i)
  end do
end subroutine ddot_general
