/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/

!> \file gather_scatter_kernels.f90
!! Fast parallel kernels to copy between different mvecs for mvec_module
!! \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
!!

subroutine dgather_1(nrows, v, w, ldw)
  integer, intent(in) :: nrows, ldw
  real(kind=8), intent(out) :: v(nrows)
  real(kind=8), intent(in) :: w(ldw,*)
  integer :: i
!dir$ assume_aligned v:64, w:8

!$omp parallel do schedule(static)
  do i = 1, nrows, 1
    v(i) = w(1,i)
  end do

end subroutine dgather_1

subroutine dgather_2(nrows, v, w1, ldw1, w2, ldw2)
  integer, intent(in) :: nrows, ldw1, ldw2
  real(kind=8), intent(out) :: v(2,nrows)
  real(kind=8), intent(in) :: w1(ldw1,*)
  real(kind=8), intent(in) :: w2(ldw2,*)
  integer :: i
!dir$ assume_aligned v:64, w1:8, w2:8

!$omp parallel do schedule(static)
  do i = 1, nrows, 1
    v(:,i) = (/w1(1,i),w2(1,i)/)
  end do

end subroutine dgather_2

subroutine dgather_4(nrows, v, w1, ldw1, w2, ldw2, w3, ldw3, w4, ldw4)
  integer, intent(in) :: nrows, ldw1, ldw2
  real(kind=8), intent(out) :: v(4,nrows)
  real(kind=8), intent(in) :: w1(ldw1,*)
  real(kind=8), intent(in) :: w2(ldw2,*)
  real(kind=8), intent(in) :: w3(ldw3,*)
  real(kind=8), intent(in) :: w4(ldw4,*)
  integer :: i
!dir$ assume_aligned v:64, w1:8, w2:8, w3:8, w4:8

!$omp parallel do schedule(static)
  do i = 1, nrows, 1
    v(:,i) = (/w1(1,i),w2(1,i),w3(1,i),w4(1,i)/)
  end do

end subroutine dgather_4


subroutine dscatter_1(nrows, v, w, ldw)
  integer, intent(in) :: nrows, ldw
  real(kind=8), intent(in) :: v(nrows)
  real(kind=8), intent(inout) :: w(ldw,*)
  integer :: i
!dir$ assume_aligned v:64, w:8

!$omp parallel do schedule(static)
  do i = 1, nrows, 1
    w(1,i) = v(i)
  end do

end subroutine dscatter_1

subroutine dscatter_2(nrows, v, w1, ldw1, w2, ldw2)
  integer, intent(in) :: nrows, ldw1, ldw2
  real(kind=8), intent(in) :: v(2,nrows)
  real(kind=8), intent(inout) :: w1(ldw1,*)
  real(kind=8), intent(inout) :: w2(ldw2,*)
  integer :: i
!dir$ assume_aligned v:64, w1:8, w2:8

!$omp parallel do schedule(static)
  do i = 1, nrows, 1
    w1(1,i) = v(1,i)
    w2(1,i) = v(2,i)
  end do

end subroutine dscatter_2

subroutine dscatter_4(nrows, v, w1, ldw1, w2, ldw2, w3, ldw3, w4, ldw4)
  integer, intent(in) :: nrows, ldw1, ldw2
  real(kind=8), intent(in) :: v(4,nrows)
  real(kind=8), intent(inout) :: w1(ldw1,*)
  real(kind=8), intent(inout) :: w2(ldw2,*)
  real(kind=8), intent(inout) :: w3(ldw3,*)
  real(kind=8), intent(inout) :: w4(ldw4,*)
  integer :: i
!dir$ assume_aligned v:64, w1:8, w2:8, w3:8, w4:8

!$omp parallel do schedule(static)
  do i = 1, nrows, 1
    w1(1,i) = v(1,i)
    w2(1,i) = v(2,i)
    w3(1,i) = v(3,i)
    w4(1,i) = v(4,i)
  end do

end subroutine dscatter_4

