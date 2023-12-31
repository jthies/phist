/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/

!> \file spmv_buff_cpy_kernels.f90
!! Fast parallel subroutines for different block sizes that copy data from a mvec to an MPI communication
!! buffer for crsMat_module
!! \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
!!

subroutine dspmv_buff_cpy_1(nlocal, nbuff, val, ind, buff)
  implicit none
  integer, intent(in) :: nlocal, nbuff
  real(kind=8), intent(in) :: val(nlocal)
  integer, intent(in) :: ind(nbuff,2)
  real(kind=8), intent(out) :: buff(nbuff)
  integer :: i
!dir$ assume_aligned val:64, ind:64, buff:64

!$omp parallel do schedule(static)
  do i = 1, nbuff, 1
    buff(ind(i,2)) = val(ind(i,1))
  end do

end subroutine dspmv_buff_cpy_1


subroutine dspmv_buff_cpy_2(nlocal, nbuff, val, ind, buff)
  implicit none
  integer, intent(in) :: nlocal, nbuff
  real(kind=8), intent(in) :: val(2,nlocal)
  integer, intent(in) :: ind(nbuff,2)
  real(kind=8), intent(out) :: buff(2,nbuff)
  integer :: i
!dir$ assume_aligned val:64, ind:64, buff:64

!$omp parallel do schedule(static)
  do i = 1, nbuff, 1
    buff(:,ind(i,2)) = val(:,ind(i,1))
  end do

end subroutine dspmv_buff_cpy_2


subroutine dspmv_buff_cpy_4(nlocal, nbuff, val, ind, buff)
  implicit none
  integer, intent(in) :: nlocal, nbuff
  real(kind=8), intent(in) :: val(4,nlocal)
  integer, intent(in) :: ind(nbuff,2)
  real(kind=8), intent(out) :: buff(4,nbuff)
  integer :: i
!dir$ assume_aligned val:64, ind:64, buff:64

!$omp parallel do schedule(static)
  do i = 1, nbuff, 1
    buff(:,ind(i,2)) = val(:,ind(i,1))
  end do

end subroutine dspmv_buff_cpy_4


subroutine dspmv_buff_cpy_8(nlocal, nbuff, val, ind, buff)
  implicit none
  integer, intent(in) :: nlocal, nbuff
  real(kind=8), intent(in) :: val(8,nlocal)
  integer, intent(in) :: ind(nbuff,2)
  real(kind=8), intent(out) :: buff(8,nbuff)
  integer :: i
!dir$ assume_aligned val:64, ind:64, buff:64

!$omp parallel do schedule(static)
  do i = 1, nbuff, 1
    buff(:,ind(i,2)) = val(:,ind(i,1))
  end do

end subroutine dspmv_buff_cpy_8


subroutine dspmv_buff_cpy_strided_1(nlocal, nbuff, ldx, val, ind, buff)
  implicit none
  integer, intent(in) :: nlocal, nbuff, ldx
  real(kind=8), intent(in) :: val(ldx,nlocal)
  integer, intent(in) :: ind(nbuff,2)
  real(kind=8), intent(out) :: buff(nbuff)
  integer :: i

!$omp parallel do schedule(static)
  do i = 1, nbuff, 1
    buff(ind(i,2)) = val(1,ind(i,1))
  end do

end subroutine dspmv_buff_cpy_strided_1


subroutine dspmv_buff_cpy_strided_2(nlocal, nbuff, ldx, val, ind, buff)
  implicit none
  integer, intent(in) :: nlocal, nbuff, ldx
  real(kind=8), intent(in) :: val(ldx,nlocal)
  integer, intent(in) :: ind(nbuff,2)
  real(kind=8), intent(out) :: buff(2,nbuff)
  integer :: i
!dir$ assume_aligned val:8, ind:64, buff:64

!$omp parallel do schedule(static)
  do i = 1, nbuff, 1
    buff(:,ind(i,2)) = val(1:2,ind(i,1))
  end do

end subroutine dspmv_buff_cpy_strided_2


subroutine dspmv_buff_cpy_strided_4(nlocal, nbuff, ldx, val, ind, buff)
  implicit none
  integer, intent(in) :: nlocal, nbuff, ldx
  real(kind=8), intent(in) :: val(ldx,nlocal)
  integer, intent(in) :: ind(nbuff,2)
  real(kind=8), intent(out) :: buff(4,nbuff)
  integer :: i
!dir$ assume_aligned val:8, ind:64, buff:64

!$omp parallel do schedule(static)
  do i = 1, nbuff, 1
    buff(:,ind(i,2)) = val(1:4,ind(i,1))
  end do

end subroutine dspmv_buff_cpy_strided_4



subroutine dspmv_buff_cpy_strided_8(nlocal, nbuff, ldx, val, ind, buff)
  implicit none
  integer, intent(in) :: nlocal, nbuff, ldx
  real(kind=8), intent(in) :: val(ldx,nlocal)
  integer, intent(in) :: ind(nbuff,2)
  real(kind=8), intent(out) :: buff(8,nbuff)
  integer :: i
!dir$ assume_aligned val:8, ind:64, buff:64

!$omp parallel do schedule(static)
  do i = 1, nbuff, 1
    buff(:,ind(i,2)) = val(1:8,ind(i,1))
  end do

end subroutine dspmv_buff_cpy_strided_8


subroutine dspmv_buff_cpy_general(nvec, nlocal, nbuff, ldx, val, ind, buff)
  implicit none
  integer, intent(in) :: nvec, nlocal, nbuff, ldx
  real(kind=8), intent(in) :: val(ldx,nlocal)
  integer, intent(in) :: ind(nbuff,2)
  real(kind=8), intent(out) :: buff(nvec,nbuff)
  integer :: i
!dir$ assume_aligned val:8, ind:64, buff:64

!$omp parallel do schedule(static)
  do i = 1, nbuff, 1
    buff(:,ind(i,2)) = val(1:nvec,ind(i,1))
  end do

end subroutine dspmv_buff_cpy_general


