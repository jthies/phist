!--------------------------------------------------------------------------------
! hopefully fast dnrm2 variants for nvec >= 1

subroutine dnrm2_strided_1(nrows, vec, ldv, vnrm)
  implicit none
  integer, intent(in)       :: nrows, ldv
  real(kind=8), intent(in)  :: vec(ldv,*)
  real(kind=8), intent(out) :: vnrm
  integer :: i

  vnrm = 0.
!$omp parallel do reduction(+:vnrm) schedule(static)
  do i = 1, nrows
    vnrm = vnrm + vec(1,i)**2
  end do
  vnrm = sqrt(vnrm)
end subroutine dnrm2_strided_1


subroutine dnrm2_strided_2(nrows, vec, ldv, vnrm)
  implicit none
  integer, intent(in)       :: nrows, ldv
  real(kind=8), intent(in)  :: vec(ldv,*)
  real(kind=8), intent(out) :: vnrm(2)
  integer :: i

  vnrm = 0.
!$omp parallel do reduction(+:vnrm) schedule(static)
  do i = 1, nrows
    vnrm = vnrm + vec(1:2,i)**2
  end do
  vnrm = sqrt(vnrm)
end subroutine dnrm2_strided_2


subroutine dnrm2_strided_4(nrows, vec, ldv, vnrm)
  implicit none
  integer, intent(in)       :: nrows, ldv
  real(kind=8), intent(in)  :: vec(ldv,*)
  real(kind=8), intent(out) :: vnrm(4)
  integer :: i

  vnrm = 0.
!$omp parallel do reduction(+:vnrm) schedule(static)
  do i = 1, nrows
    vnrm = vnrm + vec(1:4,i)**2
  end do
  vnrm = sqrt(vnrm)
end subroutine dnrm2_strided_4


subroutine dnrm2_strided_8(nrows, vec, ldv, vnrm)
  implicit none
  integer, intent(in)       :: nrows, ldv
  real(kind=8), intent(in)  :: vec(ldv,*)
  real(kind=8), intent(out) :: vnrm(8)
  integer :: i

  vnrm = 0.
!$omp parallel do reduction(+:vnrm) schedule(static)
  do i = 1, nrows
    vnrm = vnrm + vec(1:8,i)**2
  end do
  vnrm = sqrt(vnrm)
end subroutine dnrm2_strided_8


subroutine dnrm2_1(nrows, vec, vnrm)
  implicit none
  integer, intent(in)       :: nrows
  real(kind=8), intent(in)  :: vec(nrows)
  real(kind=8), intent(out) :: vnrm
  integer :: i

  vnrm = 0.
!$omp parallel do reduction(+:vnrm) schedule(static)
  do i = 1, nrows
    vnrm = vnrm + vec(i)**2
  end do
  vnrm = sqrt(vnrm)
end subroutine dnrm2_1


subroutine dnrm2_2(nrows, vec, vnrm)
  implicit none
  integer, intent(in)       :: nrows
  real(kind=8), intent(in)  :: vec(2,nrows)
  real(kind=8), intent(out) :: vnrm(2)
  integer :: i

  vnrm = 0.
!$omp parallel do reduction(+:vnrm) schedule(static)
  do i = 1, nrows
    vnrm = vnrm + vec(:,i)**2
  end do
  vnrm = sqrt(vnrm)
end subroutine dnrm2_2


subroutine dnrm2_4(nrows, vec, vnrm)
  implicit none
  integer, intent(in)       :: nrows
  real(kind=8), intent(in)  :: vec(4,nrows)
  real(kind=8), intent(out) :: vnrm(4)
  integer :: i

  vnrm = 0.
!$omp parallel do reduction(+:vnrm) schedule(static)
  do i = 1, nrows
    vnrm = vnrm + vec(:,i)**2
  end do
  vnrm = sqrt(vnrm)
end subroutine dnrm2_4


subroutine dnrm2_8(nrows, vec, vnrm)
  implicit none
  integer, intent(in)       :: nrows
  real(kind=8), intent(in)  :: vec(8,nrows)
  real(kind=8), intent(out) :: vnrm(8)
  integer :: i

  vnrm = 0.
!$omp parallel do reduction(+:vnrm) schedule(static)
  do i = 1, nrows
    vnrm = vnrm + vec(:,i)**2
  end do
  vnrm = sqrt(vnrm)
end subroutine dnrm2_8


subroutine dnrm2_general(nrows, nvec, vec, ldv, vnrm)
  implicit none
  integer, intent(in)       :: nrows, nvec, ldv
  real(kind=8), intent(in)  :: vec(ldv,*)
  real(kind=8), intent(out) :: vnrm(nvec)
  integer :: i

  vnrm = 0.
!$omp parallel do reduction(+:vnrm) schedule(static)
  do i = 1, nrows
    vnrm(1:nvec) = vnrm(1:nvec) + vec(1:nvec,i)**2
  end do
  vnrm(1:nvec) = sqrt(vnrm(1:nvec))
end subroutine dnrm2_general
