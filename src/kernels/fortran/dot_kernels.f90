! call blas
subroutine ddot_strided_1(nrows, v, ldv, w, ldw, vdot)
  implicit none
  integer, intent(in)       :: nrows, ldv, ldw
  real(kind=8), intent(in)  :: v(ldv,*), w(ldw,*)
  real(kind=8), intent(out) :: vdot
  integer :: i

  ! ddot blas interface
  interface
    function ddot(n,v,ldv,w,ldw)
      real(kind=8) :: ddot
      integer, intent(in) :: n, ldv, ldw
      real(kind=8), intent(in) :: v(ldv,*), w(ldw,*)
    end function ddot
  end interface

  vdot = ddot(nrows,v,ldv,w,ldw)

end subroutine ddot_strided_1


!--------------------------------------------------------------------------------
! hopefully fast dot variants for nvec > 1

subroutine ddot_strided_2(nrows, v, ldv, w, ldw, vdot)
  implicit none
  integer, intent(in)       :: nrows, ldv, ldw
  real(kind=8), intent(in)  :: v(ldv,*), w(ldw,*)
  real(kind=8), intent(out) :: vdot(2)
  integer :: i

  vdot = 0.
!$omp parallel do reduction(+:vdot)
  do i = 1, nrows
    vdot = vdot + v(1:2,i)*w(1:2,i)
  end do
end subroutine ddot_strided_2


subroutine ddot_strided_4(nrows, v, ldv, w, ldw, vdot)
  implicit none
  integer, intent(in)       :: nrows, ldv, ldw
  real(kind=8), intent(in)  :: v(ldv,*), w(ldw,*)
  real(kind=8), intent(out) :: vdot(4)
  integer :: i

  vdot = 0.
!$omp parallel do reduction(+:vdot)
  do i = 1, nrows
    vdot = vdot + v(1:4,i)*w(1:4,i)
  end do
end subroutine ddot_strided_4


subroutine ddot_strided_8(nrows, v, ldv, w, ldw, vdot)
  implicit none
  integer, intent(in)       :: nrows, ldv, ldw
  real(kind=8), intent(in)  :: v(ldv,*), w(ldw,*)
  real(kind=8), intent(out) :: vdot(8)
  integer :: i

  vdot = 0.
!$omp parallel do reduction(+:vdot)
  do i = 1, nrows
    vdot = vdot + v(1:8,i)*w(1:8,i)
  end do
end subroutine ddot_strided_8


subroutine ddot_2(nrows, v, w, vdot)
  implicit none
  integer, intent(in)       :: nrows
  real(kind=8), intent(in)  :: v(2,nrows), w(2,nrows)
  real(kind=8), intent(out) :: vdot(2)
  integer :: i

  vdot = 0.
!$omp parallel do reduction(+:vdot)
  do i = 1, nrows
    vdot = vdot + v(:,i)*w(:,i)
  end do
end subroutine ddot_2


subroutine ddot_4(nrows, v, w, vdot)
  implicit none
  integer, intent(in)       :: nrows
  real(kind=8), intent(in)  :: v(4,nrows), w(4,nrows)
  real(kind=8), intent(out) :: vdot(4)
  integer :: i

  vdot = 0.
!$omp parallel do reduction(+:vdot)
  do i = 1, nrows
    vdot = vdot + v(:,i)*w(:,i)
  end do
end subroutine ddot_4


subroutine ddot_8(nrows, v, w, vdot)
  implicit none
  integer, intent(in)       :: nrows
  real(kind=8), intent(in)  :: v(8,nrows), w(8,nrows)
  real(kind=8), intent(out) :: vdot(8)
  integer :: i

  vdot = 0.
!$omp parallel do reduction(+:vdot)
  do i = 1, nrows
    vdot = vdot + v(:,i)*w(:,i)
  end do
end subroutine ddot_8


subroutine ddot_general(nrows, nvec, v, ldv, w, ldw, vdot)
  implicit none
  integer, intent(in)       :: nrows, nvec, ldv, ldw
  real(kind=8), intent(in)  :: v(ldv,*), w(ldw,*)
  real(kind=8), intent(out) :: vdot(nvec)
  integer :: i

  vdot = 0.
!$omp parallel do reduction(+:vdot)
  do i = 1, nrows
    vdot(1:nvec) = vdot(1:nvec) + v(1:nvec,i)*w(1:nvec,i)
  end do
end subroutine ddot_general
