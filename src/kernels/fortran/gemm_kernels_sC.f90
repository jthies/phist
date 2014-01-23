! gemm kernels for mecT_times_mvec

subroutine dgemm_sC_1(nrows,nvecw,v,w,M)
  implicit none
  integer, intent(in) :: nrows, nvecw
  real(kind=8), intent(in) :: v(1,nrows)
  real(kind=8), intent(in) :: w(nvecw,nrows)
  real(kind=8), intent(out) :: M(1,nvecw)
  integer :: i, j

  M = 0.
!$omp parallel do reduction(+:M) schedule(static)
  do i = 1, nrows
    do j = 1, nvecw
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

  M = 0.
!$omp parallel do reduction(+:M) schedule(static)
  do i = 1, nrows
    do j = 1, nvecw
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

  M = 0.
!$omp parallel do reduction(+:M) schedule(static)
  do i = 1, nrows
    do j = 1, nvecw
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

  M = 0.
!$omp parallel do reduction(+:M) schedule(static)
  do i = 1, nrows
    do j = 1, nvecw
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

  M = 0.
!$omp parallel do reduction(+:M) schedule(static)
  do i = 1, nrows
    do j = 1, nvecw
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

  M = 0.
!$omp parallel do reduction(+:M) schedule(static)
  do i = 1, nrows
    do j = 1, nvecw
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

  M = 0.
!$omp parallel do reduction(+:M) schedule(static)
  do i = 1, nrows
    do j = 1, nvecw
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

  M = 0.
!$omp parallel do reduction(+:M) schedule(static)
  do i = 1, nrows
    do j = 1, nvecw
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

  M = 0.
!$omp parallel do reduction(+:M) schedule(static)
  do i = 1, nrows
    do j = 1, nvecw
      M(:,j) = M(:,j) + v(1:nvecv,i)*w(j,i)
    end do
  end do

end subroutine dgemm_sC_generic


