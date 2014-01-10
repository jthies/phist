subroutine dcopy_general(nrows, nvec, x, ldx, y, ldy)
  implicit none
  integer, intent(in) :: nrows, nvec, ldx, ldy
  real(kind=8), intent(in) :: x(ldx,*)
  real(kind=8), intent(out) :: y(ldy,*)
  integer i

  interface
    subroutine dcopy_general_nt_c(nrows,nvec,x,ldx,y,ldy) bind(C)
      use, intrinsic :: iso_c_binding
      integer(C_INT), value :: nrows,nvec,ldx,ldy
      real(C_DOUBLE), intent(in) :: x(*)
      real(C_DOUBLE), intent(out) :: y(*)
    end subroutine dcopy_general_nt_c
  end interface


  if( mod(loc(y),16) .eq. 0 .and. mod(ldy,2) .eq. 0 .and. mod(nvec,2) .eq. 0 .and. mod(ldx,2) .eq. 0 ) then
    call dcopy_general_nt_c(nrows, nvec, x, ldx, y, ldy)
  else
!$omp parallel do
    do i = 1, nrows
      y(1:nvec,i) = x(1:nvec,i)
    end do
  end if
end subroutine dcopy_general


subroutine dscal_2(nrows, alpha, x)
  implicit none
  integer, intent(in) :: nrows
  real(kind=8), intent(in) :: alpha(2)
  real(kind=8), intent(inout) :: x(2,*)
  integer :: i

!$omp parallel do
  do i = 1, nrows
    x(:,i) = alpha(:)*x(:,i)
  end do
end subroutine dscal_2

subroutine dscal_4(nrows, alpha, x)
  implicit none
  integer, intent(in) :: nrows
  real(kind=8), intent(in) :: alpha(4)
  real(kind=8), intent(inout) :: x(4,*)
  integer :: i

!$omp parallel do
  do i = 1, nrows
    x(:,i) = alpha(:)*x(:,i)
  end do
end subroutine dscal_4

subroutine dscal_8(nrows, alpha, x)
  implicit none
  integer, intent(in) :: nrows
  real(kind=8), intent(in) :: alpha(8)
  real(kind=8), intent(inout) :: x(8,*)
  integer :: i

!$omp parallel do
  do i = 1, nrows
    x(:,i) = alpha(:)*x(:,i)
  end do
end subroutine dscal_8


subroutine dscal_strided_2(nrows, alpha, x, ldx)
  implicit none
  integer, intent(in) :: nrows, ldx
  real(kind=8), intent(in) :: alpha(2)
  real(kind=8), intent(inout) :: x(ldx,*)
  integer :: i

!$omp parallel do
  do i = 1, nrows
    x(1:2,i) = alpha(:)*x(1:2,i)
  end do
end subroutine dscal_strided_2

subroutine dscal_strided_4(nrows, alpha, x, ldx)
  implicit none
  integer, intent(in) :: nrows, ldx
  real(kind=8), intent(in) :: alpha(4)
  real(kind=8), intent(inout) :: x(ldx,*)
  integer :: i

!$omp parallel do
  do i = 1, nrows
    x(1:4,i) = alpha(:)*x(1:4,i)
  end do
end subroutine dscal_strided_4

subroutine dscal_strided_8(nrows, alpha, x, ldx)
  implicit none
  integer, intent(in) :: nrows, ldx
  real(kind=8), intent(in) :: alpha(8)
  real(kind=8), intent(inout) :: x(ldx,*)
  integer :: i

!$omp parallel do
  do i = 1, nrows
    x(1:8,i) = alpha(:)*x(1:8,i)
  end do
end subroutine dscal_strided_8

subroutine dscal_general(nrows, nvec, alpha, x, ldx)
  implicit none
  integer, intent(in) :: nrows, nvec, ldx
  real(kind=8), intent(in) :: alpha(nvec)
  real(kind=8), intent(inout) :: x(ldx,*)
  integer :: i

!$omp parallel do
  do i = 1, nrows
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

  if( beta .eq. 0 ) then
!$omp parallel do
    do i = 1, nrows
      y(i) = alpha*x(i)
    end do
  else
!$omp parallel do
    do i = 1, nrows
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

  if( beta .eq. 0 ) then
!$omp parallel do
    do i = 1, nrows
      y(:,i) = alpha(:)*x(:,i)
    end do
  else
!$omp parallel do
    do i = 1, nrows
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

  if( beta .eq. 0 ) then
!$omp parallel do
    do i = 1, nrows
      y(:,i) = alpha(:)*x(:,i)
    end do
  else
!$omp parallel do
    do i = 1, nrows
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

  if( beta .eq. 0 ) then
!$omp parallel do
    do i = 1, nrows
      y(:,i) = alpha(:)*x(:,i) + beta*y(:,i)
    end do
  else
!$omp parallel do
    do i = 1, nrows
      y(:,i) = alpha(:)*x(:,i) + beta*y(:,i)
    end do
  end if
end subroutine daxpby_8


subroutine daxpby_strided_1(nrows, alpha, x, ldx, beta, y, ldy)
  implicit none
  integer, intent(in) :: nrows, ldx, ldy
  real(kind=8), intent(in) :: alpha
  real(kind=8), intent(in) :: beta
  real(kind=8), intent(in) :: x(ldx,nrows)
  real(kind=8), intent(inout) :: y(ldy,nrows)
  integer :: i

  if( beta .eq. 0 ) then
!$omp parallel do
    do i = 1, nrows
      y(1,i) = alpha*x(1,i)
    end do
  else
!$omp parallel do
    do i = 1, nrows
      y(1,i) = alpha*x(1,i) + beta*y(1,i)
    end do
  end if
end subroutine daxpby_strided_1

subroutine daxpby_strided_2(nrows, alpha, x, ldx, beta, y, ldy)
  implicit none
  integer, intent(in) :: nrows, ldx, ldy
  real(kind=8), intent(in) :: alpha(2)
  real(kind=8), intent(in) :: beta
  real(kind=8), intent(in) :: x(ldx,nrows)
  real(kind=8), intent(inout) :: y(ldy,nrows)
  integer :: i

  if( beta .eq. 0 ) then
!$omp parallel do
    do i = 1, nrows
      y(1:2,i) = alpha*x(1:2,i)
    end do
  else
!$omp parallel do
    do i = 1, nrows
      y(1:2,i) = alpha*x(1:2,i) + beta*y(1:2,i)
    end do
  end if
end subroutine daxpby_strided_2

subroutine daxpby_strided_4(nrows, alpha, x, ldx, beta, y, ldy)
  implicit none
  integer, intent(in) :: nrows, ldx, ldy
  real(kind=8), intent(in) :: alpha(4)
  real(kind=8), intent(in) :: beta
  real(kind=8), intent(in) :: x(ldx,nrows)
  real(kind=8), intent(inout) :: y(ldy,nrows)
  integer :: i

  if( beta .eq. 0 ) then
!$omp parallel do
    do i = 1, nrows
      y(1:4,i) = alpha*x(1:4,i)
    end do
  else
!$omp parallel do
    do i = 1, nrows
      y(1:4,i) = alpha*x(1:4,i) + beta*y(1:4,i)
    end do
  end if
end subroutine daxpby_strided_4

subroutine daxpby_strided_8(nrows, alpha, x, ldx, beta, y, ldy)
  implicit none
  integer, intent(in) :: nrows, ldx, ldy
  real(kind=8), intent(in) :: alpha(8)
  real(kind=8), intent(in) :: beta
  real(kind=8), intent(in) :: x(ldx,nrows)
  real(kind=8), intent(inout) :: y(ldy,nrows)
  integer :: i

  if( beta .eq. 0 ) then
!$omp parallel do
    do i = 1, nrows
      y(1:8,i) = alpha*x(1:8,i)
    end do
  else
!$omp parallel do
    do i = 1, nrows
      y(1:8,i) = alpha*x(1:8,i) + beta*y(1:8,i)
    end do
  end if
end subroutine daxpby_strided_8

subroutine daxpby_generic(nrows, nvec, alpha, x, ldx, beta, y, ldy)
  implicit none
  integer, intent(in) :: nrows, nvec, ldx, ldy
  real(kind=8), intent(in) :: alpha(nvec)
  real(kind=8), intent(in) :: beta
  real(kind=8), intent(in) :: x(ldx,nrows)
  real(kind=8), intent(inout) :: y(ldy,nrows)
  integer :: i

  if( beta .eq. 0 ) then
!$omp parallel do
    do i = 1, nrows
      y(1:nvec,i) = alpha*x(1:nvec,i)
    end do
  else
!$omp parallel do
    do i = 1, nrows
      y(1:nvec,i) = alpha*x(1:nvec,i) + beta*y(1:nvec,i)
    end do
  end if
end subroutine daxpby_generic

