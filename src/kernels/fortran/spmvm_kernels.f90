! hopefully fast spMVM kernels for beta != 0 (and thus without nontemporary stores)

subroutine dspmvm_NT_1(nrows, ncols, nnz, alpha, row_ptr, col_idx, val, x, y)
  integer, intent(in) :: nrows, ncols, nnz
  real(kind=8), intent(in) :: alpha
  integer, intent(in) :: row_ptr(nrows+1)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(in) :: x(1,ncols)
  real(kind=8), intent(inout) :: y(1,nrows)

  interface
    subroutine dspmvm_nt_1_c(nrows, alpha, row_ptr, col_idx, val, x, y) bind(C)
      use, intrinsic :: iso_c_binding
      integer(C_INT), value :: nrows
      real(C_DOUBLE), value :: alpha
      integer(C_INT), intent(in) :: row_ptr(*), col_idx(*)
      real(C_DOUBLE), intent(in) :: val(*), x(*)
      real(C_DOUBLE), intent(inout) :: y(*)
    end subroutine dspmvm_nt_1_c
  end interface

  call dspmvm_nt_1_c(nrows,alpha,row_ptr,col_idx,val,x,y)

end subroutine dspmvm_NT_1

subroutine dspmvm_NT_2(nrows, ncols, nnz, alpha, row_ptr, col_idx, val, x, y, ldy)
  integer, intent(in) :: nrows, ncols, nnz, ldy
  real(kind=8), intent(in) :: alpha
  integer, intent(in) :: row_ptr(nrows+1)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(in) :: x(2,ncols)
  real(kind=8), intent(inout) :: y(ldy,nrows)

  interface
    subroutine dspmvm_nt_2_c(nrows, alpha, row_ptr, col_idx, val, x, y, ldy) bind(C)
      use, intrinsic :: iso_c_binding
      integer(C_INT), value :: nrows, ldy
      real(C_DOUBLE), value :: alpha
      integer(C_INT), intent(in) :: row_ptr(*), col_idx(*)
      real(C_DOUBLE), intent(in) :: val(*), x(*)
      real(C_DOUBLE), intent(inout) :: y(*)
    end subroutine dspmvm_nt_2_c
  end interface

  call dspmvm_nt_2_c(nrows,alpha,row_ptr,col_idx,val,x,y,ldy)

end subroutine dspmvm_NT_2


subroutine dspmvm_NT_4(nrows, ncols, nnz, alpha, row_ptr, col_idx, val, x, y, ldy)
  integer, intent(in) :: nrows, ncols, nnz, ldy
  real(kind=8), intent(in) :: alpha
  integer, intent(in) :: row_ptr(nrows+1)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(in) :: x(4,ncols)
  real(kind=8), intent(inout) :: y(ldy,nrows)

  interface
    subroutine dspmvm_nt_4_c(nrows, alpha, row_ptr, col_idx, val, x, y, ldy) bind(C)
      use, intrinsic :: iso_c_binding
      integer(C_INT), value :: nrows, ldy
      real(C_DOUBLE), value :: alpha
      integer(C_INT), intent(in) :: row_ptr(*), col_idx(*)
      real(C_DOUBLE), intent(in) :: val(*), x(*)
      real(C_DOUBLE), intent(inout) :: y(*)
    end subroutine dspmvm_nt_4_c
  end interface

  call dspmvm_nt_4_c(nrows,alpha,row_ptr,col_idx,val,x,y,ldy)

end subroutine dspmvm_NT_4


subroutine dspmvm_NT_8(nrows, ncols, nnz, alpha, row_ptr, col_idx, val, x, y, ldy)
  integer, intent(in) :: nrows, ncols, nnz, ldy
  real(kind=8), intent(in) :: alpha
  integer, intent(in) :: row_ptr(nrows+1)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(in) :: x(8,ncols)
  real(kind=8), intent(inout) :: y(ldy,nrows)

  interface
    subroutine dspmvm_nt_8_c(nrows, alpha, row_ptr, col_idx, val, x, y, ldy) bind(C)
      use, intrinsic :: iso_c_binding
      integer(C_INT), value :: nrows, ldy
      real(C_DOUBLE), value :: alpha
      integer(C_INT), intent(in) :: row_ptr(*), col_idx(*)
      real(C_DOUBLE), intent(in) :: val(*), x(*)
      real(C_DOUBLE), intent(inout) :: y(*)
    end subroutine dspmvm_nt_8_c
  end interface

  call dspmvm_nt_8_c(nrows,alpha,row_ptr,col_idx,val,x,y,ldy)

end subroutine dspmvm_NT_8


subroutine dspmvm_NT_strided_2(nrows, ncols, nnz, alpha, row_ptr, col_idx, val, x, ldx, y, ldy)
  integer, intent(in) :: nrows, ncols, nnz, ldy, ldx
  real(kind=8), intent(in) :: alpha
  integer, intent(in) :: row_ptr(nrows+1)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(in) :: x(ldx,ncols)
  real(kind=8), intent(inout) :: y(ldy,nrows)

  interface
    subroutine dspmvm_nt_strided_2_c(nrows, alpha, row_ptr, col_idx, val, x, ldx, y, ldy) bind(C)
      use, intrinsic :: iso_c_binding
      integer(C_INT), value :: nrows, ldy, ldx
      real(C_DOUBLE), value :: alpha
      integer(C_INT), intent(in) :: row_ptr(*), col_idx(*)
      real(C_DOUBLE), intent(in) :: val(*), x(*)
      real(C_DOUBLE), intent(inout) :: y(*)
    end subroutine dspmvm_nt_strided_2_c
  end interface

  call dspmvm_nt_strided_2_c(nrows,alpha,row_ptr,col_idx,val,x,ldx,y,ldy)

end subroutine dspmvm_NT_strided_2


subroutine dspmvm_NT_strided_4(nrows, ncols, nnz, alpha, row_ptr, col_idx, val, x, ldx, y, ldy)
  integer, intent(in) :: nrows, ncols, nnz, ldy, ldx
  real(kind=8), intent(in) :: alpha
  integer, intent(in) :: row_ptr(nrows+1)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(in) :: x(ldx,ncols)
  real(kind=8), intent(inout) :: y(ldy,nrows)

  interface
    subroutine dspmvm_nt_strided_4_c(nrows, alpha, row_ptr, col_idx, val, x, ldx, y, ldy) bind(C)
      use, intrinsic :: iso_c_binding
      integer(C_INT), value :: nrows, ldy, ldx
      real(C_DOUBLE), value :: alpha
      integer(C_INT), intent(in) :: row_ptr(*), col_idx(*)
      real(C_DOUBLE), intent(in) :: val(*), x(*)
      real(C_DOUBLE), intent(inout) :: y(*)
    end subroutine dspmvm_nt_strided_4_c
  end interface

  call dspmvm_nt_strided_4_c(nrows,alpha,row_ptr,col_idx,val,x,ldx,y,ldy)

end subroutine dspmvm_NT_strided_4


subroutine dspmvm_NT_strided_8(nrows, ncols, nnz, alpha, row_ptr, col_idx, val, x, ldx, y, ldy)
  integer, intent(in) :: nrows, ncols, nnz, ldy, ldx
  real(kind=8), intent(in) :: alpha
  integer, intent(in) :: row_ptr(nrows+1)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(in) :: x(ldx,ncols)
  real(kind=8), intent(inout) :: y(ldy,nrows)

  interface
    subroutine dspmvm_nt_strided_8_c(nrows, alpha, row_ptr, col_idx, val, x, ldx, y, ldy) bind(C)
      use, intrinsic :: iso_c_binding
      integer(C_INT), value :: nrows, ldy, ldx
      real(C_DOUBLE), value :: alpha
      integer(C_INT), intent(in) :: row_ptr(*), col_idx(*)
      real(C_DOUBLE), intent(in) :: val(*), x(*)
      real(C_DOUBLE), intent(inout) :: y(*)
    end subroutine dspmvm_nt_strided_8_c
  end interface

  call dspmvm_nt_strided_8_c(nrows,alpha,row_ptr,col_idx,val,x,ldx,y,ldy)

end subroutine dspmvm_NT_strided_8



subroutine dspmvm_1(nrows, ncols, nnz, alpha, row_ptr, col_idx, val, x, beta, y)
  integer, intent(in) :: nrows, ncols, nnz
  real(kind=8), intent(in) :: alpha, beta
  integer, intent(in) :: row_ptr(nrows+1)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(in) :: x(1,ncols)
  real(kind=8), intent(inout) :: y(1,nrows)
  real(kind=8) :: tmp(1)
  integer :: i, j

!$omp parallel do private(tmp)
  do i = 1, nrows
    tmp(:) = 0.
    do j = row_ptr(i), row_ptr(i+1)-1, 1
      tmp(:) = tmp(:) + val(j)*x(:,col_idx(j))
    end do
    y(:,i) = alpha*tmp(:) + beta*y(:,i)
  end do
end subroutine dspmvm_1

subroutine dspmvm_2(nrows, ncols, nnz, alpha, row_ptr, col_idx, val, x, beta, y)
  integer, intent(in) :: nrows, ncols, nnz
  real(kind=8), intent(in) :: alpha, beta
  integer, intent(in) :: row_ptr(nrows+1)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(in) :: x(2,ncols)
  real(kind=8), intent(inout) :: y(2,nrows)
  real(kind=8) :: tmp(2)
  integer :: i, j

!$omp parallel do private(tmp)
  do i = 1, nrows
    tmp(:) = 0.
    do j = row_ptr(i), row_ptr(i+1)-1, 1
      tmp(:) = tmp(:) + val(j)*x(:,col_idx(j))
    end do
    y(:,i) = alpha*tmp(:) + beta*y(:,i)
  end do
end subroutine dspmvm_2

subroutine dspmvm_4(nrows, ncols, nnz, alpha, row_ptr, col_idx, val, x, beta, y)
  integer, intent(in) :: nrows, ncols, nnz
  real(kind=8), intent(in) :: alpha, beta
  integer, intent(in) :: row_ptr(nrows+1)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(in) :: x(4,ncols)
  real(kind=8), intent(inout) :: y(4,nrows)
  real(kind=8) :: tmp(4)
  integer :: i, j

!$omp parallel do private(tmp)
  do i = 1, nrows
    tmp(:) = 0.
    do j = row_ptr(i), row_ptr(i+1)-1, 1
      tmp(:) = tmp(:) + val(j)*x(:,col_idx(j))
    end do
    y(:,i) = alpha*tmp(:) + beta*y(:,i)
  end do
end subroutine dspmvm_4

subroutine dspmvm_8(nrows, ncols, nnz, alpha, row_ptr, col_idx, val, x, beta, y)
  integer, intent(in) :: nrows, ncols, nnz
  real(kind=8), intent(in) :: alpha, beta
  integer, intent(in) :: row_ptr(nrows+1)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(in) :: x(8,ncols)
  real(kind=8), intent(inout) :: y(8,nrows)
  real(kind=8) :: tmp(8)
  integer :: i, j

!$omp parallel do private(tmp)
  do i = 1, nrows
    tmp(:) = 0.
    do j = row_ptr(i), row_ptr(i+1)-1, 1
      tmp(:) = tmp(:) + val(j)*x(:,col_idx(j))
    end do
    y(:,i) = alpha*tmp(:) + beta*y(:,i)
  end do
end subroutine dspmvm_8



subroutine dspmvm_strided_1(nrows, ncols, nnz, alpha, row_ptr, col_idx, val, x, ldx, beta, y, ldy)
  integer, intent(in) :: nrows, ncols, nnz, ldx, ldy
  real(kind=8), intent(in) :: alpha, beta
  integer, intent(in) :: row_ptr(nrows+1)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(in) :: x(ldx,*)
  real(kind=8), intent(inout) :: y(ldy,*)
  real(kind=8) :: tmp(1)
  integer :: i, j

!$omp parallel do private(tmp)
  do i = 1, nrows
    tmp(:) = 0.
    do j = row_ptr(i), row_ptr(i+1)-1, 1
      tmp(:) = tmp(:) + val(j)*x(1:1,col_idx(j))
    end do
    y(1:1,i) = alpha*tmp(:) + beta*y(1:1,i)
  end do
end subroutine dspmvm_strided_1

subroutine dspmvm_strided_2(nrows, ncols, nnz, alpha, row_ptr, col_idx, val, x, ldx, beta, y, ldy)
  integer, intent(in) :: nrows, ncols, nnz, ldx, ldy
  real(kind=8), intent(in) :: alpha, beta
  integer, intent(in) :: row_ptr(nrows+1)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(in) :: x(ldx,*)
  real(kind=8), intent(inout) :: y(ldy,*)
  real(kind=8) :: tmp(2)
  integer :: i, j

!$omp parallel do private(tmp)
  do i = 1, nrows
    tmp(:) = 0.
    do j = row_ptr(i), row_ptr(i+1)-1, 1
      tmp(:) = tmp(:) + val(j)*x(1:2,col_idx(j))
    end do
    y(1:2,i) = alpha*tmp(:) + beta*y(1:2,i)
  end do
end subroutine dspmvm_strided_2

subroutine dspmvm_strided_4(nrows, ncols, nnz, alpha, row_ptr, col_idx, val, x, ldx, beta, y, ldy)
  integer, intent(in) :: nrows, ncols, nnz, ldx, ldy
  real(kind=8), intent(in) :: alpha, beta
  integer, intent(in) :: row_ptr(nrows+1)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(in) :: x(ldx,*)
  real(kind=8), intent(inout) :: y(ldy,*)
  real(kind=8) :: tmp(4)
  integer :: i, j

!$omp parallel do private(tmp)
  do i = 1, nrows
    tmp(:) = 0.
    do j = row_ptr(i), row_ptr(i+1)-1, 1
      tmp(:) = tmp(:) + val(j)*x(1:4,col_idx(j))
    end do
    y(1:4,i) = alpha*tmp(:) + beta*y(1:4,i)
  end do
end subroutine dspmvm_strided_4

subroutine dspmvm_strided_8(nrows, ncols, nnz, alpha, row_ptr, col_idx, val, x, ldx, beta, y, ldy)
  integer, intent(in) :: nrows, ncols, nnz, ldx, ldy
  real(kind=8), intent(in) :: alpha, beta
  integer, intent(in) :: row_ptr(nrows+1)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(in) :: x(ldx,*)
  real(kind=8), intent(inout) :: y(ldy,*)
  real(kind=8) :: tmp(8)
  integer :: i, j

!$omp parallel do private(tmp)
  do i = 1, nrows
    tmp(:) = 0.
    do j = row_ptr(i), row_ptr(i+1)-1, 1
      tmp(:) = tmp(:) + val(j)*x(1:8,col_idx(j))
    end do
    y(1:8,i) = alpha*tmp(:) + beta*y(1:8,i)
  end do
end subroutine dspmvm_strided_8


subroutine dspmvm_generic(nrows, nvec, ncols, nnz, alpha, row_ptr, col_idx, val, x, ldx, beta, y, ldy)
  integer, intent(in) :: nrows, nvec, ncols, nnz, ldx, ldy
  real(kind=8), intent(in) :: alpha, beta
  integer, intent(in) :: row_ptr(nrows+1)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(in) :: x(ldx,*)
  real(kind=8), intent(inout) :: y(ldy,*)
  real(kind=8) :: tmp(nvec)
  integer :: i, j

!$omp parallel do private(tmp)
  do i = 1, nrows
    tmp(:) = 0.
    do j = row_ptr(i), row_ptr(i+1)-1, 1
      tmp(:) = tmp(:) + val(j)*x(1:nvec,col_idx(j))
    end do
    y(1:nvec,i) = alpha*tmp(:) + beta*y(1:nvec,i)
  end do
end subroutine dspmvm_generic


