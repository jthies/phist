! hopefully fast spMVM kernels for beta != 0 (and thus without nontemporary stores)

subroutine dspmvm_NT_1(nlocal, nhalo, ncols, nnz, alpha, row_ptr, halo_ptr, col_idx, val, x, halo, y)
  implicit none
  integer, intent(in) :: nlocal, nhalo, ncols
  integer(kind=8), intent(in) :: nnz
  real(kind=8), intent(in) :: alpha
  integer(kind=8), intent(in) :: row_ptr(nlocal+1), halo_ptr(nlocal)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(in) :: x(1,ncols), halo(1,nhalo)
  real(kind=8), intent(inout) :: y(1,nlocal)

  interface
    subroutine dspmvm_nt_1_c(nlocal, alpha, row_ptr, halo_ptr, col_idx, val, x, halo, y) bind(C)
      use, intrinsic :: iso_c_binding
      integer(C_INT), value :: nlocal
      real(C_DOUBLE), value :: alpha
      integer(C_LONG), intent(in) :: row_ptr(*), halo_ptr(*)
      integer(C_INT), intent(in) :: col_idx(*)
      real(C_DOUBLE), intent(in) :: val(*), x(*), halo(*)
      real(C_DOUBLE), intent(inout) :: y(*)
    end subroutine dspmvm_nt_1_c
  end interface

  call dspmvm_nt_1_c(nlocal,alpha,row_ptr,halo_ptr,col_idx,val,x,halo,y)

end subroutine dspmvm_NT_1

subroutine dspmvm_NT_2(nlocal, nhalo, ncols, nnz, alpha, row_ptr, halo_ptr, col_idx, val, x, halo, y, ldy)
  implicit none
  integer, intent(in) :: nlocal, nhalo, ncols, ldy
  integer(kind=8), intent(in) :: nnz
  real(kind=8), intent(in) :: alpha
  integer(kind=8), intent(in) :: row_ptr(nlocal+1), halo_ptr(nlocal)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(in) :: x(2,ncols), halo(2,nhalo)
  real(kind=8), intent(inout) :: y(ldy,nlocal)

  interface
    subroutine dspmvm_nt_2_c(nlocal, alpha, row_ptr, halo_ptr, col_idx, val, x, halo, y, ldy) bind(C)
      use, intrinsic :: iso_c_binding
      integer(C_INT), value :: nlocal, ldy
      real(C_DOUBLE), value :: alpha
      integer(C_LONG), intent(in) :: row_ptr(*), halo_ptr(*)
      integer(C_INT), intent(in) :: col_idx(*)
      real(C_DOUBLE), intent(in) :: val(*), x(*), halo(*)
      real(C_DOUBLE), intent(inout) :: y(*)
    end subroutine dspmvm_nt_2_c
  end interface

  call dspmvm_nt_2_c(nlocal,alpha,row_ptr,halo_ptr,col_idx,val,x,halo,y,ldy)

end subroutine dspmvm_NT_2


subroutine dspmvm_NT_4(nlocal, nhalo, ncols, nnz, alpha, row_ptr, halo_ptr, col_idx, val, x, halo, y, ldy)
  implicit none
  integer, intent(in) :: nlocal, nhalo, ncols, ldy
  integer(kind=8), intent(in) :: nnz
  real(kind=8), intent(in) :: alpha
  integer(kind=8), intent(in) :: row_ptr(nlocal+1), halo_ptr(nlocal)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(in) :: x(4,ncols), halo(4,nhalo)
  real(kind=8), intent(inout) :: y(ldy,nlocal)

  interface
    subroutine dspmvm_nt_4_c(nlocal, alpha, row_ptr, halo_ptr, col_idx, val, x, halo, y, ldy) bind(C)
      use, intrinsic :: iso_c_binding
      integer(C_INT), value :: nlocal, ldy
      real(C_DOUBLE), value :: alpha
      integer(C_LONG), intent(in) :: row_ptr(*), halo_ptr(*)
      integer(C_INT), intent(in) :: col_idx(*)
      real(C_DOUBLE), intent(in) :: val(*), x(*), halo(*)
      real(C_DOUBLE), intent(inout) :: y(*)
    end subroutine dspmvm_nt_4_c
  end interface

  call dspmvm_nt_4_c(nlocal,alpha,row_ptr,halo_ptr,col_idx,val,x,halo,y,ldy)

end subroutine dspmvm_NT_4


subroutine dspmvm_NT_8(nlocal, nhalo, ncols, nnz, alpha, row_ptr, halo_ptr, col_idx, val, x, halo, y, ldy)
  implicit none
  integer, intent(in) :: nlocal, nhalo, ncols, ldy
  integer(kind=8), intent(in) :: nnz
  real(kind=8), intent(in) :: alpha
  integer(kind=8), intent(in) :: row_ptr(nlocal+1), halo_ptr(nlocal)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(in) :: x(8,ncols), halo(8,nhalo)
  real(kind=8), intent(inout) :: y(ldy,nlocal)

  interface
    subroutine dspmvm_nt_8_c(nlocal, alpha, row_ptr, halo_ptr, col_idx, val, x, halo, y, ldy) bind(C)
      use, intrinsic :: iso_c_binding
      integer(C_INT), value :: nlocal, ldy
      real(C_DOUBLE), value :: alpha
      integer(C_LONG), intent(in) :: row_ptr(*), halo_ptr(*)
      integer(C_INT), intent(in) :: col_idx(*)
      real(C_DOUBLE), intent(in) :: val(*), x(*), halo(*)
      real(C_DOUBLE), intent(inout) :: y(*)
    end subroutine dspmvm_nt_8_c
  end interface

  call dspmvm_nt_8_c(nlocal,alpha,row_ptr,halo_ptr,col_idx,val,x,halo,y,ldy)

end subroutine dspmvm_NT_8


subroutine dspmvm_NT_strided_2(nlocal, nhalo, ncols, nnz, alpha, row_ptr, halo_ptr, col_idx, val, x, ldx, halo, y, ldy)
  implicit none
  integer, intent(in) :: nlocal, nhalo, ncols, ldy, ldx
  integer(kind=8), intent(in) :: nnz
  real(kind=8), intent(in) :: alpha
  integer(kind=8), intent(in) :: row_ptr(nlocal+1), halo_ptr(nlocal)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(in) :: x(ldx,ncols), halo(2,nhalo)
  real(kind=8), intent(inout) :: y(ldy,nlocal)

  interface
    subroutine dspmvm_nt_strided_2_c(nlocal, alpha, row_ptr, halo_ptr, col_idx, val, x, ldx, halo, y, ldy) bind(C)
      use, intrinsic :: iso_c_binding
      integer(C_INT), value :: nlocal, ldy, ldx
      real(C_DOUBLE), value :: alpha
      integer(C_LONG), intent(in) :: row_ptr(*), halo_ptr(*)
      integer(C_INT), intent(in) :: col_idx(*)
      real(C_DOUBLE), intent(in) :: val(*), x(*), halo(*)
      real(C_DOUBLE), intent(inout) :: y(*)
    end subroutine dspmvm_nt_strided_2_c
  end interface

  call dspmvm_nt_strided_2_c(nlocal,alpha,row_ptr,halo_ptr,col_idx,val,x,ldx,halo,y,ldy)

end subroutine dspmvm_NT_strided_2


subroutine dspmvm_NT_strided_4(nlocal, nhalo, ncols, nnz, alpha, row_ptr, halo_ptr, col_idx, val, x, ldx, halo, y, ldy)
  implicit none
  integer, intent(in) :: nlocal, nhalo, ncols, ldy, ldx
  integer(kind=8), intent(in) :: nnz
  real(kind=8), intent(in) :: alpha
  integer(kind=8), intent(in) :: row_ptr(nlocal+1), halo_ptr(nlocal)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(in) :: x(ldx,ncols), halo(4,nhalo)
  real(kind=8), intent(inout) :: y(ldy,nlocal)

  interface
    subroutine dspmvm_nt_strided_4_c(nlocal, alpha, row_ptr, halo_ptr, col_idx, val, x, ldx, halo, y, ldy) bind(C)
      use, intrinsic :: iso_c_binding
      integer(C_INT), value :: nlocal, ldy, ldx
      real(C_DOUBLE), value :: alpha
      integer(C_LONG), intent(in) :: row_ptr(*), halo_ptr(*)
      integer(C_INT), intent(in) :: col_idx(*)
      real(C_DOUBLE), intent(in) :: val(*), x(*), halo(*)
      real(C_DOUBLE), intent(inout) :: y(*)
    end subroutine dspmvm_nt_strided_4_c
  end interface

  call dspmvm_nt_strided_4_c(nlocal,alpha,row_ptr,halo_ptr,col_idx,val,x,ldx,halo,y,ldy)

end subroutine dspmvm_NT_strided_4


subroutine dspmvm_NT_strided_8(nlocal, nhalo, ncols, nnz, alpha, row_ptr, halo_ptr, col_idx, val, x, ldx, halo, y, ldy)
  implicit none
  integer, intent(in) :: nlocal, nhalo, ncols, ldy, ldx
  integer(kind=8), intent(in) :: nnz
  real(kind=8), intent(in) :: alpha
  integer(kind=8), intent(in) :: row_ptr(nlocal+1), halo_ptr(nlocal)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(in) :: x(ldx,ncols), halo(8,nhalo)
  real(kind=8), intent(inout) :: y(ldy,nlocal)

  interface
    subroutine dspmvm_nt_strided_8_c(nlocal, alpha, row_ptr, halo_ptr, col_idx, val, x, ldx, halo, y, ldy) bind(C)
      use, intrinsic :: iso_c_binding
      integer(C_INT), value :: nlocal, ldy, ldx
      real(C_DOUBLE), value :: alpha
      integer(C_LONG), intent(in) :: row_ptr(*), halo_ptr(*)
      integer(C_INT), intent(in) :: col_idx(*)
      real(C_DOUBLE), intent(in) :: val(*), x(*), halo(*)
      real(C_DOUBLE), intent(inout) :: y(*)
    end subroutine dspmvm_nt_strided_8_c
  end interface

  call dspmvm_nt_strided_8_c(nlocal,alpha,row_ptr,halo_ptr,col_idx,val,x,ldx,halo,y,ldy)

end subroutine dspmvm_NT_strided_8



subroutine dspmvm_1(nlocal, nhalo, ncols, nnz, alpha, row_ptr, halo_ptr, col_idx, val, x, halo, beta, y)
  implicit none
  integer, intent(in) :: nlocal, nhalo, ncols
  integer(kind=8), intent(in) :: nnz
  real(kind=8), intent(in) :: alpha, beta
  integer(kind=8), intent(in) :: row_ptr(nlocal+1), halo_ptr(nlocal)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(in) :: x(1,ncols), halo(1,nhalo)
  real(kind=8), intent(inout) :: y(1,nlocal)
  real(kind=8) :: tmp(1)
  integer :: i
  integer(kind=8) :: j

!$omp parallel do private(tmp)
  do i = 1, nlocal
    tmp(:) = 0.
    do j = row_ptr(i), halo_ptr(i)-1, 1
      tmp(:) = tmp(:) + val(j)*x(:,col_idx(j))
    end do
    do j = halo_ptr(i), row_ptr(i+1)-1, 1
      tmp(:) = tmp(:) + val(j)*halo(:,col_idx(j))
    end do
    y(:,i) = alpha*tmp(:) + beta*y(:,i)
  end do
end subroutine dspmvm_1

subroutine dspmvm_2(nlocal, nhalo, ncols, nnz, alpha, row_ptr, halo_ptr, col_idx, val, x, halo, beta, y)
  implicit none
  integer, intent(in) :: nlocal, nhalo, ncols
  integer(kind=8), intent(in) :: nnz
  real(kind=8), intent(in) :: alpha, beta
  integer(kind=8), intent(in) :: row_ptr(nlocal+1), halo_ptr(nlocal)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(in) :: x(2,ncols), halo(2,nhalo)
  real(kind=8), intent(inout) :: y(2,nlocal)
  real(kind=8) :: tmp(2)
  integer :: i
  integer(kind=8) :: j

!$omp parallel do private(tmp)
  do i = 1, nlocal
    tmp(:) = 0.
    do j = row_ptr(i), halo_ptr(i)-1, 1
      tmp(:) = tmp(:) + val(j)*x(:,col_idx(j))
    end do
    do j = halo_ptr(i), row_ptr(i+1)-1, 1
      tmp(:) = tmp(:) + val(j)*halo(:,col_idx(j))
    end do
    y(:,i) = alpha*tmp(:) + beta*y(:,i)
  end do
end subroutine dspmvm_2

subroutine dspmvm_4(nlocal, nhalo, ncols, nnz, alpha, row_ptr, halo_ptr, col_idx, val, x, halo, beta, y)
  implicit none
  integer, intent(in) :: nlocal, nhalo, ncols
  integer(kind=8), intent(in) :: nnz
  real(kind=8), intent(in) :: alpha, beta
  integer(kind=8), intent(in) :: row_ptr(nlocal+1), halo_ptr(nlocal)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(in) :: x(4,ncols), halo(4,nhalo)
  real(kind=8), intent(inout) :: y(4,nlocal)
  real(kind=8) :: tmp(4)
  integer :: i
  integer(kind=8) :: j

!$omp parallel do private(tmp)
  do i = 1, nlocal
    tmp(:) = 0.
    do j = row_ptr(i), halo_ptr(i)-1, 1
      tmp(:) = tmp(:) + val(j)*x(:,col_idx(j))
    end do
    do j = halo_ptr(i), row_ptr(i+1)-1, 1
      tmp(:) = tmp(:) + val(j)*halo(:,col_idx(j))
    end do
    y(:,i) = alpha*tmp(:) + beta*y(:,i)
  end do
end subroutine dspmvm_4

subroutine dspmvm_8(nlocal, nhalo, ncols, nnz, alpha, row_ptr, halo_ptr, col_idx, val, x, halo, beta, y)
  implicit none
  integer, intent(in) :: nlocal, nhalo, ncols
  integer(kind=8), intent(in) :: nnz
  real(kind=8), intent(in) :: alpha, beta
  integer(kind=8), intent(in) :: row_ptr(nlocal+1), halo_ptr(nlocal)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(in) :: x(8,ncols), halo(8,nhalo)
  real(kind=8), intent(inout) :: y(8,nlocal)
  real(kind=8) :: tmp(8)
  integer :: i
  integer(kind=8) :: j

!$omp parallel do private(tmp)
  do i = 1, nlocal
    tmp(:) = 0.
    do j = row_ptr(i), halo_ptr(i)-1, 1
      tmp(:) = tmp(:) + val(j)*x(:,col_idx(j))
    end do
    do j = halo_ptr(i), row_ptr(i+1)-1, 1
      tmp(:) = tmp(:) + val(j)*halo(:,col_idx(j))
    end do
    y(:,i) = alpha*tmp(:) + beta*y(:,i)
  end do
end subroutine dspmvm_8



subroutine dspmvm_strided_1(nlocal, nhalo, ncols, nnz, alpha, row_ptr, halo_ptr, col_idx, val, x, ldx, halo, beta, y, ldy)
  implicit none
  integer, intent(in) :: nlocal, nhalo, ncols, ldx, ldy
  integer(kind=8), intent(in) :: nnz
  real(kind=8), intent(in) :: alpha, beta
  integer(kind=8), intent(in) :: row_ptr(nlocal+1), halo_ptr(nlocal)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(in) :: x(ldx,*), halo(1,nhalo)
  real(kind=8), intent(inout) :: y(ldy,*)
  real(kind=8) :: tmp(1)
  integer :: i
  integer(kind=8) :: j

!$omp parallel do private(tmp)
  do i = 1, nlocal
    tmp(:) = 0.
    do j = row_ptr(i), halo_ptr(i)-1, 1
      tmp(:) = tmp(:) + val(j)*x(1:1,col_idx(j))
    end do
    do j = halo_ptr(i), row_ptr(i+1)-1, 1
      tmp(:) = tmp(:) + val(j)*halo(:,col_idx(j))
    end do
    y(1:1,i) = alpha*tmp(:) + beta*y(1:1,i)
  end do
end subroutine dspmvm_strided_1

subroutine dspmvm_strided_2(nlocal, nhalo, ncols, nnz, alpha, row_ptr, halo_ptr, col_idx, val, x, ldx, halo, beta, y, ldy)
  implicit none
  integer, intent(in) :: nlocal, nhalo, ncols, ldx, ldy
  integer(kind=8), intent(in) :: nnz
  real(kind=8), intent(in) :: alpha, beta
  integer(kind=8), intent(in) :: row_ptr(nlocal+1), halo_ptr(nlocal)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(in) :: x(ldx,*), halo(2,nhalo)
  real(kind=8), intent(inout) :: y(ldy,*)
  real(kind=8) :: tmp(2)
  integer :: i
  integer(kind=8) :: j

!$omp parallel do private(tmp)
  do i = 1, nlocal
    tmp(:) = 0.
    do j = row_ptr(i), halo_ptr(i)-1, 1
      tmp(:) = tmp(:) + val(j)*x(1:2,col_idx(j))
    end do
    do j = halo_ptr(i), row_ptr(i+1)-1, 1
      tmp(:) = tmp(:) + val(j)*halo(:,col_idx(j))
    end do
    y(1:2,i) = alpha*tmp(:) + beta*y(1:2,i)
  end do
end subroutine dspmvm_strided_2

subroutine dspmvm_strided_4(nlocal, nhalo, ncols, nnz, alpha, row_ptr, halo_ptr, col_idx, val, x, ldx, halo, beta, y, ldy)
  implicit none
  integer, intent(in) :: nlocal, nhalo, ncols, ldx, ldy
  integer(kind=8), intent(in) :: nnz
  real(kind=8), intent(in) :: alpha, beta
  integer(kind=8), intent(in) :: row_ptr(nlocal+1), halo_ptr(nlocal)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(in) :: x(ldx,*), halo(4,nhalo)
  real(kind=8), intent(inout) :: y(ldy,*)
  real(kind=8) :: tmp(4)
  integer :: i
  integer(kind=8) :: j

!$omp parallel do private(tmp)
  do i = 1, nlocal
    tmp(:) = 0.
    do j = row_ptr(i), halo_ptr(i)-1, 1
      tmp(:) = tmp(:) + val(j)*x(1:4,col_idx(j))
    end do
    do j = halo_ptr(i), row_ptr(i+1)-1, 1
      tmp(:) = tmp(:) + val(j)*halo(:,col_idx(j))
    end do
    y(1:4,i) = alpha*tmp(:) + beta*y(1:4,i)
  end do
end subroutine dspmvm_strided_4

subroutine dspmvm_strided_8(nlocal, nhalo, ncols, nnz, alpha, row_ptr, halo_ptr, col_idx, val, x, ldx, halo, beta, y, ldy)
  implicit none
  integer, intent(in) :: nlocal, nhalo, ncols, ldx, ldy
  integer(kind=8), intent(in) :: nnz
  real(kind=8), intent(in) :: alpha, beta
  integer(kind=8), intent(in) :: row_ptr(nlocal+1), halo_ptr(nlocal)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(in) :: x(ldx,*), halo(8,nhalo)
  real(kind=8), intent(inout) :: y(ldy,*)
  real(kind=8) :: tmp(8)
  integer :: i
  integer(kind=8) :: j

!$omp parallel do private(tmp)
  do i = 1, nlocal
    tmp(:) = 0.
    do j = row_ptr(i), halo_ptr(i)-1, 1
      tmp(:) = tmp(:) + val(j)*x(1:8,col_idx(j))
    end do
    do j = halo_ptr(i), row_ptr(i+1)-1, 1
      tmp(:) = tmp(:) + val(j)*halo(:,col_idx(j))
    end do
    y(1:8,i) = alpha*tmp(:) + beta*y(1:8,i)
  end do
end subroutine dspmvm_strided_8


subroutine dspmvm_generic(nvec, nlocal, nhalo, ncols, nnz, alpha, row_ptr, halo_ptr, col_idx, val, x, ldx, halo, beta, y, ldy)
  implicit none
  integer, intent(in) :: nvec, nlocal, nhalo, ncols, ldx, ldy
  integer(kind=8), intent(in) :: nnz
  real(kind=8), intent(in) :: alpha, beta
  integer(kind=8), intent(in) :: row_ptr(nlocal+1), halo_ptr(nlocal)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(in) :: x(ldx,*), halo(nvec,nhalo)
  real(kind=8), intent(inout) :: y(ldy,*)
  real(kind=8) :: tmp(nvec)
  integer :: i
  integer(kind=8) :: j

!$omp parallel do private(tmp)
  do i = 1, nlocal
    tmp(:) = 0.
    do j = row_ptr(i), halo_ptr(i)-1, 1
      tmp(:) = tmp(:) + val(j)*x(1:nvec,col_idx(j))
    end do
    do j = halo_ptr(i), row_ptr(i+1)-1, 1
      tmp(:) = tmp(:) + val(j)*halo(:,col_idx(j))
    end do
    y(1:nvec,i) = alpha*tmp(:) + beta*y(1:nvec,i)
  end do
end subroutine dspmvm_generic


