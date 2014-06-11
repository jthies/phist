!! Kaczmarz sweep, node local, OpenMP to be added.
!! isweep: +1 -> forward,
!!         -1 -> backward Kaczmarz sweep

!! non-zero rhs b, B=I, generic kernel for any nvec
subroutine dkswp_I_generic(isweep, nvec, nlocal, nhalo, ncols, nnz, row_ptr, 
halo_ptr, col_idx, val, shifts, x, ldx, halo, b, ldb)
  implicit none

  integer, intent(in) :: isweep, nvec, nlocal, nhalo, ncols, ldx, ldy
  integer(kind=8), intent(in) :: nnz
  complex(kind=8), intent(in) :: shifts(nvec)
  integer(kind=8), intent(in) :: row_ptr(nlocal+1), halo_ptr(nlocal)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(in) :: x(ldx,*), halo(nvec,nhalo)
  real(kind=8), intent(inout) :: y(ldy,*)
  real(kind=8) :: tmp(nvec)
  integer :: i
  integer(kind=8) :: j
  
  integer :: istart, iend
  
  if (isweep=1) then
    istart=1
    iend=nlocal
  else if (isweep=-1) then
    istart=nlocal
    iend=1
  end if

!!!$omp parallel do private(tmp) schedule(static)
  do i = istart, iend, isweep
    tmp(:) = shifts*x(1:nvec,i)
    do j = row_ptr(i), halo_ptr(i)-1, 1
      tmp(:) = tmp(:) + val(j)*x(1:nvec,col_idx(j))
    end do
    do j = halo_ptr(i), row_ptr(i+1)-1, 1
      tmp(:) = tmp(:) + val(j)*halo(:,col_idx(j))
    end do
    y(1:nvec,i) = alpha*tmp(:) + beta*y(1:nvec,i)
  end do
end subroutine dzkswp_I_generic

subroutine dzkswp_I_bzero_generic(isweep,nvec, nlocal, nhalo, ncols, nnz, row_ptr, 
halo_ptr, col_idx, val, shifts, x, ldx, halo, beta)
  implicit none
end subroutine dkswp_I_bzero_generic



