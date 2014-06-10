!! compute for a given shift sigma=sigma_r+sigma_i the 2-norm of each row
!! of the matrix shift[j]*I-A and store the inverse of the results in the
!! columns j of the block vector nrmsai2i.
subroutine dkacz_norms_ai2i(nvec, nlocal, nhalo, ncols, nnz, row_ptr, halo_ptr, 
col_idx, val, shifts_r,shifts_i, 
nrmsai2i,istart)
  implicit none
  integer, intent(in) :: nvec, nlocal, nhalo, ncols
  integer(kind=8), intent(in) :: nnz
  c(kind=8), intent(in) :: shifts_r(nvec), shifts_i(nvec)
  integer(kind=8), intent(in) :: row_ptr(nlocal+1), halo_ptr(nlocal)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(out) :: nrmsai2(nlocal,nvec)
  integer :: i
  integer(kind=8) :: j

!TODO $omp parallel do private(tmp) schedule(static)
  do i = 1,nlocal
    nrms_ai2i(i,:)=0.d0
    !TODO - continue here!
    do j = row_ptr(i), halo_ptr(i)-1, 1
      tmp_r(:) = tmp_r(:) - val(j)*x_r(1:nvec,col_idx(j))
      tmp_i(:) = tmp_i(:) - val(j)*x_i(1:nvec,col_idx(j))
    end do
    do j = halo_ptr(i), row_ptr(i+1)-1, 1
      tmp_r(:) = tmp_r(:) - val(j)*halo_r(:,col_idx(j))
      tmp_i(:) = tmp_i(:) - val(j)*halo_i(:,col_idx(j))
    end do
    x_r(1:nvec,i) = rhs(i)-tmp_r(:)*nrmsai2i
    x_i(1:nvec,i) =       -tmp_i(:)*nrmsai2i
  end do
end subroutine dkacz_generic



subroutine dkacz_generic(nvec, nlocal, nhalo, ncols, nnz, row_ptr, halo_ptr, 
col_idx, val, shifts_r,shifts_i, x_r,x_i, ldx, halo_r, halo_i,nrmsai2i,istart,iend,istep)
  implicit none
  integer, intent(in) :: nvec, nlocal, nhalo, ncols, ldx
  integer(kind=8), intent(in) :: nnz
  c(kind=8), intent(in) :: shift_r, shift_i
  integer(kind=8), intent(in) :: row_ptr(nlocal+1), halo_ptr(nlocal)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(inout) :: x_r(ldx,*), x_i(ldx,*)
  real(kind=8), intend(inout) :: halo_r(nvec,nhalo),halo_i(nvec,nhalo)
  real(kind=8), intent(in) :: nrmsai2(nlocal)
  integer, intent(in) :: istart,iend,istep
  real(kind=8) :: tmp_r(nvec), tmp_i(nvec)
  integer :: i
  integer(kind=8) :: j

!TODO $omp parallel do private(tmp) schedule(static)
  do i = istart, iend,istep
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! compute (sigma_j I - A)_i*x
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    tmp_r(1:nvec) = shifts_r(1:nvec)*x_r(1:nvec,i) - &
                    shifts_i(1:nvec)*x_i(1:nvec,i)
    tmp_i(1:nvec) = shifts_i(1:nvec)*x_r(1:nvec,i) + &
                    shifts_r(1:nvec)*x_i(1:nvec,i)
    do j = row_ptr(i), halo_ptr(i)-1, 1
      tmp_r(:) = tmp_r(:) - val(j)*x_r(1:nvec,col_idx(j))
      tmp_i(:) = tmp_i(:) - val(j)*x_i(1:nvec,col_idx(j))
    end do
    do j = halo_ptr(i), row_ptr(i+1)-1, 1
      tmp_r(:) = tmp_r(:) - val(j)*halo_r(:,col_idx(j))
      tmp_i(:) = tmp_i(:) - val(j)*halo_i(:,col_idx(j))
    end do
    x_r(1:nvec,i) = rhs(i)-tmp_r(:)*nrmsai2i
    x_i(1:nvec,i) =       -tmp_i(:)*nrmsai2i
  end do
end subroutine dkacz_generic


