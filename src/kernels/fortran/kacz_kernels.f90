!! compute for a given shift sigma=shift_r+shift_i the 2-norm of each row
!! of the matrix shift[j]*I-A and store the inverse of the results in the
!! columns j of the block vector nrms_ai2i. On entry, nrms_ai2i(:,1) should
!! contain the diagonal elements of A, aii, to make things easier here (!).
subroutine crsmat_norms_ai2i(nshifts, nlocal, nnz, row_ptr, &
         val, shifts_r,shifts_i, nrms_ai2i)
  implicit none
  integer, intent(in) :: nshifts, nlocal
  integer(kind=8), intent(in) :: nnz
  real(kind=8), intent(in) :: shifts_r(nshifts), shifts_i(nshifts)
  integer(kind=8), intent(in) :: row_ptr(nlocal+1)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(out) :: nrms_ai2i(nlocal,nshifts)
  ! local variables
  integer :: i
  integer(kind=8) :: j
  real(kind=8) :: tmp

!$omp parallel do private(tmp) schedule(static)
  do i = 1,nlocal
    ! for off-diagonal elements, add aij^2 to tmp.
    ! for diagonal element, the term we need is
    ! (aii-s[j])^2 = (aii-(sr[j]+i*si[j]))(aii-(sr[j]-i*si[j]))
    !              = aii^2-2*aii*sr+sr^2+si^2. 
    ! aii^2 is added in the regular loop below.
    tmp=nrms_ai2i(i,1) ! =aii, see comment in the function header
    nrms_ai2i(i,:)=shifts_r(1:nshifts)*shifts_r(1:nshifts) + &
                   shifts_i(1:nshifts)*shifts_i(1:nshifts) - &
                   2.d0*tmp*shifts_r(1:nshifts)
    tmp=0.d0
    do j = row_ptr(i), row_ptr(i+1)-1, 1
      tmp = tmp + val(j)*val(j)
    end do
    nrms_ai2i(i,:)=1.d0/(nrms_ai2i(i,:)+tmp)
  end do
end subroutine crsmat_norms_ai2i

!! general implementation of forward or backward Kaczmarz sweep for a single shift
!! and possibly multiple vector columns in X and B
subroutine dkacz_generic(nvec, nlocal, nhalo, ncols, nnz, &
row_ptr, halo_ptr, col_idx, val, &
shift_r,shift_i, b, ldb, &
x_r,x_i, ldx, halo_r, halo_i,nrms_ai2i,omega,istart,iend,istep)
  implicit none
  integer, intent(in) :: nvec, nlocal, nhalo, ncols, ldx, ldb
  integer(kind=8), intent(in) :: nnz
  real(kind=8), intent(in) :: shift_r, shift_i
  integer(kind=8), intent(in) :: row_ptr(nlocal+1), halo_ptr(nlocal)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  real(kind=8), intent(inout) :: x_r(ldx,*), x_i(ldx,*),b(ldb,*)
  real(kind=8), intent(inout) :: halo_r(nvec,nhalo),halo_i(nvec,nhalo)
  real(kind=8), intent(in) :: nrms_ai2i(nlocal)
  real(kind=8), intent(in) :: omega
  integer, intent(in) :: istart,iend,istep
  real(kind=8) :: tmp_r(nvec), tmp_i(nvec)
  integer :: i
  integer(kind=8) :: j

!TODO - OpenMP, coloring
  do i = istart, iend,istep
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! compute (shift_j I - A)_i*x
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    tmp_r(1:nvec) = shift_r*x_r(1:nvec,i) - &
                    shift_i*x_i(1:nvec,i)
    tmp_i(1:nvec) = shift_i*x_r(1:nvec,i) + &
                    shift_r*x_i(1:nvec,i)
    do j = row_ptr(i), halo_ptr(i)-1, 1
      tmp_r(:) = tmp_r(:) - val(j)*x_r(1:nvec,col_idx(j))
      tmp_i(:) = tmp_i(:) - val(j)*x_i(1:nvec,col_idx(j))
    end do
    do j = halo_ptr(i), row_ptr(i+1)-1, 1
      tmp_r(:) = tmp_r(:) - val(j)*halo_r(1:nvec,col_idx(j))
      tmp_i(:) = tmp_i(:) - val(j)*halo_i(1:nvec,col_idx(j))
    end do
    ! Kaczmarz update of X

    ! a) scaling factors
    tmp_r(:)=tmp_r(:)*omega*nrms_ai2i(i)
    tmp_i(:)=tmp_i(:)*omega*nrms_ai2i(i)

    ! b) projection step
    x_r(1:nvec,i)=x_r(1:nvec,i) + (tmp_r(:)*shift_r+tmp_i(:)*shift_i)
    x_i(1:nvec,i)=x_i(1:nvec,i) + (tmp_i(:)*shift_r-tmp_r(:)*shift_i)

    do j = row_ptr(i), halo_ptr(i)-1, 1
      x_r(1:nvec,col_idx(j)) = x_r(1:nvec,col_idx(j)) - &
                               tmp_r(:)*val(j)
      x_i(1:nvec,col_idx(j)) = x_i(1:nvec,col_idx(j)) - &
                               tmp_i(:)*val(j)
    end do
    do j = halo_ptr(i), row_ptr(i+1)-1, 1
      halo_r(1:nvec,col_idx(j)) = halo_r(1:nvec,col_idx(j)) - &
                               tmp_r(:)*val(j)
      halo_i(1:nvec,col_idx(j)) = halo_i(1:nvec,col_idx(j)) - &
                               tmp_i(:)*val(j)
    end do
  end do
end subroutine dkacz_generic


