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
        row_ptr, halo_ptr, col_idx, val, map, &
        shift_r,shift_i, b, ldb, &
        x_r,x_i, ldx, halo_r, halo_i,nrms_ai2i,omega,istart,iend,istep)

  use :: omp_lib
  use :: map_module

  implicit none

  integer, intent(in) :: nvec, nlocal, nhalo, ncols, ldx, ldb
  integer(kind=8), intent(in) :: nnz
  real(kind=8), intent(in) :: shift_r, shift_i
  integer(kind=8), intent(in) :: row_ptr(nlocal+1), halo_ptr(nlocal)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  TYPE(Map_t), intent(in) :: map
  real(kind=8), intent(inout) :: x_r(ldx,*), x_i(ldx,*),b(ldb,*)
  real(kind=8), intent(inout) :: halo_r(nvec,nhalo),halo_i(nvec,nhalo)
  real(kind=8), intent(in) :: nrms_ai2i(nlocal)
  real(kind=8), intent(in) :: omega
  integer, intent(in) :: istart,iend,istep
  ! locals
  real(kind=8) :: tmp_r(nvec), tmp_i(nvec)
  integer :: i
  integer(kind=8) :: j
  logical :: use_clr_kernel
  integer istart_clr, iend_clr

  ! if there is a dist-2 coloring available, call
  ! an alternative kernel.
  use_clr_kernel= (map%nColors>0) .and. &
                  allocated(map%color_offset) .and. &
                  allocated(map%color_idx) .and. &
                  (map%coloringType==2) .and. &
                  ((istep==1) .or. (istep==-1))

  write(*,*) 'use_clr_kernel=',use_clr_kernel
  write(*,*) 'map%nColors=',map%nColors
  write(*,*) 'map%coloringType=',map%coloringType
  write(*,*) 'istep=',istep
  write(*,*) 'istart=',istart
  write(*,*) 'iend=',iend

  if (use_clr_kernel) then
    if (istep==1) then
      istart_clr = 1
      iend_clr= map%nColors
    else if (istep==-1) then
      istart_clr= map%nColors+1
      iend_clr = 2
    end if
    write(*,*) 'calling generic OpenMP/Coloring kernel'
    call dkacz_generic_clr(nvec, nlocal, nhalo, ncols, nnz, &
        row_ptr, halo_ptr, col_idx, val, map, &
        shift_r,shift_i, b, ldb, &
        x_r,x_i, ldx, halo_r, halo_i,nrms_ai2i,omega,istart_clr,iend_clr,istep)
    return  
  end if

  ! sequential implementation ignoring coloring information
  ! (lexicographic or given ordering)
  write(*,*) 'using generic sequential/lexicographic kernel'

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

    tmp_r(:)=tmp_r(:)-b(1:nvec,i)

    ! Kaczmarz update of X

    ! a) scaling factors
    tmp_r(:)=tmp_r(:)*omega*nrms_ai2i(i)
    tmp_i(:)=tmp_i(:)*omega*nrms_ai2i(i)

    ! b) projection step
    x_r(1:nvec,i)=x_r(1:nvec,i) - (tmp_r(:)*shift_r+tmp_i(:)*shift_i)
    x_i(1:nvec,i)=x_i(1:nvec,i) - (tmp_i(:)*shift_r-tmp_r(:)*shift_i)

    do j = row_ptr(i), halo_ptr(i)-1, 1
      x_r(1:nvec,col_idx(j)) = x_r(1:nvec,col_idx(j)) + &
                               tmp_r(:)*val(j)
      x_i(1:nvec,col_idx(j)) = x_i(1:nvec,col_idx(j)) + &
                               tmp_i(:)*val(j)
    end do
    do j = halo_ptr(i), row_ptr(i+1)-1, 1
      halo_r(1:nvec,col_idx(j)) = halo_r(1:nvec,col_idx(j)) + &
                               tmp_r(:)*val(j)
      halo_i(1:nvec,col_idx(j)) = halo_i(1:nvec,col_idx(j)) + &
                               tmp_i(:)*val(j)
    end do
  end do

end subroutine dkacz_generic

!! variant of kacz_genric kernel that exploits coloring information
!! in the map to generate (OpenMP) parallelism
subroutine dkacz_generic_clr(nvec, nlocal, nhalo, ncols, nnz, &
        row_ptr, halo_ptr, col_idx, val, map, &
        shift_r,shift_i, b, ldb, &
        x_r,x_i, ldx, halo_r, halo_i,nrms_ai2i,omega,istart,iend,istep)

  use :: map_module

  implicit none

  integer, intent(in) :: nvec, nlocal, nhalo, ncols, ldx, ldb
  integer(kind=8), intent(in) :: nnz
  real(kind=8), intent(in) :: shift_r, shift_i
  integer(kind=8), intent(in) :: row_ptr(nlocal+1), halo_ptr(nlocal)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  TYPE(Map_t), intent(in) :: map
  real(kind=8), intent(inout) :: x_r(ldx,*), x_i(ldx,*),b(ldb,*)
  real(kind=8), intent(inout) :: halo_r(nvec,nhalo),halo_i(nvec,nhalo)
  real(kind=8), intent(in) :: nrms_ai2i(nlocal)
  real(kind=8), intent(in) :: omega
  integer, intent(in) :: istart,iend,istep
  ! locals
  real(kind=8) :: tmp_r(nvec), tmp_i(nvec)
  integer :: i, ic, jc,j0,j1
  integer(kind=8) :: j
  
  if (istep==1) then
    j0=0
    j1=-1
  else if (istep==-1) then
    j0=-1
    j1=0
  end if

  ! shared-memory parallel implementation using coloring info in the map
  do ic = istart, iend,istep
! note: it doesn't make sense to try and use the same scheduling as elsewhere
! (to avoid NUMA problems) because the ordering here is diffrent. For now, we
! will stick with 1 MPI process per socket and think about NUMA later (we could
! 'first touch' all vectors using the coloring if it is defined in the map, but
! that would infringe the spMVM performance...)
!$omp parallel do private(tmp_r,tmp_i,i,j) schedule(runtime)
  do jc = map%color_offset(ic)+j0,map%color_offset(ic+istep)+j1,istep
    i=map%color_idx(jc)    
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

    tmp_r(:)=tmp_r(:)-b(1:nvec,i)

    ! Kaczmarz update of X

    ! a) scaling factors
    tmp_r(:)=tmp_r(:)*omega*nrms_ai2i(i)
    tmp_i(:)=tmp_i(:)*omega*nrms_ai2i(i)

    ! b) projection step
    x_r(1:nvec,i)=x_r(1:nvec,i) - (tmp_r(:)*shift_r+tmp_i(:)*shift_i)
    x_i(1:nvec,i)=x_i(1:nvec,i) - (tmp_i(:)*shift_r-tmp_r(:)*shift_i)

    do j = row_ptr(i), halo_ptr(i)-1, 1
      x_r(1:nvec,col_idx(j)) = x_r(1:nvec,col_idx(j)) + &
                               tmp_r(:)*val(j)
      x_i(1:nvec,col_idx(j)) = x_i(1:nvec,col_idx(j)) + &
                               tmp_i(:)*val(j)
    end do
    do j = halo_ptr(i), row_ptr(i+1)-1, 1
      halo_r(1:nvec,col_idx(j)) = halo_r(1:nvec,col_idx(j)) + &
                               tmp_r(:)*val(j)
      halo_i(1:nvec,col_idx(j)) = halo_i(1:nvec,col_idx(j)) + &
                               tmp_i(:)*val(j)
    end do
  end do
  end do
end subroutine dkacz_generic_clr

!! implementation of kacz_generic for the case b=0. Once kacz_generic
!! is stable, we could copy/paste/adjust it into this function rather
!! than creating an empty b array and calling the other kernel.
!!
!! We want a specialized kernel for this situation (b=0) as it occurs
!! in every iteration of CGNM/CARP-CG except on the initial call.
subroutine dkacz_bzero_generic(nvec, nlocal, nhalo, ncols, nnz, &
row_ptr, halo_ptr, col_idx, val, map, &
shift_r,shift_i,&
x_r,x_i, ldx, halo_r, halo_i,nrms_ai2i,omega,istart,iend,istep)
  use :: map_module
  implicit none
  integer, intent(in) :: nvec, nlocal, nhalo, ncols, ldx
  integer(kind=8), intent(in) :: nnz
  real(kind=8), intent(in) :: shift_r, shift_i
  integer(kind=8), intent(in) :: row_ptr(nlocal+1), halo_ptr(nlocal)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  TYPE(Map_t), intent(in) :: map
  real(kind=8), intent(inout) :: x_r(ldx,*), x_i(ldx,*)
  real(kind=8), intent(inout) :: halo_r(nvec,nhalo),halo_i(nvec,nhalo)
  real(kind=8), intent(in) :: nrms_ai2i(nlocal)
  real(kind=8), intent(in) :: omega
  integer, intent(in) :: istart,iend,istep
  !-------------------------------------------------------------------!
  real(kind=8) :: tmp_r(nvec), tmp_i(nvec)
  integer :: i
  integer(kind=8) :: j
  real(kind=8) :: b(nvec,nlocal)

!write(*,*) 'ldx=',ldx

  b(:,:)=0.d0

  call dkacz_generic(nvec, nlocal, nhalo, ncols, nnz, &
row_ptr, halo_ptr, col_idx, val, map, &
shift_r,shift_i, b, nvec, &
x_r,x_i, ldx, halo_r, halo_i,nrms_ai2i,omega,istart,iend,istep)

end subroutine dkacz_bzero_generic
