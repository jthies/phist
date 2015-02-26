#include "phist_config.h"

!! pass all subroutine args to specialized kernels:
#define _KACZ_LOOP_ARGS_RC_ nvec,nlocal,nhalo,ncols,nnz,row_ptr,halo_ptr,col_idx,val,map,\
shift_r,shift_i,b,ldb,x_r,x_i,ldx,halo_r,halo_i,nrms_ai2i,omega,istart,iend,istep,i0,i1,j0,j1

!! pass dummy args for imaginary parts in case of real matrix and real shifts
#define _KACZ_LOOP_ARGS_  nvec,nlocal,nhalo,ncols,nnz,row_ptr,halo_ptr,col_idx,val,map,\
shift_r,shift_r,b,ldb,x_r,x_r,ldx,halo_r,halo_r,nrms_ai2i,omega,istart,iend,istep,i0,i1,j0,j1

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
  do i = 1,nlocal, 1
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

!! variant with real matrix and rhs but complex shift and x vector
#define KACZ_RC_VARIANT

!! general implementation of forward or backward Kaczmarz sweep for a single shift
!! and possibly multiple vector columns in X and B
subroutine dkacz_selector(nvec, nlocal, nhalo, ncols, nnz, &
        row_ptr, halo_ptr, col_idx, val, map, &
        shift_r,shift_i, b, ldb, &
        x_r,x_i, ldx, halo_r, halo_i,nrms_ai2i,omega,istart_in,iend_in,istep)
#ifdef PHIST_HAVE_OPENMP
  use :: omp_lib
#endif
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
  integer, intent(in) :: istart_in,iend_in,istep
  ! locals
  real(kind=8) :: tmp_r(nvec), tmp_i(nvec)
  integer :: i, ic, jc, i0, i1, j0, j1, istart, iend
  integer(kind=8) :: j
  logical :: use_clr_kernel
  logical :: use_bzero_kernel
  integer istart_clr, iend_clr

  ! if there is a dist-2 coloring available, call
  ! an alternative kernel.
  use_clr_kernel= (map%nColors>0) .and. &
                  allocated(map%color_offset) .and. &
                  allocated(map%color_idx) .and. &
                  (map%coloringType==2) .and. &
                  ((istep==1) .or. (istep==-1))

  use_bzero_kernel=(ldb==0)
  
  if (use_clr_kernel) then
    if (istep==1) then
      istart = 1
      iend= map%nColors
      j0=0
      j1=-1
    else if (istep==-1) then
      istart= map%nColors+1
      iend = 2
      j0=-1
      j1=0
    end if
  else
    istart=istart_in
    iend=iend_in
  end if



  if (use_clr_kernel) then
    if (use_bzero_kernel) then
      if (nvec==1) then
        call kacz_loop_RC_clr_b0_1(_KACZ_LOOP_ARGS_RC_)
      else if (nvec==2) then
        call kacz_loop_RC_clr_b0_2(_KACZ_LOOP_ARGS_RC_)
      else if (nvec==4) then
        call kacz_loop_RC_clr_b0_4(_KACZ_LOOP_ARGS_RC_)
      else if (nvec==8) then
        call kacz_loop_RC_clr_b0_8(_KACZ_LOOP_ARGS_RC_)
      else if (nvec==16) then
        call kacz_loop_RC_clr_b0_16(_KACZ_LOOP_ARGS_RC_)
      else if (nvec==32) then
        call kacz_loop_RC_clr_b0_32(_KACZ_LOOP_ARGS_RC_)
      else if (nvec==64) then
        call kacz_loop_RC_clr_b0_64(_KACZ_LOOP_ARGS_RC_)
      else
        call kacz_loop_RC_clr_b0_generic(_KACZ_LOOP_ARGS_RC_)
      end if ! block size
    else
      if (nvec==1) then
        call kacz_loop_RC_clr_1(_KACZ_LOOP_ARGS_RC_)
      else if (nvec==2) then
        call kacz_loop_RC_clr_2(_KACZ_LOOP_ARGS_RC_)
      else if (nvec==4) then
        call kacz_loop_RC_clr_4(_KACZ_LOOP_ARGS_RC_)
      else if (nvec==8) then
        call kacz_loop_RC_clr_8(_KACZ_LOOP_ARGS_RC_)
      else if (nvec==16) then
        call kacz_loop_RC_clr_16(_KACZ_LOOP_ARGS_RC_)
      else if (nvec==32) then
        call kacz_loop_RC_clr_32(_KACZ_LOOP_ARGS_RC_)
      else if (nvec==64) then
        call kacz_loop_RC_clr_64(_KACZ_LOOP_ARGS_RC_)
      else
        call kacz_loop_RC_clr_generic(_KACZ_LOOP_ARGS_RC_)
      end if ! block size
    end if ! rhs=0
  else ! no coloring - use lexicographic kernel
    if (use_bzero_kernel) then
      if (nvec==1) then
        call kacz_loop_RC_seq_b0_1(_KACZ_LOOP_ARGS_RC_)
      else if (nvec==2) then
        call kacz_loop_RC_seq_b0_2(_KACZ_LOOP_ARGS_RC_)
      else if (nvec==4) then
        call kacz_loop_RC_seq_b0_4(_KACZ_LOOP_ARGS_RC_)
      else if (nvec==8) then
        call kacz_loop_RC_seq_b0_8(_KACZ_LOOP_ARGS_RC_)
      else if (nvec==16) then
        call kacz_loop_RC_seq_b0_16(_KACZ_LOOP_ARGS_RC_)
      else if (nvec==32) then
        call kacz_loop_RC_seq_b0_32(_KACZ_LOOP_ARGS_RC_)
      else if (nvec==64) then
        call kacz_loop_RC_seq_b0_64(_KACZ_LOOP_ARGS_RC_)
      else
        call kacz_loop_RC_seq_b0_generic(_KACZ_LOOP_ARGS_RC_)
      end if ! block size
    else
      if (nvec==1) then
        call kacz_loop_RC_seq_1(_KACZ_LOOP_ARGS_RC_)
      else if (nvec==2) then
        call kacz_loop_RC_seq_2(_KACZ_LOOP_ARGS_RC_)
      else if (nvec==4) then
        call kacz_loop_RC_seq_4(_KACZ_LOOP_ARGS_RC_)
      else if (nvec==8) then
        call kacz_loop_RC_seq_8(_KACZ_LOOP_ARGS_RC_)
      else if (nvec==16) then
        call kacz_loop_RC_seq_16(_KACZ_LOOP_ARGS_RC_)
      else if (nvec==32) then
        call kacz_loop_RC_seq_32(_KACZ_LOOP_ARGS_RC_)
      else if (nvec==64) then
        call kacz_loop_RC_seq_64(_KACZ_LOOP_ARGS_RC_)
      else
        call kacz_loop_RC_seq_generic(_KACZ_LOOP_ARGS_RC_)
      end if ! block size
    end if ! rhs=0
  end if ! use coloring for parallelization

end subroutine dkacz_selector

! TODO - is this still an issue?
#if defined(__INTEL_COMPILER)&&defined(TESTING)&&defined(PHIST_HAVE_COLPACK)
#warning "certain ifort debugging flags ('-check all', I think) break the real-only kacz kernel with coloring"
#endif

!! general implementation of forward or backward Kaczmarz sweep for a single shift
!! and possibly multiple vector columns in X and B
subroutine dkacz_selector_real(nvec, nlocal, nhalo, ncols, nnz, &
        row_ptr, halo_ptr, col_idx, val, map, &
        shift_r, b, ldb, &
        x_r, ldx, halo_r, nrms_ai2i,omega,istart_in,iend_in,istep)

  use :: omp_lib
  use :: map_module

  implicit none

  integer, intent(in) :: nvec, nlocal, nhalo, ncols, ldx, ldb
  integer(kind=8), intent(in) :: nnz
  real(kind=8), intent(in) :: shift_r
  integer(kind=8), intent(in) :: row_ptr(nlocal+1), halo_ptr(nlocal)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  TYPE(Map_t), intent(in) :: map
  real(kind=8), intent(inout) :: x_r(ldx,*), b(ldb,*)
  real(kind=8), intent(inout) :: halo_r(nvec,nhalo)
  real(kind=8), intent(in) :: nrms_ai2i(nlocal)
  real(kind=8), intent(in) :: omega
  integer, intent(in) :: istart_in,iend_in,istep
  ! locals
  real(kind=8) :: tmp_r(nvec),tmp_i(nvec)
  integer :: i, ic, jc, i0, i1, j0, j1, istart, iend
  integer(kind=8) :: j
  logical :: use_clr_kernel
  logical :: use_bzero_kernel
  integer istart_clr, iend_clr

  ! if there is a dist-2 coloring available, call
  ! an alternative kernel.
  use_clr_kernel= (map%nColors>0) .and. &
                  allocated(map%color_offset) .and. &
                  allocated(map%color_idx) .and. &
                  (map%coloringType==2) .and. &
                  ((istep==1) .or. (istep==-1))

  use_bzero_kernel=(ldb==0)

  if (use_clr_kernel) then
    if (istep==1) then
      istart = 1
      iend= map%nColors
      j0=0
      j1=-1
    else if (istep==-1) then
      istart= map%nColors+1
      iend = 2
      j0=-1
      j1=0
    end if
  else
    istart=istart_in
    iend=iend_in
  end if

  if (use_clr_kernel) then
    if (use_bzero_kernel) then
      if (nvec==1) then
        call kacz_loop_clr_b0_1(_KACZ_LOOP_ARGS_)
      else if (nvec==2) then
        call kacz_loop_clr_b0_2(_KACZ_LOOP_ARGS_)
      else if (nvec==4) then
        call kacz_loop_clr_b0_4(_KACZ_LOOP_ARGS_)
      else if (nvec==8) then
        call kacz_loop_clr_b0_8(_KACZ_LOOP_ARGS_)
      else if (nvec==16) then
        call kacz_loop_clr_b0_16(_KACZ_LOOP_ARGS_)
      else if (nvec==32) then
        call kacz_loop_clr_b0_32(_KACZ_LOOP_ARGS_)
      else if (nvec==64) then
        call kacz_loop_clr_b0_64(_KACZ_LOOP_ARGS_)
      else
        call kacz_loop_clr_b0_generic(_KACZ_LOOP_ARGS_)
      end if ! block size
    else
      if (nvec==1) then
        call kacz_loop_clr_1(_KACZ_LOOP_ARGS_)
      else if (nvec==2) then
        call kacz_loop_clr_2(_KACZ_LOOP_ARGS_)
      else if (nvec==4) then
        call kacz_loop_clr_4(_KACZ_LOOP_ARGS_)
      else if (nvec==8) then
        call kacz_loop_clr_8(_KACZ_LOOP_ARGS_)
      else if (nvec==16) then
        call kacz_loop_clr_16(_KACZ_LOOP_ARGS_)
      else if (nvec==32) then
        call kacz_loop_clr_32(_KACZ_LOOP_ARGS_)
      else if (nvec==64) then
        call kacz_loop_clr_64(_KACZ_LOOP_ARGS_)
      else
        call kacz_loop_clr_generic(_KACZ_LOOP_ARGS_)
      end if ! block size
    end if ! rhs=0
  else ! no coloring - use lexicographic kernel
    if (use_bzero_kernel) then
      if (nvec==1) then
        call kacz_loop_seq_b0_1(_KACZ_LOOP_ARGS_)
      else if (nvec==2) then
        call kacz_loop_seq_b0_2(_KACZ_LOOP_ARGS_)
      else if (nvec==4) then
        call kacz_loop_seq_b0_4(_KACZ_LOOP_ARGS_)
      else if (nvec==8) then
        call kacz_loop_seq_b0_8(_KACZ_LOOP_ARGS_)
      else if (nvec==16) then
        call kacz_loop_seq_b0_16(_KACZ_LOOP_ARGS_)
      else if (nvec==32) then
        call kacz_loop_seq_b0_32(_KACZ_LOOP_ARGS_)
      else if (nvec==64) then
        call kacz_loop_seq_b0_64(_KACZ_LOOP_ARGS_)
      else
        call kacz_loop_seq_b0_generic(_KACZ_LOOP_ARGS_)
      end if ! block size
    else
      if (nvec==1) then
        call kacz_loop_seq_1(_KACZ_LOOP_ARGS_)
      else if (nvec==2) then
        call kacz_loop_seq_2(_KACZ_LOOP_ARGS_)
      else if (nvec==4) then
        call kacz_loop_seq_4(_KACZ_LOOP_ARGS_)
      else if (nvec==8) then
        call kacz_loop_seq_8(_KACZ_LOOP_ARGS_)
      else if (nvec==16) then
        call kacz_loop_seq_16(_KACZ_LOOP_ARGS_)
      else if (nvec==32) then
        call kacz_loop_seq_32(_KACZ_LOOP_ARGS_)
      else if (nvec==64) then
        call kacz_loop_seq_64(_KACZ_LOOP_ARGS_)
      else
        call kacz_loop_seq_generic(_KACZ_LOOP_ARGS_)
      end if ! block size
    end if ! rhs=0
  end if ! use coloring for parallelization



end subroutine dkacz_selector_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! special subroutines (experimental, fo rbetter performance)    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! REAL-ONLY                     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define KACZ_CLR
#define KACZ_BZERO

#define NVEC 1
#define SUB_NAME kacz_loop_clr_b0_1
#include "kacz_sub_def.h"

#define NVEC 2
#define SUB_NAME kacz_loop_clr_b0_2
#include "kacz_sub_def.h"

#define NVEC 4
#define SUB_NAME kacz_loop_clr_b0_4
#include "kacz_sub_def.h"

#define NVEC 8
#define SUB_NAME kacz_loop_clr_b0_8
#include "kacz_sub_def.h"

#define NVEC 16
#define SUB_NAME kacz_loop_clr_b0_16
#include "kacz_sub_def.h"

#define NVEC 32
#define SUB_NAME kacz_loop_clr_b0_32
#include "kacz_sub_def.h"

#define NVEC 64
#define SUB_NAME kacz_loop_clr_b0_64
#include "kacz_sub_def.h"

! NVEC is now not defined, so it will be
! taken from the first arg to the subroutine
#define SUB_NAME kacz_loop_clr_b0_generic
#include "kacz_sub_def.h"

#undef KACZ_BZERO

#define NVEC 1
#define SUB_NAME kacz_loop_clr_1
#include "kacz_sub_def.h"

#define NVEC 2
#define SUB_NAME kacz_loop_clr_2
#include "kacz_sub_def.h"

#define NVEC 4
#define SUB_NAME kacz_loop_clr_4
#include "kacz_sub_def.h"

#define NVEC 8
#define SUB_NAME kacz_loop_clr_8
#include "kacz_sub_def.h"

#define NVEC 16
#define SUB_NAME kacz_loop_clr_16
#include "kacz_sub_def.h"

#define NVEC 32
#define SUB_NAME kacz_loop_clr_32
#include "kacz_sub_def.h"

#define NVEC 64
#define SUB_NAME kacz_loop_clr_64
#include "kacz_sub_def.h"

#define SUB_NAME kacz_loop_clr_generic
#include "kacz_sub_def.h"

#undef KACZ_CLR
#define KACZ_BZERO

#define NVEC 1
#define SUB_NAME kacz_loop_seq_b0_1
#include "kacz_sub_def.h"

#define NVEC 2
#define SUB_NAME kacz_loop_seq_b0_2
#include "kacz_sub_def.h"

#define NVEC 4
#define SUB_NAME kacz_loop_seq_b0_4
#include "kacz_sub_def.h"

#define NVEC 8
#define SUB_NAME kacz_loop_seq_b0_8
#include "kacz_sub_def.h"

#define NVEC 16
#define SUB_NAME kacz_loop_seq_b0_16
#include "kacz_sub_def.h"

#define NVEC 32
#define SUB_NAME kacz_loop_seq_b0_32
#include "kacz_sub_def.h"

#define NVEC 64
#define SUB_NAME kacz_loop_seq_b0_64
#include "kacz_sub_def.h"

#define SUB_NAME kacz_loop_seq_b0_generic
#include "kacz_sub_def.h"

#undef KACZ_BZERO

#define NVEC 1
#define SUB_NAME kacz_loop_seq_1
#include "kacz_sub_def.h"

#define NVEC 2
#define SUB_NAME kacz_loop_seq_2
#include "kacz_sub_def.h"

#define NVEC 4
#define SUB_NAME kacz_loop_seq_4
#include "kacz_sub_def.h"

#define NVEC 8
#define SUB_NAME kacz_loop_seq_8
#include "kacz_sub_def.h"

#define NVEC 16
#define SUB_NAME kacz_loop_seq_16
#include "kacz_sub_def.h"

#define NVEC 32
#define SUB_NAME kacz_loop_seq_32
#include "kacz_sub_def.h"

#define NVEC 64
#define SUB_NAME kacz_loop_seq_64
#include "kacz_sub_def.h"

#define SUB_NAME kacz_loop_seq_generic
#include "kacz_sub_def.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! REAL with complex shifts      !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define KACZ_RC_VARIANT
#define KACZ_CLR
#define KACZ_BZERO

#define NVEC 1
#define SUB_NAME kacz_loop_RC_clr_b0_1
#include "kacz_sub_def.h"

#define NVEC 2
#define SUB_NAME kacz_loop_RC_clr_b0_2
#include "kacz_sub_def.h"

#define NVEC 4
#define SUB_NAME kacz_loop_RC_clr_b0_4
#include "kacz_sub_def.h"

#define NVEC 8
#define SUB_NAME kacz_loop_RC_clr_b0_8
#include "kacz_sub_def.h"

#define NVEC 16
#define SUB_NAME kacz_loop_RC_clr_b0_16
#include "kacz_sub_def.h"

#define NVEC 32
#define SUB_NAME kacz_loop_RC_clr_b0_32
#include "kacz_sub_def.h"

#define NVEC 64
#define SUB_NAME kacz_loop_RC_clr_b0_64
#include "kacz_sub_def.h"

#undef KACZ_BZERO

#define NVEC 1
#define SUB_NAME kacz_loop_RC_clr_1
#include "kacz_sub_def.h"

#define NVEC 2
#define SUB_NAME kacz_loop_RC_clr_2
#include "kacz_sub_def.h"

#define NVEC 4
#define SUB_NAME kacz_loop_RC_clr_4
#include "kacz_sub_def.h"

#define NVEC 8
#define SUB_NAME kacz_loop_RC_clr_8
#include "kacz_sub_def.h"

#define NVEC 16
#define SUB_NAME kacz_loop_RC_clr_16
#include "kacz_sub_def.h"

#define NVEC 32
#define SUB_NAME kacz_loop_RC_clr_32
#include "kacz_sub_def.h"

#define NVEC 64
#define SUB_NAME kacz_loop_RC_clr_64
#include "kacz_sub_def.h"

#define SUB_NAME kacz_loop_RC_clr_generic
#include "kacz_sub_def.h"

#undef KACZ_CLR
#define KACZ_BZERO

#define NVEC 1
#define SUB_NAME kacz_loop_RC_seq_b0_1
#include "kacz_sub_def.h"

#define NVEC 2
#define SUB_NAME kacz_loop_RC_seq_b0_2
#include "kacz_sub_def.h"

#define NVEC 4
#define SUB_NAME kacz_loop_RC_seq_b0_4
#include "kacz_sub_def.h"

#define NVEC 8
#define SUB_NAME kacz_loop_RC_seq_b0_8
#include "kacz_sub_def.h"

#define NVEC 16
#define SUB_NAME kacz_loop_RC_seq_b0_16
#include "kacz_sub_def.h"

#define NVEC 32
#define SUB_NAME kacz_loop_RC_seq_b0_32
#include "kacz_sub_def.h"

#define NVEC 64
#define SUB_NAME kacz_loop_RC_seq_b0_64
#include "kacz_sub_def.h"

#define SUB_NAME kacz_loop_RC_seq_b0_generic
#include "kacz_sub_def.h"

#undef KACZ_BZERO

#define NVEC 1
#define SUB_NAME kacz_loop_RC_seq_1
#include "kacz_sub_def.h"

#define NVEC 2
#define SUB_NAME kacz_loop_RC_seq_2
#include "kacz_sub_def.h"

#define NVEC 4
#define SUB_NAME kacz_loop_RC_seq_4
#include "kacz_sub_def.h"

#define NVEC 8
#define SUB_NAME kacz_loop_RC_seq_8
#include "kacz_sub_def.h"

#define NVEC 16
#define SUB_NAME kacz_loop_RC_seq_16
#include "kacz_sub_def.h"

#define NVEC 32
#define SUB_NAME kacz_loop_RC_seq_32
#include "kacz_sub_def.h"

#define NVEC 64
#define SUB_NAME kacz_loop_RC_seq_64
#include "kacz_sub_def.h"

#define SUB_NAME kacz_loop_RC_seq_generic
#include "kacz_sub_def.h"


