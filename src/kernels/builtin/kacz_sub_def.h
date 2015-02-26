subroutine SUB_NAME(nvec,nlocal, nhalo, ncols, nnz, &
        row_ptr, halo_ptr, col_idx, val, map, &
        shift_r,shift_i, b, ldb, &
        x_r,x_i, ldx, halo_r, halo_i,nrms_ai2i,omega,&
        istart,iend,istep,i0,i1,j0,j1)

#ifdef PHIST_HAVE_OPENMP
  use :: omp_lib
#endif
  use :: map_module

  implicit none

  integer, intent(in) :: nlocal, nhalo, ncols, ldx, ldb
  integer(kind=8), intent(in) :: nnz
  real(kind=8), intent(in) :: shift_r, shift_i
  integer(kind=8), intent(in) :: row_ptr(nlocal+1), halo_ptr(nlocal)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  TYPE(Map_t), intent(in) :: map
  real(kind=8), intent(inout) :: x_r(ldx,*), x_i(ldx,*),b(ldb,*)
  real(kind=8), intent(inout) :: halo_r(NVEC,nhalo),halo_i(NVEC,nhalo)
  real(kind=8), intent(in) :: nrms_ai2i(nlocal)
  real(kind=8), intent(in) :: omega
  integer, intent(in) :: i0, i1, j0, j1, istart, iend
  ! locals
  real(kind=8) :: tmp_r(NVEC), tmp_i(NVEC)
  integer :: i, ic, jc
  integer(kind=8) :: j
  integer istart_clr, iend_clr

#ifdef KACZ_CLR
#include "kacz_loop_clr_def.h"
#else
#include "kacz_loop_seq_def.h"
#endif
end subroutine SUB_NAME

#undef SUB_NAME
