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

  integer, intent(in) :: nvec, nlocal, nhalo, ncols, ldx, ldb
  integer(kind=8), intent(in) :: nnz
  real(kind=8), intent(in) :: shift_r, shift_i
  integer(kind=8), intent(in) :: row_ptr(nlocal+1), halo_ptr(nlocal)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  TYPE(Map_t), intent(in) :: map
  real(kind=8), intent(inout) :: x_r(NVEC,nlocal), x_i(NVEC,nlocal),b(NVEC,nlocal)
  real(kind=8), intent(inout) :: halo_r(NVEC,nhalo),halo_i(NVEC,nhalo)
  real(kind=8), intent(in) :: nrms_ai2i(nlocal)
  real(kind=8), intent(in) :: omega
  integer, intent(in) :: i0, i1, j0, j1, istart, iend, istep
  ! locals
  real(kind=8) :: tmp_r(NVEC), tmp_i(NVEC)
  integer :: i, ic, jc
  integer(kind=8) :: j
  integer istart_clr, iend_clr
  real(kind=8) :: row_norm
! TODO - we don't actually check the memory alignment before calling these subroutines!
#if 0
#ifdef NVEC
!dir$ assume_aligned row_ptr:64, halo_ptr:64, col_idx:64, val:64, x_r:64, halo_r:64, x_i:64, halo_i:64, b:64
#else
!dir$ assume_aligned row_ptr:64, halo_ptr:64, col_idx:64, val:64, x_r:8, halo_r:64, x_i:8, halo_i:64, b:8
#endif
#endif
#ifdef KACZ_CLR
#include "kacz_loop_clr_def.h"
#else
#include "kacz_loop_seq_def.h"
#endif
end subroutine SUB_NAME

#undef SUB_NAME
