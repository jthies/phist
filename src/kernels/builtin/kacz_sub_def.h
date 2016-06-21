subroutine SUB_NAME(nvec,nlocal, nhalo, ncols, nnz, &
        row_ptr, halo_ptr, col_idx, val, map, &
        shift_r,shift_i, b, ldb, &
        x_r,x_i, ldx, halo_r, halo_i,omega,&
        istart,iend,istep,i0,i1,j0,j1)

#ifdef PHIST_HAVE_OPENMP
  use :: omp_lib
#endif
  use :: map_module

  implicit none

  integer, intent(in) :: nvec, nlocal, nhalo, ncols, ldx, ldb
  integer(kind=8), intent(in) :: nnz
  real(kind=8), intent(in) :: shift_r(NVEC), shift_i(NVEC)
  integer(kind=8), intent(in) :: row_ptr(nlocal+1), halo_ptr(nlocal)
  integer, intent(in) :: col_idx(nnz)
  real(kind=8), intent(in) :: val(nnz)
  TYPE(Map_t), intent(in) :: map
  real(kind=8), intent(inout) :: x_r(NVEC,nlocal), x_i(NVEC,nlocal),b(NVEC,nlocal)
  real(kind=8), intent(inout) :: halo_r(NVEC,nhalo),halo_i(NVEC,nhalo)
  real(kind=8), intent(in) :: omega(NVEC)
  integer, intent(in) :: i0, i1, j0, j1, istart, iend, istep
  ! locals
  real(kind=8) :: tmp_r(NVEC), tmp_i(NVEC)
  integer :: i, ic, jc
  integer(kind=8) :: j
  integer istart_clr, iend_clr
  real(kind=8) :: row_norm(NVEC)
  real(kind=8) :: d
#ifdef NVEC
!dir$ assume_aligned row_ptr:64, halo_ptr:64, col_idx:64, val:64, x_r:64, halo_r:64, x_i:64, halo_i:8
#else
!dir$ assume_aligned row_ptr:64, halo_ptr:64, col_idx:64, val:64, x_r:8, halo_r:64, x_i:8, halo_i:8
#endif

! include file that contains the outer (i-) loop over the matrix rows

#ifdef KACZ_CLR
#include "kacz_loop_clr_def.h"
#else
#include "kacz_loop_seq_def.h"
#endif

    !! in local row i

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! compute (A-shift_j I)_i*x
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    d=0.0_8
    row_norm(1:NVEC)=0.0_8 ! compute row 2-norm on-the-fly
#ifdef KACZ_NO_SHIFT
    tmp_r(1:NVEC) = 0.0_8
#else
# ifndef KACZ_RC_VARIANT
    tmp_r(1:NVEC) = -shift_r(1:NVEC)*x_r(1:NVEC,i)
# else
    tmp_r(1:NVEC) = -shift_r(1:NVEC)*x_r(1:NVEC,i) &
                    +shift_i(1:NVEC)*x_i(1:NVEC,i)
    tmp_i(1:NVEC) = -shift_i(1:NVEC)*x_r(1:NVEC,i) &
                    -shift_r(1:NVEC)*x_i(1:NVEC,i)
# endif
#endif

    do j = row_ptr(i), halo_ptr(i)-1, 1
      tmp_r(1:NVEC) = tmp_r(1:NVEC) + val(j)*x_r(1:NVEC,col_idx(j))
#ifdef KACZ_RC_VARIANT
      tmp_i(1:NVEC) = tmp_i(1:NVEC) + val(j)*x_i(1:NVEC,col_idx(j))
#endif
      row_norm(1)=row_norm(1)+val(j)*val(j)
#ifndef KACZ_NO_SHIFT
      if (col_idx(j)==i) then
        d=val(j)
      end if
#endif
    end do
    do j = halo_ptr(i), row_ptr(i+1)-1, 1
      tmp_r(1:NVEC) = tmp_r(1:NVEC) + val(j)*halo_r(1:NVEC,col_idx(j))
#ifdef KACZ_RC_VARIANT
      tmp_i(1:NVEC) = tmp_i(1:NVEC) + val(j)*halo_i(1:NVEC,col_idx(j))
#endif
      row_norm(1)=row_norm(1)+val(j)*val(j)
    end do

#ifndef KACZ_BZERO
    tmp_r(1:NVEC)=tmp_r(1:NVEC)-b(1:NVEC,i)
#endif

#ifdef KACZ_NO_SHIFT
  row_norm(1:NVEC) = row_norm(1)
#else
  ! correct for real or complex shift
  row_norm(1:NVEC) = row_norm(1) - 2*d*shift_r(1:NVEC)+shift_r(1:NVEC)*shift_r(1:NVEC)
# ifdef KACZ_RC_VARIANT
  row_norm(1:NVEC) = row_norm(1:NVEC)+shift_i(1:NVEC)*shift_i(1:NVEC)
# endif
#endif

    ! Kaczmarz update of X. tmp is now (A-s_jI)_i x -b

    ! a) scaling factors
    WHERE(row_norm(1:NVEC)==0.0_8) row_norm=1.0_8
    tmp_r(1:NVEC)=tmp_r(1:NVEC)*omega(1:NVEC)/row_norm(1:NVEC)
#ifdef KACZ_RC_VARIANT
    tmp_i(1:NVEC)=tmp_i(1:NVEC)*omega(1:NVEC)/row_norm(1:NVEC)
#endif
    ! b) projection step
#ifndef KACZ_NO_SHIFT
# ifdef KACZ_RC_VARIANT
    ! note that in the complex case we need to take the conjugate somewhere (AA^Hx=b),
    ! we do it here in the projection, resulting in a sign change of shift_i
    x_r(1:NVEC,i)=x_r(1:NVEC,i) + (tmp_r(1:NVEC)*shift_r(1:NVEC)+tmp_i(1:NVEC)*shift_i(1:NVEC))
    x_i(1:NVEC,i)=x_i(1:NVEC,i) + (tmp_i(1:NVEC)*shift_r(1:NVEC)-tmp_r(1:NVEC)*shift_i(1:NVEC))
# else
    x_r(1:NVEC,i)=x_r(1:NVEC,i) + tmp_r(1:NVEC)*shift_r(1:NVEC)
# endif
#endif
    do j = row_ptr(i), halo_ptr(i)-1, 1
      x_r(1:NVEC,col_idx(j)) = x_r(1:NVEC,col_idx(j)) - &
                               tmp_r(1:NVEC)*val(j)
#ifdef KACZ_RC_VARIANT
      x_i(1:NVEC,col_idx(j)) = x_i(1:NVEC,col_idx(j)) - &
                               tmp_i(1:NVEC)*val(j)
#endif
    end do
    do j = halo_ptr(i), row_ptr(i+1)-1, 1
      halo_r(1:NVEC,col_idx(j)) = halo_r(1:NVEC,col_idx(j)) - &
                               tmp_r(1:NVEC)*val(j)
#ifdef KACZ_RC_VARIANT
      halo_i(1:NVEC,col_idx(j)) = halo_i(1:NVEC,col_idx(j)) - &
                               tmp_i(1:NVEC)*val(j)
#endif
    end do
  end do
#ifdef KACZ_CLR
  end do
#endif
end subroutine SUB_NAME

#undef SUB_NAME
#ifdef NVEC
#undef NVEC
#endif

