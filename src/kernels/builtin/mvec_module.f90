#include "fdebug.h"
!> \file mvec_module.f90
!! Defines mvec_module, the phist builtin implementation of phist_Dmvec_*
!! \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
!!

#include "phist_config.h"
#include "phist_kernel_flags.h"
#include "phist_defs.h"
!> Implementations of phist_Dmvec_* for phist builtin kernels which uses row-major mvecs
!! 
!! Actual work is delegated to fast, parallel implementations tuned for LARGE
!! (e.g. much larger than CPU cache) block vectors
module mvec_module
  use env_module,   only: allocate_aligned, deallocate_aligned
  use map_module,   only: Map_t, map_compatible_map
  use sdmat_module, only: SDMat_t
  implicit none
  private

!>@todo duplicated code from crsmat_module.f90
#ifndef PHIST_HAVE_GHOST
#define G_LIDX_T C_INT32_T
#define G_GIDX_T C_INT64_T
#else
#ifdef GHOST_IDX64_LOCAL
#define G_LIDX_T C_INT64_T
#else
#define G_LIDX_T C_INT32_T
#endif
#ifdef GHOST_IDX64_GLOBAL
#define G_GIDX_T C_INT64_T
#else
#define G_GIDX_T C_INT32_T
#endif
#endif

#if PHIST_HIGH_PRECISION_KERNELS_FORCE
#define FORCE_HIGH_PRECISION .true.
#else
#define FORCE_HIGH_PRECISION .false.
#endif

  public :: MVec_t
  !public :: phist_Dmvec_create
  !public :: phist_Dmvec_delete
  !public :: phist_Dmvec_extract_view
  !public :: phist_Dmvec_my_length
  !public :: phist_Dmvec_get_map
  !public :: phist_Dmvec_num_vectors
  !public :: phist_Dmvec_view_block
  !public :: phist_Dmvec_get_block
  !public :: phist_Dmvec_set_block
  !public :: phist_Dmvec_gather_mvecs
  !public :: phist_Dmvec_scatter_mvecs
  public :: mvec_gather_mvecs
  public :: mvec_scatter_mvecs
  !public :: phist_Dmvec_put_value
  !public :: phist_Dmvec_random
  public :: mvec_random
  !public :: phist_Dmvec_print
  !public :: phist_Dmvec_norm2
  public :: mvec_norm2
  !public :: phist_Dmvec_scale
  !public :: phist_Dmvec_vscale
  public :: mvec_scale
  public :: mvec_vscale
  !public :: phist_Dmvec_add_mvec
  !public :: phist_Dmvec_vadd_mvec
  public :: mvec_add_mvec
  public :: mvec_vadd_mvec
  !public :: phist_Dmvec_dot_mvec
  public :: mvec_dot_mvec
  !public :: phist_Dmvec_times_sdMat
  !public :: phist_Dmvec_times_sdMat_inplace
  public :: mvec_times_sdmat
  public :: mvec_times_sdmat_inplace
  !public :: phist_DmvecT_times_mvec
  !public :: phist_DmvecT_times_mvec_times_sdMat_inplace
  public :: mvecT_times_mvec
  public :: mvecT_times_mvec_times_sdMat_inplace
  !public :: phist_Dmvec_to_mvec
  !public :: phist_Dmvec_QR
  public :: mvec_QR


#ifndef PHIST_HAVE_GHOST
#define G_LIDX_T C_INT32_T
#define G_GIDX_T C_INT64_T
#else
#ifdef GHOST_IDX64_LOCAL
#define G_LIDX_T C_INT64_T
#else
#define G_LIDX_T C_INT32_T
#endif
#ifdef GHOST_IDX64_GLOBAL
#define G_GIDX_T C_INT64_T
#else
#define G_GIDX_T C_INT32_T
#endif
#endif


  !==================================================================================
  !> mvec with row-wise layout
  type MVec_t
    !--------------------------------------------------------------------------------
    integer     :: jmin, jmax
    integer :: paddedN
    type(Map_t) :: map
    real(kind=8), contiguous, pointer :: val(:,:) => null()
    logical     :: is_view
    !--------------------------------------------------------------------------------
  end type MVec_t


  !==================================================================================
  ! interfaces to low-level kernel functions
  interface
    !void drandom_1(int nrows, double *restrict y, int64_t pre_skip, int64_t post_skip);
    subroutine drandom_1(nrows, y, pre_skip, post_skip) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE, C_INT64_T
      integer(kind=C_INT), value :: nrows
      real(kind=C_DOUBLE), intent(out) :: y
      integer(kind=C_INT64_T), value :: pre_skip, post_skip
    end subroutine
    !void drandom_general(int nvec, int nrows, double *restrict v, int ldv, int64_t pre_skip, int64_t post_skip)
    subroutine drandom_general(nvec, nrows, y, ldv, pre_skip, post_skip) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE, C_INT64_T
      integer(kind=C_INT), value :: nvec, nrows, ldv
      real(kind=C_DOUBLE), intent(out) :: y
      integer(kind=C_INT64_T), value :: pre_skip, post_skip
    end subroutine
    !void ddot_self_prec_1(int nrows, const double *restrict x, double *restrict res, double *restrict resC)
    subroutine ddot_self_prec_1(nrows, x, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: nrows
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(out) :: res(*), resC(*)
    end subroutine
    !void ddot_fused_scale_self_prec_1(int nrows, double *restrict x, double scal, double scalC, double *restrict res, double *restrict resC)
    subroutine ddot_fused_scale_self_prec_1(nrows, x, scal, scalC, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: nrows
      real(kind=C_DOUBLE), intent(inout) :: x
      real(kind=C_DOUBLE), value :: scal, scalC
      real(kind=C_DOUBLE), intent(out) :: res(*), resC(*)
    end subroutine
    !void ddot_self_prec_2(int nrows, const double *restrict x, double *restrict res, double *restrict resC)
    subroutine ddot_self_prec_2(nrows, x, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: nrows
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(out) :: res(*), resC(*)
    end subroutine
    !void ddot_self_prec_4(int nrows, const double *restrict x, double *restrict res, double *restrict resC)
    subroutine ddot_self_prec_4(nrows, x, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: nrows
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(out) :: res(*), resC(*)
    end subroutine
    !void ddot_prec_1(int nrows, const double *restrict x, const double *restrict y, double *restrict res, double *restrict resC)
    subroutine ddot_prec_1(nrows, x, y, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: nrows
      real(kind=C_DOUBLE), intent(in) :: x, y
      real(kind=C_DOUBLE), intent(out) :: res(*), resC(*)
    end subroutine
    !void ddot_fused_scale_prec_1(int nrows, const double *restrict x, const double *restrict y, double scal, double scalC, double *restrict res, double *restrict resC)
    subroutine ddot_fused_scale_prec_1(nrows, x, y, scal, scalC, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: nrows
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), value :: scal, scalC
      real(kind=C_DOUBLE), intent(out) :: res(*), resC(*)
    end subroutine
    !void ddot_prec_2(int nrows, const double *restrict x, const double *restrict y, double *restrict res, double *restrict resC)
    subroutine ddot_prec_2(nrows, x, y, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: nrows
      real(kind=C_DOUBLE), intent(in) :: x, y
      real(kind=C_DOUBLE), intent(out) :: res(*), resC(*)
    end subroutine
    !void ddot_prec_4(int nrows, const double *restrict x, const double *restrict y, double *restrict res, double *restrict resC)
    subroutine ddot_prec_4(nrows, x, y, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: nrows
      real(kind=C_DOUBLE), intent(in) :: x, y
      real(kind=C_DOUBLE), intent(out) :: res(*), resC(*)
    end subroutine
    !void daxpby_prec(int n, double alpha, const double *restrict a, const double *restrict aC,
    !                        double beta,        double *restrict b,       double *restrict bC)
    subroutine daxpby_prec(n,alpha,a,aC,beta,b,bC) bind(C)
      use, intrinsic :: iso_c_binding
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), value :: alpha, beta
      real(kind=C_DOUBLE), intent(in) :: a(*), aC(*)
      real(kind=C_DOUBLE), intent(inout) :: b(*), bC(*)
    end subroutine
    !void prec_reduction_1(int n, const double *restrict s_, const double *restrict c_, double *restrict r, double *restrict rC)
    subroutine prec_reduction_1(n, s, c, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), intent(in) :: s(*), c(*)
      real(kind=C_DOUBLE), intent(out) :: res(*), resC(*)
    end subroutine
    !void prec_reduction_2(int n, const double *restrict s_, const double *restrict c_, double *restrict r, double *restrict rC)
    subroutine prec_reduction_2(n, s, c, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), intent(in) :: s(*), c(*)
      real(kind=C_DOUBLE), intent(out) :: res(*), resC(*)
    end subroutine
    !void prec_reduction_2k(int n, int k, const double *restrict s_, const double *restrict c_, double *restrict r, double *restrict rC)
    subroutine prec_reduction_2k(n, k, s, c, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k
      real(kind=C_DOUBLE), intent(in) :: s(*), c(*)
      real(kind=C_DOUBLE), intent(out) :: res(*), resC(*)
    end subroutine
    !void prec_reduction_4(int n, const double *restrict s_, const double *restrict c_, double *restrict r, double *restrict rC)
    subroutine prec_reduction_4(n, s, c, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), intent(in) :: s(*), c(*)
      real(kind=C_DOUBLE), intent(out) :: res(*), resC(*)
    end subroutine
    !void prec_reduction_4k(int n, int k, const double *restrict s_, const double *restrict c_, double *restrict r, double *restrict rC)
    subroutine prec_reduction_4k(n, k, s, c, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k
      real(kind=C_DOUBLE), intent(in) :: s(*), c(*)
      real(kind=C_DOUBLE), intent(out) :: res(*), resC(*)
    end subroutine
    !void prec_reduction_k(int n, int k, const double *restrict s_, const double *restrict c_, double *restrict r, double *restrict rC)
    subroutine prec_reduction_k(n, k, s, c, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k
      real(kind=C_DOUBLE), intent(in) :: s(*), c(*)
      real(kind=C_DOUBLE), intent(out) :: res(*), resC(*)
    end subroutine
    !void dgemm_sC_self_prec_4(int nrows, const double *restrict x, double *restrict res, double *restrict resC)
    subroutine dgemm_sC_self_prec_4(n, x, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(out) :: res, resC
    end subroutine
    !void dgemm_sC_self_prec_2(int nrows, const double *restrict x, double *restrict res, double *restrict resC)
    subroutine dgemm_sC_self_prec_2(n, x, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(out) :: res, resC
    end subroutine
    !void dgemm_sC_prec_4_4(int nrows, const double *restrict x, const double *restrict y, double *restrict res, double *restrict resC)
    subroutine dgemm_sC_prec_4_4(n, x, y, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), intent(in) :: x, y
      real(kind=C_DOUBLE), intent(out) :: res(*), resC(*)
    end subroutine
    !void dgemm_sC_prec_2_2(int nrows, const double *restrict x, const double *restrict y, double *restrict res, double *restrict resC)
    subroutine dgemm_sC_prec_2_2(n, x, y, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), intent(in) :: x, y
      real(kind=C_DOUBLE), intent(out) :: res(*), resC(*)
    end subroutine
    !void dgemm_sC_prec_2_1(int nrows, const double *restrict x, const double *restrict y, double *restrict res, double *restrict resC)
    subroutine dgemm_sC_prec_2_1(n, x, y, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), intent(in) :: x, y
      real(kind=C_DOUBLE), intent(out) :: res(*), resC(*)
    end subroutine
    !void dgemm_sC_prec_4_2(int nrows, const double *restrict x, const double *restrict y, double *restrict res, double *restrict resC)
    subroutine dgemm_sC_prec_4_2(n, x, y, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), intent(in) :: x, y
      real(kind=C_DOUBLE), intent(out) :: res(*), resC(*)
    end subroutine
    !void dgemm_sC_prec_4_1(int nrows, const double *restrict x, const double *restrict y, double *restrict res, double *restrict resC)
    subroutine dgemm_sC_prec_4_1(n, x, y, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), intent(in) :: x, y
      real(kind=C_DOUBLE), intent(out) :: res(*), resC(*)
    end subroutine
    !void dgemm_sC_prec_4_k(int nrows, int k, const double *restrict x, const double *restrict y, double *restrict res, double *restrict resC)
    subroutine dgemm_sC_prec_4_k(n, k, x, y, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k
      real(kind=C_DOUBLE), intent(in) :: x, y
      real(kind=C_DOUBLE), intent(out) :: res(*), resC(*)
    end subroutine
    !void dgemm_sC_prec_2_k(int nrows, int k, const double *restrict x, const double *restrict y, double *restrict res, double *restrict resC)
    subroutine dgemm_sC_prec_2_k(n, k, x, y, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k
      real(kind=C_DOUBLE), intent(in) :: x, y
      real(kind=C_DOUBLE), intent(out) :: res(*), resC(*)
    end subroutine
    !void dgemm_sC_prec_1_k(int nrows, int k, const double *restrict x, const double *restrict y, double *restrict res, double *restrict resC)
    subroutine dgemm_sC_prec_1_k(n, k, x, y, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k
      real(kind=C_DOUBLE), intent(in) :: x, y
      real(kind=C_DOUBLE), intent(out) :: res(*), resC(*)
    end subroutine

    !void dgemm_fused_sCD_self_prec_4(int nrows, double *restrict x, const double *restrict scal, const double *restrict scalC, double *restrict res, double *restrict resC)
    subroutine dgemm_fused_sCD_self_prec_4(n, x, scal, scalC, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), intent(inout) :: x
      real(kind=C_DOUBLE), intent(in) :: scal(*), scalC(*)
      real(kind=C_DOUBLE), intent(out) :: res, resC
    end subroutine
    !void dgemm_fused_sCD_self_prec_2(int nrows, double *restrict x, const double *restrict scal, const double *restrict scalC, double *restrict res, double *restrict resC)
    subroutine dgemm_fused_sCD_self_prec_2(n, x, scal, scalC, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), intent(inout) :: x
      real(kind=C_DOUBLE), intent(in) :: scal(*), scalC(*)
      real(kind=C_DOUBLE), intent(out) :: res, resC
    end subroutine
    !void dgemm_fused_sCD_prec_4_4(int nrows, const double *restrict x, double *restrict y, const double *restrict scal, const double *restrict scalC, double *restrict res, double *restrict resC)
    subroutine dgemm_fused_sCD_prec_4_4(n, x, y, scal, scalC, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(in) :: scal(*), scalC(*)
      real(kind=C_DOUBLE), intent(out) :: res(*), resC(*)
    end subroutine
    !void dgemm_fused_sCD_prec_2_2(int nrows, const double *restrict x, double *restrict y, const double *restrict scal, const double *restrict scalC, double *restrict res, double *restrict resC)
    subroutine dgemm_fused_sCD_prec_2_2(n, x, y, scal, scalC, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(in) :: scal(*), scalC(*)
      real(kind=C_DOUBLE), intent(out) :: res(*), resC(*)
    end subroutine
    !void dgemm_fused_sCD_prec_1_2(int nrows, const double *restrict x, double *restrict y, const double *restrict scal, const double *restrict scalC, double *restrict res, double *restrict resC)
    subroutine dgemm_fused_sCD_prec_1_2(n, x, y, scal, scalC, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(in) :: scal(*), scalC(*)
      real(kind=C_DOUBLE), intent(out) :: res(*), resC(*)
    end subroutine
    !void dgemm_fused_sCD_prec_2_1(int nrows, const double *restrict x, double *restrict y, const double *restrict scal, const double *restrict scalC, double *restrict res, double *restrict resC)
    subroutine dgemm_fused_sCD_prec_2_1(n, x, y, scal, scalC, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(in) :: scal(*), scalC(*)
      real(kind=C_DOUBLE), intent(out) :: res(*), resC(*)
    end subroutine
    !void dgemm_fused_sCD_prec_2_4(int nrows, const double *restrict x, double *restrict y, const double *restrict scal, const double *restrict scalC, double *restrict res, double *restrict resC)
    subroutine dgemm_fused_sCD_prec_2_4(n, x, y, scal, scalC, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(in) :: scal(*), scalC(*)
      real(kind=C_DOUBLE), intent(out) :: res(*), resC(*)
    end subroutine
    !void dgemm_fused_sCD_prec_4_2(int nrows, const double *restrict x, double *restrict y, const double *restrict scal, const double *restrict scalC, double *restrict res, double *restrict resC)
    subroutine dgemm_fused_sCD_prec_4_2(n, x, y, scal, scalC, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(in) :: scal(*), scalC(*)
      real(kind=C_DOUBLE), intent(out) :: res(*), resC(*)
    end subroutine
    !void dgemm_fused_sCD_prec_1_4(int nrows, const double *restrict x, double *restrict y, const double *restrict scal, const double *restrict scalC, double *restrict res, double *restrict resC)
    subroutine dgemm_fused_sCD_prec_1_4(n, x, y, scal, scalC, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(in) :: scal(*), scalC(*)
      real(kind=C_DOUBLE), intent(out) :: res(*), resC(*)
    end subroutine
    !void dgemm_fused_sCD_prec_4_1(int nrows, const double *restrict x, double *restrict y, const double *restrict scal, const double *restrict scalC, double *restrict res, double *restrict resC)
    subroutine dgemm_fused_sCD_prec_4_1(n, x, y, scal, scalC, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(in) :: scal(*), scalC(*)
      real(kind=C_DOUBLE), intent(out) :: res(*), resC(*)
    end subroutine
    !void dgemm_fused_sCD_prec_k_4(int nrows, int k, const double *restrict x, double *restrict y, const double *restrict scal, const double *restrict scalC, double *restrict res, double *restrict resC)
    subroutine dgemm_fused_sCD_prec_k_4(n, k, x, y, scal, scalC, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(in) :: scal(*), scalC(*)
      real(kind=C_DOUBLE), intent(out) :: res(*), resC(*)
    end subroutine
    !void dgemm_fused_sCD_prec_k_2(int nrows, int k, const double *restrict x, double *restrict y, const double *restrict scal, const double *restrict scalC, double *restrict res, double *restrict resC)
    subroutine dgemm_fused_sCD_prec_k_2(n, k, x, y, scal, scalC, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(in) :: scal(*), scalC(*)
      real(kind=C_DOUBLE), intent(out) :: res(*), resC(*)
    end subroutine
    !void dgemm_fused_sCD_prec_k_1(int nrows, int k, const double *restrict x, double *restrict y, const double *restrict scal, const double *restrict scalC, double *restrict res, double *restrict resC)
    subroutine dgemm_fused_sCD_prec_k_1(n, k, x, y, scal, scalC, res, resC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(in) :: scal, scalC
      real(kind=C_DOUBLE), intent(out) :: res(*), resC(*)
    end subroutine
    !void dgemm_sb_inplace_prec_4(int nrows, double *restrict x, const double *restrict r, const double *restrict rC)
    subroutine dgemm_sb_inplace_prec_4(n,x,r,rC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), intent(inout) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
    end subroutine
    !void dgemm_sb_inplace_prec_2(int nrows, double *restrict x, const double *restrict r, const double *restrict rC)
    subroutine dgemm_sb_inplace_prec_2(n,x,r,rC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), intent(inout) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
    end subroutine
    !void dgemm_sb_inplace_prec_1(int nrows, double *restrict x, const double *restrict r, const double *restrict rC)
    subroutine dgemm_sb_inplace_prec_1(n,x,r,rC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), intent(inout) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
    end subroutine
    !void dgemm_sb_prec_k_4(int nrows, int k, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y)
    subroutine dgemm_sb_prec_k_4(n,k,alpha,x,r,rC,beta,y) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k
      real(kind=C_DOUBLE), value :: alpha, beta
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
    end subroutine
    !void dgemm_sb_prec_k_2(int nrows, int k, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y)
    subroutine dgemm_sb_prec_k_2(n,k,alpha,x,r,rC,beta,y) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k
      real(kind=C_DOUBLE), value :: alpha, beta
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
    end subroutine
    !void dgemm_sb_prec_k_1(int nrows, int k, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y)
    subroutine dgemm_sb_prec_k_1(n,k,alpha,x,r,rC,beta,y) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k
      real(kind=C_DOUBLE), value :: alpha, beta
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
    end subroutine
    !void dgemm_sb_prec_k_4_nt(int nrows, int k, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y)
    subroutine dgemm_sb_prec_k_4_nt(n,k,alpha,x,r,rC,y) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k
      real(kind=C_DOUBLE), value :: alpha
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(out):: y
    end subroutine
    !void dgemm_sb_prec_k_2_nt(int nrows, int k, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y)
    subroutine dgemm_sb_prec_k_2_nt(n,k,alpha,x,r,rC,y) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k
      real(kind=C_DOUBLE), value :: alpha
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(out):: y
    end subroutine
    !void dgemm_sb_prec_k_1_nt(int nrows, int k, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y)
    subroutine dgemm_sb_prec_k_1_nt(n,k,alpha,x,r,rC,y) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k
      real(kind=C_DOUBLE), value :: alpha
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(out):: y
    end subroutine
    !void dgemm_sb_prec_k_strided_4(int nrows, int k, double alpha, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double beta, double *restrict y)
    subroutine dgemm_sb_prec_k_strided_4(n,k,alpha,x,ldx,r,rC,beta,y) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k, ldx
      real(kind=C_DOUBLE), value :: alpha, beta
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
    end subroutine
    !void dgemm_sb_prec_k_strided_2(int nrows, int k, double alpha, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double beta, double *restrict y)
    subroutine dgemm_sb_prec_k_strided_2(n,k,alpha,x,ldx,r,rC,beta,y) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k, ldx
      real(kind=C_DOUBLE), value :: alpha, beta
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
    end subroutine
    !void dgemm_sb_prec_k_strided_1(int nrows, int k, double alpha, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double beta, double *restrict y)
    subroutine dgemm_sb_prec_k_strided_1(n,k,alpha,x,ldx,r,rC,beta,y) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k, ldx
      real(kind=C_DOUBLE), value :: alpha, beta
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
    end subroutine
    !void dgemm_sb_prec_k_strided_4_nt(int nrows, int k, double alpha, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double *restrict y)
    subroutine dgemm_sb_prec_k_strided_4_nt(n,k,alpha,x,ldx,r,rC,y) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k, ldx
      real(kind=C_DOUBLE), value :: alpha
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(out):: y
    end subroutine
    !void dgemm_sb_prec_k_strided_2_nt(int nrows, int k, double alpha, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double *restrict y)
    subroutine dgemm_sb_prec_k_strided_2_nt(n,k,alpha,x,ldx,r,rC,y) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k, ldx
      real(kind=C_DOUBLE), value :: alpha
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(out):: y
    end subroutine
    !void dgemm_sb_prec_k_strided_1_nt(int nrows, int k, double alpha, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double *restrict y)
    subroutine dgemm_sb_prec_k_strided_1_nt(n,k,alpha,x,ldx,r,rC,y) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k, ldx
      real(kind=C_DOUBLE), value :: alpha
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(out):: y
    end subroutine
    !void dgemm_sb_augmented_prec_strided_k_4(int nrows, int k, double alpha, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
    subroutine dgemm_sb_augmented_prec_strided_k_4(n,k,alpha,x,ldx,r,rC,beta,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k, ldx
      real(kind=C_DOUBLE), value :: alpha, beta
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(out) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_augmented_prec_k_4(int nrows, int k, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
    subroutine dgemm_sb_augmented_prec_k_4(n,k,alpha,x,r,rC,beta,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k
      real(kind=C_DOUBLE), value :: alpha, beta
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(out) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_augmented_prec_4_4(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
    subroutine dgemm_sb_augmented_prec_4_4(n,alpha,x,r,rC,beta,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), value :: alpha, beta
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(out) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_augmented_prec_2_4(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
    subroutine dgemm_sb_augmented_prec_2_4(n,alpha,x,r,rC,beta,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), value :: alpha, beta
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(out) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_augmented_prec_1_4(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
    subroutine dgemm_sb_augmented_prec_1_4(n,alpha,x,r,rC,beta,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), value :: alpha, beta
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(out) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_augmented_prec_strided_k_2(int nrows, int k, double alpha, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
    subroutine dgemm_sb_augmented_prec_strided_k_2(n,k,alpha,x,ldx,r,rC,beta,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k, ldx
      real(kind=C_DOUBLE), value :: alpha, beta
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(out) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_augmented_prec_k_2(int nrows, int k, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
    subroutine dgemm_sb_augmented_prec_k_2(n,k,alpha,x,r,rC,beta,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k
      real(kind=C_DOUBLE), value :: alpha, beta
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(out) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_augmented_prec_4_2(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
    subroutine dgemm_sb_augmented_prec_4_2(n,alpha,x,r,rC,beta,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), value :: alpha, beta
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(out) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_augmented_prec_2_2(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
    subroutine dgemm_sb_augmented_prec_2_2(n,alpha,x,r,rC,beta,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), value :: alpha, beta
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(out) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_augmented_prec_1_2(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
    subroutine dgemm_sb_augmented_prec_1_2(n,alpha,x,r,rC,beta,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), value :: alpha, beta
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(out) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_augmented_prec_strided_k_1(int nrows, int k, double alpha, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
    subroutine dgemm_sb_augmented_prec_strided_k_1(n,k,alpha,x,ldx,r,rC,beta,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k, ldx
      real(kind=C_DOUBLE), value :: alpha, beta
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(out) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_augmented_prec_k_1(int nrows, int k, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
    subroutine dgemm_sb_augmented_prec_k_1(n,k,alpha,x,r,rC,beta,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k
      real(kind=C_DOUBLE), value :: alpha, beta
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(out) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_augmented_prec_4_1(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
    subroutine dgemm_sb_augmented_prec_4_1(n,alpha,x,r,rC,beta,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), value :: alpha, beta
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(out) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_augmented_prec_2_1(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
    subroutine dgemm_sb_augmented_prec_2_1(n,alpha,x,r,rC,beta,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), value :: alpha, beta
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(out) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_augmented_prec_1_1(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
    subroutine dgemm_sb_augmented_prec_1_1(n,alpha,x,r,rC,beta,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), value :: alpha, beta
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(out) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_augmented_prec_strided_k_4_nt(int nrows, int k, double alpha, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
    subroutine dgemm_sb_augmented_prec_strided_k_4_nt(n,k,alpha,x,ldx,r,rC,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k, ldx
      real(kind=C_DOUBLE), value :: alpha
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(out):: y
      real(kind=C_DOUBLE), intent(out) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_augmented_prec_k_4_nt(int nrows, int k, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
    subroutine dgemm_sb_augmented_prec_k_4_nt(n,k,alpha,x,r,rC,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k
      real(kind=C_DOUBLE), value :: alpha
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(out):: y
      real(kind=C_DOUBLE), intent(out) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_augmented_prec_4_4_nt(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
    subroutine dgemm_sb_augmented_prec_4_4_nt(n,alpha,x,r,rC,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), value :: alpha
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(out):: y
      real(kind=C_DOUBLE), intent(out) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_augmented_prec_2_4_nt(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
    subroutine dgemm_sb_augmented_prec_2_4_nt(n,alpha,x,r,rC,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), value :: alpha
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(out):: y
      real(kind=C_DOUBLE), intent(out) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_augmented_prec_1_4_nt(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
    subroutine dgemm_sb_augmented_prec_1_4_nt(n,alpha,x,r,rC,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), value :: alpha
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(out):: y
      real(kind=C_DOUBLE), intent(out) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_augmented_prec_strided_k_2_nt(int nrows, int k, double alpha, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
    subroutine dgemm_sb_augmented_prec_strided_k_2_nt(n,k,alpha,x,ldx,r,rC,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k, ldx
      real(kind=C_DOUBLE), value :: alpha
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(out):: y
      real(kind=C_DOUBLE), intent(out) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_augmented_prec_k_2_nt(int nrows, int k, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
    subroutine dgemm_sb_augmented_prec_k_2_nt(n,k,alpha,x,r,rC,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k
      real(kind=C_DOUBLE), value :: alpha
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(out):: y
      real(kind=C_DOUBLE), intent(out) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_augmented_prec_4_2_nt(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
    subroutine dgemm_sb_augmented_prec_4_2_nt(n,alpha,x,r,rC,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), value :: alpha
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(out):: y
      real(kind=C_DOUBLE), intent(out) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_augmented_prec_2_2_nt(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
    subroutine dgemm_sb_augmented_prec_2_2_nt(n,alpha,x,r,rC,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), value :: alpha
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(out):: y
      real(kind=C_DOUBLE), intent(out) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_augmented_prec_1_2_nt(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
    subroutine dgemm_sb_augmented_prec_1_2_nt(n,alpha,x,r,rC,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), value :: alpha
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(out):: y
      real(kind=C_DOUBLE), intent(out) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_augmented_prec_strided_k_1_nt(int nrows, int k, double alpha, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
    subroutine dgemm_sb_augmented_prec_strided_k_1_nt(n,k,alpha,x,ldx,r,rC,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k, ldx
      real(kind=C_DOUBLE), value :: alpha
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(out):: y
      real(kind=C_DOUBLE), intent(out) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_augmented_prec_k_1_nt(int nrows, int k, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
    subroutine dgemm_sb_augmented_prec_k_1_nt(n,k,alpha,x,r,rC,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k
      real(kind=C_DOUBLE), value :: alpha
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(out):: y
      real(kind=C_DOUBLE), intent(out) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_augmented_prec_4_1_nt(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
    subroutine dgemm_sb_augmented_prec_4_1_nt(n,alpha,x,r,rC,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), value :: alpha
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(out):: y
      real(kind=C_DOUBLE), intent(out) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_augmented_prec_2_1_nt(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
    subroutine dgemm_sb_augmented_prec_2_1_nt(n,alpha,x,r,rC,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), value :: alpha
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(out):: y
      real(kind=C_DOUBLE), intent(out) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_augmented_prec_1_1_nt(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
    subroutine dgemm_sb_augmented_prec_1_1_nt(n,alpha,x,r,rC,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), value :: alpha
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(out):: y
      real(kind=C_DOUBLE), intent(out) :: d(*), dC(*)
    end subroutine

    !void dgemm_sb_add_sd_prec_strided_k_4(int nrows, int k, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double *restrict y, const double *restrict d, const double *restrict dC)
    subroutine dgemm_sb_add_sd_prec_strided_k_4(n,k,x,ldx,r,rC,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k, ldx
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(in) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_add_sd_prec_k_4(int nrows, int k, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, const double *restrict d, const double *restrict dC)
    subroutine dgemm_sb_add_sd_prec_k_4(n,k,x,r,rC,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(in) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_add_sd_prec_4_4(int nrows, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, const double *restrict d, const double *restrict dC)
    subroutine dgemm_sb_add_sd_prec_4_4(n,x,r,rC,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(in) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_add_sd_prec_2_4(int nrows, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, const double *restrict d, const double *restrict dC)
    subroutine dgemm_sb_add_sd_prec_2_4(n,x,r,rC,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(in) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_add_sd_prec_1_4(int nrows, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, const double *restrict d, const double *restrict dC)
    subroutine dgemm_sb_add_sd_prec_1_4(n,x,r,rC,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(in) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_add_sd_prec_strided_k_2(int nrows, int k, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double *restrict y, const double *restrict d, const double *restrict dC)
    subroutine dgemm_sb_add_sd_prec_strided_k_2(n,k,x,ldx,r,rC,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k, ldx
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(in) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_add_sd_prec_k_2(int nrows, int k, double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, const double *restrict d, const double *restrict dC)
    subroutine dgemm_sb_add_sd_prec_k_2(n,k,x,r,rC,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(in) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_add_sd_prec_4_2(int nrows, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, const double *restrict d, const double *restrict dC)
    subroutine dgemm_sb_add_sd_prec_4_2(n,x,r,rC,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(in) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_add_sd_prec_2_2(int nrows, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, const double *restrict d, const double *restrict dC)
    subroutine dgemm_sb_add_sd_prec_2_2(n,x,r,rC,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(in) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_add_sd_prec_1_2(int nrows, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, const double *restrict d, const double *restrict dC)
    subroutine dgemm_sb_add_sd_prec_1_2(n,x,r,rC,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(in) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_add_sd_prec_strided_k_1(int nrows, int k, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double *restrict y, const double *restrict d, const double *restrict dC)
    subroutine dgemm_sb_add_sd_prec_strided_k_1(n,k,x,ldx,r,rC,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k, ldx
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(in) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_add_sd_prec_k_1(int nrows, int k, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, const double *restrict d, const double *restrict dC)
    subroutine dgemm_sb_add_sd_prec_k_1(n,k,x,r,rC,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n, k
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(in) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_add_sd_prec_4_1(int nrows, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, const double *restrict d, const double *restrict dC)
    subroutine dgemm_sb_add_sd_prec_4_1(n,x,r,rC,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(in) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_add_sd_prec_2_1(int nrows, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, const double *restrict d, const double *restrict dC)
    subroutine dgemm_sb_add_sd_prec_2_1(n,x,r,rC,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(in) :: d(*), dC(*)
    end subroutine
    !void dgemm_sb_add_sd_prec_1_1(int nrows, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, const double *restrict d, const double *restrict dC)
    subroutine dgemm_sb_add_sd_prec_1_1(n,x,r,rC,y,d,dC) bind(C)
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), intent(in) :: x
      real(kind=C_DOUBLE), intent(in) :: r(*),rC(*)
      real(kind=C_DOUBLE), intent(inout) :: y
      real(kind=C_DOUBLE), intent(in) :: d(*), dC(*)
    end subroutine
  end interface

  !> interface of function-ptr for mvec_put_func
  abstract interface
    function mvecElemFunc(row, col, val_ptr,last_arg) bind(C) result(ierr)
      use, intrinsic :: iso_c_binding
      integer(G_GIDX_T), value :: row
      integer(G_LIDX_T), value :: col
      TYPE(C_PTR), value :: val_ptr
      TYPE(C_PTR), value :: last_arg
      integer(C_INT) :: ierr
    end function mvecElemFunc
  end interface

contains

  !==================================================================================
  !> copy data from a list of other mvecs into one mvec
  subroutine mvec_gather_mvecs(mvec, block_list)
    !--------------------------------------------------------------------------------
    type(MVec_t), intent(inout) :: mvec
    type(MVec_t), intent(in)    :: block_list(1:)
    !--------------------------------------------------------------------------------
    integer :: nvec, nblocks, nrows, ldm
    logical :: single_vector_gather
    integer :: i, j, nvec_i, off, jmin_i, jmax_i
    logical :: strided
    !--------------------------------------------------------------------------------

    ! determine data layout
    nvec = mvec%jmax-mvec%jmin+1
    nblocks = size(block_list)
    nrows = mvec%map%nlocal(mvec%map%me)
    ldm = size(mvec%val,1)
    if( .not. mvec%is_view .or. &
      & ( mvec%jmin .eq. lbound(mvec%val,1) .and. &
      &   mvec%jmax .eq. ubound(mvec%val,1)       ) ) then
      strided = .false.
    else
      strided = .true.
    end if


    ! fast versions for fixed nvec and gathering single vectors
    single_vector_gather = .true.
    if( nvec .ne. nblocks ) single_vector_gather = .false.
    do i = 1, nblocks, 1
      if( block_list(i)%jmin .ne. block_list(i)%jmax ) then
        single_vector_gather = .false.
      end if
    end do

    if( single_vector_gather .and. .not. strided ) then
      if( nvec .eq. 1 ) then
        call dgather_1(nrows, mvec%val(mvec%jmin,1), &
          &            block_list(1)%val(block_list(1)%jmin,1), size(block_list(1)%val,1) )
        return
      else if( nvec .eq. 2 ) then
        call dgather_2(nrows, mvec%val(mvec%jmin,1), &
          &            block_list(1)%val(block_list(1)%jmin,1), size(block_list(1)%val,1), &
          &            block_list(2)%val(block_list(2)%jmin,1), size(block_list(2)%val,1) )
        return
      else if( nvec .eq. 4 ) then
        call dgather_4(nrows, mvec%val(mvec%jmin,1), &
          &            block_list(1)%val(block_list(1)%jmin,1), size(block_list(1)%val,1), &
          &            block_list(2)%val(block_list(2)%jmin,1), size(block_list(2)%val,1), &
          &            block_list(3)%val(block_list(3)%jmin,1), size(block_list(3)%val,1), &
          &            block_list(4)%val(block_list(4)%jmin,1), size(block_list(4)%val,1) )
        return
      end if
    end if

!$omp parallel do private(off) schedule(static)
    do j = 1, nrows, 1
      off = mvec%jmin
      do i = 1, nblocks, 1
        jmin_i = block_list(i)%jmin
        jmax_i = block_list(i)%jmax
        nvec_i = jmax_i-jmin_i+1
        mvec%val(off:off+nvec_i-1,j) = block_list(i)%val(jmin_i:jmax_i,j)
        off = off +nvec_i
      end do
    end do

    !--------------------------------------------------------------------------------
  end subroutine mvec_gather_mvecs


  ! currently unused
  !==================================================================================
  !> copy data from a list of other mvecs into one mvec
  subroutine mvec_scatter_mvecs(mvec, block_list)
    !--------------------------------------------------------------------------------
    type(MVec_t), intent(inout) :: mvec
    type(MVec_t), intent(in)    :: block_list(1:)
    !--------------------------------------------------------------------------------
    integer :: nvec, nblocks, nrows, ldm
    logical :: single_vector_scatter
    integer :: i, j, nvec_i, off, jmin_i, jmax_i
    logical :: strided
    !--------------------------------------------------------------------------------

    ! determine data layout
    nvec = mvec%jmax-mvec%jmin+1
    nblocks = size(block_list)
    nrows = mvec%map%nlocal(mvec%map%me)
    ldm = size(mvec%val,1)
    if( .not. mvec%is_view .or. &
      & ( mvec%jmin .eq. lbound(mvec%val,1) .and. &
      &   mvec%jmax .eq. ubound(mvec%val,1)       ) ) then
      strided = .false.
    else
      strided = .true.
    end if


    ! fast versions for fixed nvec and scattering single vectors
    single_vector_scatter = .true.
    if( nvec .ne. nblocks ) single_vector_scatter = .false.
    do i = 1, nblocks, 1
      if( block_list(i)%jmin .ne. block_list(i)%jmax ) then
        single_vector_scatter = .false.
      end if
    end do

    if( single_vector_scatter .and. .not. strided ) then
      if( nvec .eq. 1 ) then
        call dscatter_1(nrows, mvec%val(mvec%jmin,1), &
          &            block_list(1)%val(block_list(1)%jmin,1), size(block_list(1)%val,1) )
        return
      else if( nvec .eq. 2 ) then
        call dscatter_2(nrows, mvec%val(mvec%jmin,1), &
          &            block_list(1)%val(block_list(1)%jmin,1), size(block_list(1)%val,1), &
          &            block_list(1)%val(block_list(2)%jmin,1), size(block_list(2)%val,1) )
        return
      else if( nvec .eq. 4 ) then
        call dscatter_4(nrows, mvec%val(mvec%jmin,1), &
          &            block_list(1)%val(block_list(1)%jmin,1), size(block_list(1)%val,1), &
          &            block_list(2)%val(block_list(2)%jmin,1), size(block_list(2)%val,1), &
          &            block_list(3)%val(block_list(3)%jmin,1), size(block_list(3)%val,1), &
          &            block_list(4)%val(block_list(4)%jmin,1), size(block_list(4)%val,1) )
        return
      end if
    end if

!$omp parallel do private(off) schedule(static)
    do j = 1, nrows, 1
      off = mvec%jmin
      do i = 1, nblocks, 1
        jmin_i = block_list(i)%jmin
        jmax_i = block_list(i)%jmax
        nvec_i = jmax_i-jmin_i+1
        block_list(i)%val(jmin_i:jmax_i,j) = mvec%val(off:off+nvec_i-1,j)
        off = off +nvec_i
      end do
    end do

    !--------------------------------------------------------------------------------
  end subroutine mvec_scatter_mvecs


  subroutine mvec_random(mvec)
    !--------------------------------------------------------------------------------
    type(MVec_t), intent(inout) :: mvec
    !--------------------------------------------------------------------------------
    integer :: nvec, nrows, ldx
    integer(kind=8) :: pre_skip, post_skip
    !--------------------------------------------------------------------------------

    nrows = mvec%map%nlocal(mvec%map%me)
    nvec = mvec%jmax-mvec%jmin+1
    ldx = size(mvec%val,1)
    pre_skip = mvec%map%distrib(mvec%map%me)*nvec-1
    post_skip = (mvec%map%distrib(mvec%map%nProcs)-mvec%map%distrib(mvec%map%me+1))*nvec

    if( nvec .eq. ldx ) then
      call drandom_1(nvec*nrows, mvec%val(1,1), pre_skip, post_skip)
    else
      call drandom_general(nvec, nrows, mvec%val(mvec%jmin,1), ldx, pre_skip, post_skip)
    end if

  end subroutine mvec_random


  !==================================================================================
  !> calculate 2-norm of the vectors in a multivector
  subroutine mvec_norm2(mvec, vnrm, iflag)
    use mpi
    !--------------------------------------------------------------------------------
    type(MVec_t), intent(in)    :: mvec
    real(kind=8), intent(out)   :: vnrm(mvec%jmin:mvec%jmax)
    integer,      intent(inout) :: iflag
    !--------------------------------------------------------------------------------
    integer :: nvec, nrows, lda
    logical :: strided
    real(kind=8) :: localNrm(mvec%jmin:mvec%jmax)
    real(kind=8) :: localDot_prec(mvec%jmin:mvec%jmax,2)
    real(kind=8) :: globalDot_prec(mvec%jmin:mvec%jmax,2,mvec%map%nProcs)
    real(kind=8) :: globalDot_prec_(mvec%jmin:mvec%jmax,mvec%map%nProcs,2)
    !--------------------------------------------------------------------------------

    ! determine data layout
    if( .not. mvec%is_view .or. &
      & ( mvec%jmin .eq. lbound(mvec%val,1) .and. &
      &   mvec%jmax .eq. ubound(mvec%val,1)       ) ) then
      strided = .false.
    else
      strided = .true.
    end if

    nvec = mvec%jmax-mvec%jmin+1
    nrows = mvec%map%nlocal(mvec%map%me)
    lda = size(mvec%val,1)


    ! check if we need higher precision
    if(  FORCE_HIGH_PRECISION .or. iand(iflag,PHIST_ROBUST_REDUCTIONS) .gt. 0 ) then
#ifndef PHIST_HIGH_PRECISION_KERNELS
      iflag = PHIST_NOT_IMPLEMENTED
      return
#else
      ! check if we can do it
      if( strided ) then
        iflag = PHIST_NOT_IMPLEMENTED
        return
      end if
      if( nvec .ne. 1 .and. nvec .ne. 2 .and. nvec .ne. 4 ) then
        iflag = PHIST_NOT_IMPLEMENTED
        return
      end if

      if( nvec .eq. 1 ) then
        call ddot_self_prec_1(mvec%paddedN, mvec%val(mvec%jmin,1), localDot_prec(mvec%jmin,1), localDot_prec(mvec%jmin,2))
      else if( nvec .eq. 2 ) then
        call ddot_self_prec_2(mvec%paddedN, mvec%val(mvec%jmin,1), localDot_prec(mvec%jmin,1), localDot_prec(mvec%jmin,2))
      else if( nvec .eq. 4 ) then
        call ddot_self_prec_4(mvec%paddedN, mvec%val(mvec%jmin,1), localDot_prec(mvec%jmin,1), localDot_prec(mvec%jmin,2))
      end if

      call MPI_Allgather(localDot_prec,2*nvec,MPI_DOUBLE_PRECISION, &
        &                globalDot_prec,2*nvec,MPI_DOUBLE_PRECISION, &
        &                mvec%map%comm, iflag)

      ! rearrange received data
      globalDot_prec_(:,:,1) = globalDot_prec(:,1,:)
      globalDot_prec_(:,:,2) = globalDot_prec(:,2,:)
      if( nvec .eq. 1 ) then
        call prec_reduction_1(mvec%map%nProcs, globalDot_prec_(mvec%jmin,1,1), globalDot_prec_(mvec%jmin,1,2), &
          &                                    localDot_prec(mvec%jmin,1), localDot_prec(mvec%jmin,2))
      else if( nvec .eq. 2 ) then
        call prec_reduction_2(mvec%map%nProcs, globalDot_prec_(mvec%jmin,1,1), globalDot_prec_(mvec%jmin,1,2), &
          &                                    localDot_prec(mvec%jmin,1), localDot_prec(mvec%jmin,2))
      else if( nvec .eq. 4 ) then
        call prec_reduction_4(mvec%map%nProcs, globalDot_prec_(mvec%jmin,1,1), globalDot_prec_(mvec%jmin,1,2), &
          &                                    localDot_prec(mvec%jmin,1), localDot_prec(mvec%jmin,2))
      end if
      vnrm = sqrt(localDot_prec(:,1) + localDot_prec(:,2))

      return
#endif
    end if

    ! for single vectors call appropriate blas
    if( nvec .eq. 1 ) then
      if( strided ) then
        call dnrm2_strided_1(nrows, mvec%val(mvec%jmin,1), lda, vnrm)
      else
        call dnrm2_1(nrows, mvec%val, vnrm)
      end if
    else if( nvec .eq. 2 ) then
      if( strided ) then
        call dnrm2_strided_2(nrows, mvec%val(mvec%jmin,1), lda, vnrm)
      else
        call dnrm2_2(nrows, mvec%val, vnrm)
      end if
    else if( nvec .eq. 4 ) then
      if( strided ) then
        call dnrm2_strided_4(nrows, mvec%val(mvec%jmin,1), lda, vnrm)
      else
        call dnrm2_4(nrows, mvec%val, vnrm)
      end if
    else if( nvec .eq. 8 ) then
      if( strided ) then
        call dnrm2_strided_8(nrows, mvec%val(mvec%jmin,1), lda, vnrm)
      else
        call dnrm2_8(nrows, mvec%val, vnrm)
      end if
    else
      call dnrm2_general(nrows, nvec, mvec%val(mvec%jmin,1), lda, vnrm)
    end if

    localNrm = vnrm*vnrm
    call MPI_Allreduce(localNrm,vnrm,nvec,MPI_DOUBLE_PRECISION,MPI_SUM,mvec%map%comm,iflag)
    vnrm = sqrt(vnrm)

    !--------------------------------------------------------------------------------
  end subroutine mvec_norm2


  !==================================================================================
  ! scale mvec
  subroutine mvec_scale(x,alpha)
    !--------------------------------------------------------------------------------
    type(MVec_t), intent(inout) :: x
    real(kind=8), intent(in)    :: alpha
    !--------------------------------------------------------------------------------
    integer :: nvec
    real(kind=8), allocatable :: alpha_vec(:)
    !--------------------------------------------------------------------------------

    if( alpha .eq. 1 ) return

    nvec = x%jmax-x%jmin+1

    allocate(alpha_vec(nvec))
    alpha_vec = alpha
    call mvec_vscale(x,alpha_vec)

    !--------------------------------------------------------------------------------
  end subroutine mvec_scale


  !==================================================================================
  ! scale mvec
  subroutine mvec_vscale(x,alpha)
    !--------------------------------------------------------------------------------
    type(MVec_t), intent(inout) :: x
    real(kind=8), intent(in)    :: alpha(*)
    !--------------------------------------------------------------------------------
    integer :: nvec, nrows, lda
    logical :: strided
    !--------------------------------------------------------------------------------

    ! determine data layout
    if( .not. x%is_view .or. &
      & ( x%jmin .eq. lbound(x%val,1) .and. &
      &   x%jmax .eq. ubound(x%val,1)       ) ) then
      strided = .false.
    else
      strided = .true.
    end if

    nvec = x%jmax-x%jmin+1
    nrows = x%map%nlocal(x%map%me)
    lda = size(x%val,1)

    !call dlascl2(nvec,nrows,alpha(1),x%val(x%jmin,1),lda)

    if( nvec .eq. 1 ) then
      if( strided ) then
        call dscal_strided_1(nrows,alpha,x%val(x%jmin,1),lda)
      else
        call dscal_1(nrows,alpha,x%val)
      end if
    else if( nvec .eq. 2 ) then
      if( strided ) then
        call dscal_strided_2(nrows,alpha,x%val(x%jmin,1),lda)
      else
        call dscal_2(nrows,alpha,x%val)
      end if
    else if( nvec .eq. 4 ) then
      if( strided ) then
        call dscal_strided_4(nrows,alpha,x%val(x%jmin,1),lda)
      else
        call dscal_4(nrows,alpha,x%val)
      end if
    else if( nvec .eq. 8 ) then
      if( strided ) then
        call dscal_strided_8(nrows,alpha,x%val(x%jmin,1),lda)
      else
        call dscal_8(nrows,alpha,x%val)
      end if
    else
      call dscal_general(nrows, nvec, alpha, x%val(x%jmin,1), lda)
    end if
    !--------------------------------------------------------------------------------
  end subroutine mvec_vscale


  !==================================================================================
  ! daxpy routine for mvecs, delegate to mvec_vadd_mvec
  subroutine mvec_add_mvec(alpha,x,beta,y)
    !--------------------------------------------------------------------------------
    real(kind=8), intent(in)    :: alpha
    type(MVec_t), intent(in)    :: x
    real(kind=8), intent(in)    :: beta
    type(MVec_t), intent(inout) :: y
    !--------------------------------------------------------------------------------
    real(kind=8), allocatable :: alpha_vec(:)
    !--------------------------------------------------------------------------------

    if( alpha .eq. 0 ) then
      call mvec_scale(y,beta)
    else
      allocate(alpha_vec(y%jmin:y%jmax))
      alpha_vec = alpha
      call mvec_vadd_mvec(alpha_vec,x,beta,y)
    end if

    !--------------------------------------------------------------------------------
  end subroutine mvec_add_mvec


  !==================================================================================
  ! daxpy routine for mvecs
  subroutine mvec_vadd_mvec(alpha,x,beta,y)
    !--------------------------------------------------------------------------------
    real(kind=8), intent(in)    :: alpha(*)
    type(MVec_t), intent(in)    :: x
    real(kind=8), intent(in)    :: beta
    type(MVec_t), intent(inout) :: y
    !--------------------------------------------------------------------------------
    integer :: nvec, nrows, ldx, ldy
    logical :: strided_x, strided_y, strided
    logical :: only_scale, only_copy, one_alpha, y_aligned
    integer :: i
    !--------------------------------------------------------------------------------

    ! first gather some data to decide what we want
    nvec = y%jmax-y%jmin+1
    nrows = y%map%nlocal(y%map%me)
    ldy = size(y%val,1)

    if( .not. y%is_view .or. &
      & ( y%jmin .eq. lbound(y%val,1) .and. &
      &   y%jmax .eq. ubound(y%val,1)       ) ) then
      strided_y = .false.
    else
      strided_y = .true.
    end if


    only_scale = .true.
    do i = 1, nvec, 1
      if( alpha(i) .ne. 0 ) only_scale = .false.
    end do
    if( only_scale ) then
      call mvec_scale(y,beta)
      return
    end if

    only_copy = .true.
    if( beta .ne. 0 ) only_copy = .false.
    do i = 1, nvec, 1
      if( alpha(i) .ne. 1 ) only_copy = .false.
    end do

    ldx = size(x%val,1)

    if( .not. x%is_view .or. &
      & ( x%jmin .eq. lbound(x%val,1) .and. &
      &   x%jmax .eq. ubound(x%val,1)       ) ) then
      strided_x = .false.
    else
      strided_x = .true.
    end if

    one_alpha = .true.
    do i = 1, nvec-1, 1
      if( alpha(i) .ne. alpha(i+1) ) one_alpha = .false.
    end do

    y_aligned = .true.
    if( mod(loc(y%val(y%jmin,1)),16) .ne. 0 .or. mod(ldy,2) .ne. 0 ) then
      y_aligned = .false.
    end if

    strided = strided_x .or. strided_y

    if( only_copy ) then
      if( .not. strided ) then
        call dcopy_1(nvec*nrows, x%val, y%val)
      else
        call dcopy_general(nvec,nrows,x%val(x%jmin,1), ldx, y%val(y%jmin,1), ldy)
      end if
      return
    end if

#ifdef PHIST_HAVE_SSE
    if( beta .eq. 0 ) then
      if( nvec .eq. 2 .and. y_aligned ) then
        if( strided ) then
          call daxpy_NT_strided_2(nrows, alpha(1), x%val(x%jmin,1), ldx, y%val(y%jmin,1), ldy)
        else
          call daxpy_NT_2(nrows, alpha(1), x%val(x%jmin,1), y%val(y%jmin,1))
        end if
        return
      else if( nvec .eq. 4 .and. y_aligned ) then
        if( strided ) then
          call daxpy_NT_strided_4(nrows, alpha(1), x%val(x%jmin,1), ldx, y%val(y%jmin,1), ldy)
        else
          call daxpy_NT_4(nrows, alpha(1), x%val(x%jmin,1), y%val(y%jmin,1))
        end if
        return
      else if( nvec .eq. 8 .and. y_aligned ) then
        if( strided ) then
          call daxpy_NT_strided_8(nrows, alpha(1), x%val(x%jmin,1), ldx, y%val(y%jmin,1), ldy)
        else
          call daxpy_NT_8(nrows, alpha(1), x%val(x%jmin,1), y%val(y%jmin,1))
        end if
        return
      end if
    end if
#endif
    if( nvec .eq. 1 ) then
      if( strided ) then
        call daxpby_strided_1(nrows, alpha(1), x%val(x%jmin,1), ldx, beta, y%val(y%jmin,1), ldy)
      else
        call daxpby_1(nrows, alpha(1), x%val(x%jmin,1), beta, y%val(y%jmin,1))
      end if
      return
    else if( nvec .eq. 2 ) then
      if( strided ) then
        call daxpby_strided_2(nrows, alpha(1), x%val(x%jmin,1), ldx, beta, y%val(y%jmin,1), ldy)
      else
        call daxpby_2(nrows, alpha(1), x%val(x%jmin,1), beta, y%val(y%jmin,1))
      end if
      return
    else if( nvec .eq. 4 ) then
      if( strided ) then
        call daxpby_strided_4(nrows, alpha(1), x%val(x%jmin,1), ldx, beta, y%val(y%jmin,1), ldy)
      else
        call daxpby_4(nrows, alpha(1), x%val(x%jmin,1), beta, y%val(y%jmin,1))
      end if
      return
    else if( nvec .eq. 8 ) then
      if( strided ) then
        call daxpby_strided_8(nrows, alpha(1), x%val(x%jmin,1), ldx, beta, y%val(y%jmin,1), ldy)
      else
        call daxpby_8(nrows, alpha(1), x%val(x%jmin,1), beta, y%val(y%jmin,1))
      end if
      return
    end if

    call daxpby_generic(nrows, nvec, alpha(1), x%val(x%jmin,1), ldx, beta, y%val(y%jmin,1), ldy)
    !--------------------------------------------------------------------------------
  end subroutine mvec_vadd_mvec


  !==================================================================================
  ! dot product for mvecs
  subroutine mvec_dot_mvec(x,y,dot,iflag)
    use mpi
    !--------------------------------------------------------------------------------
    type(MVec_t), intent(in)    :: x, y
    real(kind=8), intent(out)   :: dot(x%jmin:x%jmax)
    integer,      intent(inout) :: iflag
    !--------------------------------------------------------------------------------
    integer :: nvec, nrows, ldx, ldy, ierr
    logical :: strided_x, strided_y, strided
    real(kind=8) :: localDot(x%jmin:x%jmax)
    real(kind=8) :: localDot_prec(x%jmin:x%jmax,2)
    real(kind=8) :: globalDot_prec(x%jmin:x%jmax,2,x%map%nProcs)
    real(kind=8) :: globalDot_prec_(x%jmin:x%jmax,x%map%nProcs,2)
    !--------------------------------------------------------------------------------

    ! determine data layout
    if( .not. x%is_view .or. &
      & ( x%jmin .eq. lbound(x%val,1) .and. &
      &   x%jmax .eq. ubound(x%val,1)       ) ) then
      strided_x = .false.
    else
      strided_x = .true.
    end if

    if( .not. y%is_view .or. &
      & ( y%jmin .eq. lbound(y%val,1) .and. &
      &   y%jmax .eq. ubound(y%val,1)       ) ) then
      strided_y = .false.
    else
      strided_y = .true.
    end if

    strided = strided_x .or. strided_y

    nvec = x%jmax-x%jmin+1
    nrows = x%map%nlocal(x%map%me)
    ldx = size(x%val,1)
    ldy = size(y%val,1)

    ! check if we need higher precision
    if(  FORCE_HIGH_PRECISION .or.iand(iflag,PHIST_ROBUST_REDUCTIONS) .gt. 0 ) then
#ifndef PHIST_HIGH_PRECISION_KERNELS
      iflag = PHIST_NOT_IMPLEMENTED
      return
#else
      ! check if we can do it
      if( strided ) then
        iflag = PHIST_NOT_IMPLEMENTED
        return
      end if
      if( nvec .ne. 1 .and. nvec .ne. 2 .and. nvec .ne. 4 ) then
        iflag = PHIST_NOT_IMPLEMENTED
        return
      end if

      if( nvec .eq. 1 ) then
        call ddot_prec_1(x%paddedN, x%val(x%jmin,1), y%val(y%jmin,1), localDot_prec(x%jmin,1), localDot_prec(x%jmin,2))
      else if( nvec .eq. 2 ) then
        call ddot_prec_2(x%paddedN, x%val(x%jmin,1), y%val(y%jmin,1), localDot_prec(x%jmin,1), localDot_prec(x%jmin,2))
      else if( nvec .eq. 4 ) then
        call ddot_prec_4(x%paddedN, x%val(x%jmin,1), y%val(y%jmin,1), localDot_prec(x%jmin,1), localDot_prec(x%jmin,2))
      end if

      call MPI_Allgather(localDot_prec,2*nvec,MPI_DOUBLE_PRECISION, &
        &                globalDot_prec,2*nvec,MPI_DOUBLE_PRECISION, &
        &                x%map%comm, iflag)

      ! rearrange received data
      globalDot_prec_(:,:,1) = globalDot_prec(:,1,:)
      globalDot_prec_(:,:,2) = globalDot_prec(:,2,:)
      if( nvec .eq. 1 ) then
        call prec_reduction_1(x%map%nProcs, globalDot_prec_(x%jmin,1,1), globalDot_prec_(x%jmin,1,2), &
          &                                 localDot_prec(x%jmin,1), localDot_prec(x%jmin,2))
      else if( nvec .eq. 2 ) then
        call prec_reduction_2(x%map%nProcs, globalDot_prec_(x%jmin,1,1), globalDot_prec_(x%jmin,1,2), &
          &                                 localDot_prec(x%jmin,1), localDot_prec(x%jmin,2))
      else if( nvec .eq. 4 ) then
        call prec_reduction_4(x%map%nProcs, globalDot_prec_(x%jmin,1,1), globalDot_prec_(x%jmin,1,2), &
          &                                 localDot_prec(x%jmin,1), localDot_prec(x%jmin,2))
      end if
      dot = localDot_prec(:,1) + localDot_prec(:,2)

      return
#endif
    end if

    ! for single vectors call appropriate blas
    if( nvec .eq. 1 ) then
      if( strided ) then
        call ddot_strided_1(nrows, x%val(x%jmin,1), ldx, y%val(y%jmin,1), ldy, dot)
      else
        call ddot_1(nrows, x%val, y%val, dot)
      end if
    else if( nvec .eq. 2 ) then
      if( strided ) then
        call ddot_strided_2(nrows, x%val(x%jmin,1), ldx, y%val(y%jmin,1), ldy, dot)
      else
        call ddot_2(nrows, x%val, y%val, dot)
      end if
    else if( nvec .eq. 4 ) then
      if( strided ) then
        call ddot_strided_4(nrows, x%val(x%jmin,1), ldx, y%val(y%jmin,1), ldy, dot)
      else
        call ddot_4(nrows, x%val, y%val, dot)
      end if
    else if( nvec .eq. 8 ) then
      if( strided ) then
        call ddot_strided_8(nrows, x%val(x%jmin,1), ldx, y%val(y%jmin,1), ldy, dot)
      else
        call ddot_8(nrows, x%val, y%val, dot)
      end if
    else
      call ddot_general(nrows, nvec, x%val(x%jmin,1), ldx, y%val(y%jmin,1), ldy, dot)
    end if

    localDot = dot
    call MPI_Allreduce(localDot,dot,nvec,MPI_DOUBLE_PRECISION,MPI_SUM,x%map%comm,iflag)

    !--------------------------------------------------------------------------------
  end subroutine mvec_dot_mvec


  !==================================================================================
  ! special gemm routine for mvec times sdmat
  subroutine mvec_times_sdmat(alpha,v,M,beta,w,iflag)
    !--------------------------------------------------------------------------------
    real(kind=8),  intent(in)    :: alpha, beta
    type(MVec_t),  intent(in)    :: v
    type(SDMat_t), intent(in)    :: M
    type(Mvec_t),  intent(inout) :: w
    integer,       intent(inout) :: iflag
    !--------------------------------------------------------------------------------
    integer :: nrows, nvecv, nvecw, ldv, ldw
    logical :: strided_v, strided_w
    real(kind=8), allocatable :: Mtmp(:,:), MCtmp(:,:)
    !--------------------------------------------------------------------------------

    ! check if we only need to scale
    if( alpha .eq. 0 ) then
      call mvec_scale(w,beta)
      iflag = 0
      return
    end if

    ! determine data layout
    nrows = v%map%nlocal(v%map%me)
    nvecv = v%jmax-v%jmin+1
    nvecw = w%jmax-w%jmin+1
    ldv = size(v%val,1)
    ldw = size(w%val,1)
    if( .not. v%is_view .or. &
      & ( v%jmin .eq. lbound(v%val,1) .and. &
      &   v%jmax .eq. ubound(v%val,1)       ) ) then
      strided_v = .false.
    else
      strided_v = .true.
    end if

    if( .not. w%is_view .or. &
      & ( w%jmin .eq. lbound(w%val,1) .and. &
      &   w%jmax .eq. ubound(w%val,1)       ) ) then
      strided_w = .false.
    else
      strided_w = .true.
    end if


    if(  FORCE_HIGH_PRECISION .or.iand(iflag,PHIST_ROBUST_REDUCTIONS) .gt. 0 ) then
#ifndef PHIST_HIGH_PRECISION_KERNELS
      iflag = PHIST_NOT_IMPLEMENTED
      return
#else
      ! check if we can do it
      if( strided_w .or. (nvecw .ne. 1 .and. nvecw .ne. 2 .and. nvecw .ne. 4) ) then
        iflag = PHIST_NOT_IMPLEMENTED
        return
      end if

      ! copy M to buffer
      allocate(Mtmp(nvecv,nvecw),MCtmp(nvecv,nvecw))
      Mtmp = M%val(M%imin:M%imax,M%jmin:M%jmax)
      MCtmp = M%err(M%imin:M%imax,M%jmin:M%jmax)

      ! for beta=0 we can use nontemporal stores
      if( beta .eq. 0. ) then
        if( nvecw .eq. 4 ) then
          if( strided_v ) then
            call dgemm_sb_prec_k_strided_4_nt(v%paddedN,nvecv,alpha,v%val(v%jmin,1),ldv,Mtmp,MCtmp,w%val(w%jmin,1))
          else
            call dgemm_sb_prec_k_4_nt(v%paddedN,nvecv,alpha,v%val(v%jmin,1),Mtmp,MCtmp,w%val(w%jmin,1))
          end if
        else if( nvecw .eq. 2 ) then
          if( strided_v ) then
            call dgemm_sb_prec_k_strided_2_nt(v%paddedN,nvecv,alpha,v%val(v%jmin,1),ldv,Mtmp,MCtmp,w%val(w%jmin,1))
          else
            call dgemm_sb_prec_k_2_nt(v%paddedN,nvecv,alpha,v%val(v%jmin,1),Mtmp,MCtmp,w%val(w%jmin,1))
          end if
        else ! nvecw .eq. 1
          if( strided_v ) then
            call dgemm_sb_prec_k_strided_1_nt(v%paddedN,nvecv,alpha,v%val(v%jmin,1),ldv,Mtmp,MCtmp,w%val(w%jmin,1))
          else
            call dgemm_sb_prec_k_1_nt(v%paddedN,nvecv,alpha,v%val(v%jmin,1),Mtmp,MCtmp,w%val(w%jmin,1))
          end if
        end if
      else ! beta .ne. 0
        if( nvecw .eq. 4 ) then
          if( strided_v ) then
            call dgemm_sb_prec_k_strided_4(v%paddedN,nvecv,alpha,v%val(v%jmin,1),ldv,Mtmp,MCtmp,beta,w%val(w%jmin,1))
          else
            call dgemm_sb_prec_k_4(v%paddedN,nvecv,alpha,v%val(v%jmin,1),Mtmp,MCtmp,beta,w%val(w%jmin,1))
          end if
        else if( nvecw .eq. 2 ) then
          if( strided_v ) then
            call dgemm_sb_prec_k_strided_2(v%paddedN,nvecv,alpha,v%val(v%jmin,1),ldv,Mtmp,MCtmp,beta,w%val(w%jmin,1))
          else
            call dgemm_sb_prec_k_2(v%paddedN,nvecv,alpha,v%val(v%jmin,1),Mtmp,MCtmp,beta,w%val(w%jmin,1))
          end if
        else ! nvecw .eq. 1
          if( strided_v ) then
            call dgemm_sb_prec_k_strided_1(v%paddedN,nvecv,alpha,v%val(v%jmin,1),ldv,Mtmp,MCtmp,beta,w%val(w%jmin,1))
          else
            call dgemm_sb_prec_k_1(v%paddedN,nvecv,alpha,v%val(v%jmin,1),Mtmp,MCtmp,beta,w%val(w%jmin,1))
          end if
        end if
      end if
      iflag = 0
      return
#endif
    end if

    iflag = 0

    allocate(Mtmp(nvecw,nvecv))
    Mtmp = transpose(M%val(M%imin:M%imax,M%jmin:M%jmax))
    ! recognize small block mvecs
    if( .not. strided_w ) then
      if( nvecw .eq. 1 ) then
        if( strided_v ) then
          call dgemm_sB_strided_1_k(nrows, nvecv, alpha, v%val(v%jmin,1), ldv, Mtmp, beta, w%val)
        else
          call dgemm_sB_1_k        (nrows, nvecv, alpha, v%val, Mtmp, beta, w%val)
        end if
        return
      else if( nvecw .eq. 2 ) then
        if( strided_v ) then
          call dgemm_sB_strided_2_k(nrows, nvecv, alpha, v%val(v%jmin,1), ldv, Mtmp, beta, w%val)
        else
          call dgemm_sB_2_k        (nrows, nvecv, alpha, v%val, Mtmp, beta, w%val)
        end if
        return
      else if( nvecw .eq. 4 ) then
        if( strided_v ) then
          call dgemm_sB_strided_4_k(nrows, nvecv, alpha, v%val(v%jmin,1), ldv, Mtmp, beta, w%val)
        else
          call dgemm_sB_4_k        (nrows, nvecv, alpha, v%val, Mtmp, beta, w%val)
        end if
        return
      else if( nvecw .eq. 8 ) then
        if( strided_v ) then
          call dgemm_sB_strided_8_k(nrows, nvecv, alpha, v%val(v%jmin,1), ldv, Mtmp, beta, w%val)
        else
          call dgemm_sB_8_k        (nrows, nvecv, alpha, v%val, Mtmp, beta, w%val)
        end if
        return
      end if
    end if
    deallocate(Mtmp)
    allocate(Mtmp(nvecv,nvecw))
    Mtmp = M%val(M%imin:M%imax,M%jmin:M%jmax)

    if( .not. strided_v ) then
      if( nvecv .eq. 1 ) then
        if( strided_w ) then
          call dgemm_sB_strided_k_1(nrows, nvecw, alpha, v%val, Mtmp, beta, w%val(w%jmin,1), ldw)
        else
          call dgemm_sB_k_1        (nrows, nvecw, alpha, v%val, Mtmp, beta, w%val)
        end if
        return
      else if( nvecv .eq. 2 ) then
        if( strided_w ) then
          call dgemm_sB_strided_k_2(nrows, nvecw, alpha, v%val, Mtmp, beta, w%val(w%jmin,1), ldw)
        else
          call dgemm_sB_k_2        (nrows, nvecw, alpha, v%val, Mtmp, beta, w%val)
        end if
        return
      else if( nvecv .eq. 4 ) then
        if( strided_w ) then
          call dgemm_sB_strided_k_4(nrows, nvecw, alpha, v%val, Mtmp, beta, w%val(w%jmin,1), ldw)
        else
          call dgemm_sB_k_4        (nrows, nvecw, alpha, v%val, Mtmp, beta, w%val)
        end if
        return
      else if( nvecv .eq. 8 ) then
        if( strided_w ) then
          call dgemm_sB_strided_k_8(nrows, nvecw, alpha, v%val, Mtmp, beta, w%val(w%jmin,1), ldw)
        else
          call dgemm_sB_k_8        (nrows, nvecw, alpha, v%val, Mtmp, beta, w%val)
        end if
        return
      end if
    end if


    call dgemm_sB_generic(nrows,nvecw,nvecv,alpha,v%val(v%jmin,1),ldv, Mtmp, beta, w%val(w%jmin,1),ldw)


    !--------------------------------------------------------------------------------
  end subroutine mvec_times_sdmat

  !==================================================================================
  ! special gemm routine for the augmented mvec times sdmat variant
  subroutine mvec_times_sdmat_augmented(alpha,v,M,beta,w,N,iflag)
    use mpi
    !--------------------------------------------------------------------------------
    real(kind=8),  intent(in)    :: alpha, beta
    type(MVec_t),  intent(in)    :: v
    type(SDMat_t), intent(in)    :: M
    type(Mvec_t),  intent(inout) :: w
    type(SDMat_t), intent(inout) :: N
    integer,       intent(inout) :: iflag
    !--------------------------------------------------------------------------------
    integer :: nrows, nvecv, nvecw, ldv, ldw
    logical :: strided_v, strided_w, handled
    real(kind=8), allocatable :: Mtmp(:,:), MCtmp(:,:)
    real(kind=8), allocatable :: Ntmp(:,:), NCtmp(:,:)
    real(kind=8), allocatable :: localBuff(:,:,:), globalBuff(:,:,:,:), globalBuff_(:,:,:,:)
    !--------------------------------------------------------------------------------

    ! determine data layout
    nrows = v%map%nlocal(v%map%me)
    nvecv = v%jmax-v%jmin+1
    nvecw = w%jmax-w%jmin+1
    ldv = size(v%val,1)
    ldw = size(w%val,1)
    if( .not. v%is_view .or. &
      & ( v%jmin .eq. lbound(v%val,1) .and. &
      &   v%jmax .eq. ubound(v%val,1)       ) ) then
      strided_v = .false.
    else
      strided_v = .true.
    end if

    if( .not. w%is_view .or. &
      & ( w%jmin .eq. lbound(w%val,1) .and. &
      &   w%jmax .eq. ubound(w%val,1)       ) ) then
      strided_w = .false.
    else
      strided_w = .true.
    end if


    if(  FORCE_HIGH_PRECISION .or.iand(iflag,PHIST_ROBUST_REDUCTIONS) .gt. 0 ) then
#ifndef PHIST_HIGH_PRECISION_KERNELS
      iflag = PHIST_NOT_IMPLEMENTED
      return
#else
      ! check if we can do it
      if( strided_w .or. (nvecw .ne. 1 .and. nvecw .ne. 2 .and. nvecw .ne. 4) ) then
        iflag = PHIST_NOT_IMPLEMENTED
        return
      end if

      ! copy M to buffer
      allocate(Mtmp(nvecv,nvecw),MCtmp(nvecv,nvecw))
      Mtmp = M%val(M%imin:M%imax,M%jmin:M%jmax)
      MCtmp = M%err(M%imin:M%imax,M%jmin:M%jmax)
      allocate(Ntmp(nvecw,nvecw),NCtmp(nvecw,nvecw))

      ! for beta=0 we can use nontemporal stores
      if( beta .eq. 0. ) then
        if( nvecw .eq. 4 ) then
          if( strided_v ) then
            call dgemm_sb_augmented_prec_strided_k_4_nt(v%paddedN,nvecv,alpha,v%val(v%jmin,1),ldv,Mtmp,MCtmp,w%val(w%jmin,1),Ntmp,NCtmp)
          else
            if( nvecv .eq. 1 ) then
              call dgemm_sb_augmented_prec_1_4_nt(v%paddedN,alpha,v%val(v%jmin,1),Mtmp,MCtmp,w%val(w%jmin,1),Ntmp,NCtmp)
            else if( nvecv .eq. 2 ) then
              call dgemm_sb_augmented_prec_2_4_nt(v%paddedN,alpha,v%val(v%jmin,1),Mtmp,MCtmp,w%val(w%jmin,1),Ntmp,NCtmp)
            else if( nvecv .eq. 4 ) then
              call dgemm_sb_augmented_prec_4_4_nt(v%paddedN,alpha,v%val(v%jmin,1),Mtmp,MCtmp,w%val(w%jmin,1),Ntmp,NCtmp)
            else
              call dgemm_sb_augmented_prec_k_4_nt(v%paddedN,nvecv,alpha,v%val(v%jmin,1),Mtmp,MCtmp,w%val(w%jmin,1),Ntmp,NCtmp)
            end if
          end if
        else if( nvecw .eq. 2 ) then
          if( strided_v ) then
            call dgemm_sb_augmented_prec_strided_k_2_nt(v%paddedN,nvecv,alpha,v%val(v%jmin,1),ldv,Mtmp,MCtmp,w%val(w%jmin,1),Ntmp,NCtmp)
          else
            if( nvecv .eq. 1 ) then
              call dgemm_sb_augmented_prec_1_2_nt(v%paddedN,alpha,v%val(v%jmin,1),Mtmp,MCtmp,w%val(w%jmin,1),Ntmp,NCtmp)
            else if( nvecv .eq. 2 ) then
              call dgemm_sb_augmented_prec_2_2_nt(v%paddedN,alpha,v%val(v%jmin,1),Mtmp,MCtmp,w%val(w%jmin,1),Ntmp,NCtmp)
            else if( nvecv .eq. 4 ) then
              call dgemm_sb_augmented_prec_4_2_nt(v%paddedN,alpha,v%val(v%jmin,1),Mtmp,MCtmp,w%val(w%jmin,1),Ntmp,NCtmp)
            else
              call dgemm_sb_augmented_prec_k_2_nt(v%paddedN,nvecv,alpha,v%val(v%jmin,1),Mtmp,MCtmp,w%val(w%jmin,1),Ntmp,NCtmp)
            end if
          end if
        else ! nvecw .eq. 1
          if( strided_v ) then
            call dgemm_sb_augmented_prec_strided_k_1_nt(v%paddedN,nvecv,alpha,v%val(v%jmin,1),ldv,Mtmp,MCtmp,w%val(w%jmin,1),Ntmp,NCtmp)
          else
            if( nvecv .eq. 1 ) then
              call dgemm_sb_augmented_prec_1_1_nt(v%paddedN,alpha,v%val(v%jmin,1),Mtmp,MCtmp,w%val(w%jmin,1),Ntmp,NCtmp)
            else if( nvecv .eq. 2 ) then
              call dgemm_sb_augmented_prec_2_1_nt(v%paddedN,alpha,v%val(v%jmin,1),Mtmp,MCtmp,w%val(w%jmin,1),Ntmp,NCtmp)
            else if( nvecv .eq. 4 ) then
              call dgemm_sb_augmented_prec_4_1_nt(v%paddedN,alpha,v%val(v%jmin,1),Mtmp,MCtmp,w%val(w%jmin,1),Ntmp,NCtmp)
            else
              call dgemm_sb_augmented_prec_k_1_nt(v%paddedN,nvecv,alpha,v%val(v%jmin,1),Mtmp,MCtmp,w%val(w%jmin,1),Ntmp,NCtmp)
            end if
          end if
        end if
      else ! beta .ne. 0
        if( nvecw .eq. 4 ) then
          if( strided_v ) then
            call dgemm_sb_augmented_prec_strided_k_4(v%paddedN,nvecv,alpha,v%val(v%jmin,1),ldv,Mtmp,MCtmp,beta,w%val(w%jmin,1),Ntmp,NCtmp)
          else
            if( nvecv .eq. 1 ) then
              call dgemm_sb_augmented_prec_1_4(v%paddedN,alpha,v%val(v%jmin,1),Mtmp,MCtmp,beta,w%val(w%jmin,1),Ntmp,NCtmp)
            else if( nvecv .eq. 2 ) then
              call dgemm_sb_augmented_prec_2_4(v%paddedN,alpha,v%val(v%jmin,1),Mtmp,MCtmp,beta,w%val(w%jmin,1),Ntmp,NCtmp)
            else if( nvecv .eq. 4 ) then
              call dgemm_sb_augmented_prec_4_4(v%paddedN,alpha,v%val(v%jmin,1),Mtmp,MCtmp,beta,w%val(w%jmin,1),Ntmp,NCtmp)
            else
              call dgemm_sb_augmented_prec_k_4(v%paddedN,nvecv,alpha,v%val(v%jmin,1),Mtmp,MCtmp,beta,w%val(w%jmin,1),Ntmp,NCtmp)
            end if
          end if
        else if( nvecw .eq. 2 ) then
          if( strided_v ) then
            call dgemm_sb_augmented_prec_strided_k_2(v%paddedN,nvecv,alpha,v%val(v%jmin,1),ldv,Mtmp,MCtmp,beta,w%val(w%jmin,1),Ntmp,NCtmp)
          else
            if( nvecv .eq. 1 ) then
              call dgemm_sb_augmented_prec_1_2(v%paddedN,alpha,v%val(v%jmin,1),Mtmp,MCtmp,beta,w%val(w%jmin,1),Ntmp,NCtmp)
            else if( nvecv .eq. 2 ) then
              call dgemm_sb_augmented_prec_2_2(v%paddedN,alpha,v%val(v%jmin,1),Mtmp,MCtmp,beta,w%val(w%jmin,1),Ntmp,NCtmp)
            else if( nvecv .eq. 4 ) then
              call dgemm_sb_augmented_prec_4_2(v%paddedN,alpha,v%val(v%jmin,1),Mtmp,MCtmp,beta,w%val(w%jmin,1),Ntmp,NCtmp)
            else
              call dgemm_sb_augmented_prec_k_2(v%paddedN,nvecv,alpha,v%val(v%jmin,1),Mtmp,MCtmp,beta,w%val(w%jmin,1),Ntmp,NCtmp)
            end if
          end if
        else ! nvecw .eq. 1
          if( strided_v ) then
            call dgemm_sb_augmented_prec_strided_k_1(v%paddedN,nvecv,alpha,v%val(v%jmin,1),ldv,Mtmp,MCtmp,beta,w%val(w%jmin,1),Ntmp,NCtmp)
          else
            if( nvecv .eq. 1 ) then
              call dgemm_sb_augmented_prec_1_1(v%paddedN,alpha,v%val(v%jmin,1),Mtmp,MCtmp,beta,w%val(w%jmin,1),Ntmp,NCtmp)
            else if( nvecv .eq. 2 ) then
              call dgemm_sb_augmented_prec_2_1(v%paddedN,alpha,v%val(v%jmin,1),Mtmp,MCtmp,beta,w%val(w%jmin,1),Ntmp,NCtmp)
            else if( nvecv .eq. 4 ) then
              call dgemm_sb_augmented_prec_4_1(v%paddedN,alpha,v%val(v%jmin,1),Mtmp,MCtmp,beta,w%val(w%jmin,1),Ntmp,NCtmp)
            else
              call dgemm_sb_augmented_prec_k_1(v%paddedN,nvecv,alpha,v%val(v%jmin,1),Mtmp,MCtmp,beta,w%val(w%jmin,1),Ntmp,NCtmp)
            end if
          end if
        end if

      end if
      ! gather results
      allocate(globalBuff(nvecw,nvecw,2,v%map%nProcs))
      allocate(globalBuff_(nvecw,nvecw,v%map%nProcs,2))
      allocate(localBuff(nvecw,nvecw,2))
      localBuff(:,:,1) = Ntmp
      localBuff(:,:,2) = NCtmp
      call MPI_Allgather(localBuff,2*nvecw*nvecw,MPI_DOUBLE_PRECISION, &
        &                globalBuff,2*nvecw*nvecw,MPI_DOUBLE_PRECISION, &
        &                v%map%comm, iflag)
      globalBuff_(:,:,:,1) = globalBuff(:,:,1,:)
      globalBuff_(:,:,:,2) = globalBuff(:,:,2,:)
!write(*,*) 'here globalBuff', globalBuff_(:,:,:,1)
!write(*,*) 'here globalBuffC', globalBuff_(:,:,:,2)
      ! MPI reduction
      if( nvecw*nvecw .eq. 1 ) then
        call prec_reduction_1(v%map%nProcs, globalBuff_(1,1,1,1), globalBuff_(1,1,1,2), &
          &                                 localBuff(1,1,1), localBuff(1,1,2)          )
      else if( nvecw*nvecw .eq. 4 ) then
        call prec_reduction_4(v%map%nProcs, globalBuff_(1,1,1,1), globalBuff_(1,1,1,2), &
          &                                 localBuff(1,1,1), localBuff(1,1,2)          )
      else if( mod(nvecw*nvecw,4) .eq. 0 ) then
        call prec_reduction_4k(v%map%nProcs, nvecw*nvecw, globalBuff_(1,1,1,1), globalBuff_(1,1,1,2), &
          &                                               localBuff(1,1,1), localBuff(1,1,2)          )
      else
        call prec_reduction_k(v%map%nProcs, nvecw*nvecw, globalBuff_(1,1,1,1), globalBuff_(1,1,1,2), &
          &                                              localBuff(1,1,1), localBuff(1,1,2)          )
      end if
!write(*,*) 'here vTw', localBuff(:,:,1)
!write(*,*) 'here vTw _C', localBuff(:,:,2)

      ! set result
      N%val(N%imin:N%imax,N%jmin:N%jmax) = localBuff(:,:,1)
      N%err(N%imin:N%imax,N%jmin:N%jmax) = localBuff(:,:,2)

      return
#endif
    end if


    handled = .false.
    allocate(Ntmp(nvecw,nvecw),NCtmp(nvecw,nvecw))
    if( .not. handled ) then
      allocate(Mtmp(nvecw,nvecv))
      Mtmp = transpose(M%val(M%imin:M%imax,M%jmin:M%jmax))
      ! recognize small block mvecs
      if( .not. strided_w ) then
        if( nvecw .eq. 1 ) then
          if( strided_v ) then
            call dgemm_sB_augmented_1_strided_k(nrows, nvecv, alpha, v%val(v%jmin,1), ldv, Mtmp, beta, w%val,Ntmp)
          else
            call dgemm_sB_augmented_1_k        (nrows, nvecv, alpha, v%val, Mtmp, beta, w%val,Ntmp)
          end if
          handled = .true.
        else if( nvecw .eq. 2 ) then
          if( strided_v ) then
            call dgemm_sB_augmented_2_strided_k(nrows, nvecv, alpha, v%val(v%jmin,1), ldv, Mtmp, beta, w%val,Ntmp)
          else
            call dgemm_sB_augmented_2_k        (nrows, nvecv, alpha, v%val, Mtmp, beta, w%val,Ntmp)
          end if
          handled = .true.
        else if( nvecw .eq. 4 ) then
          if( strided_v ) then
            call dgemm_sB_augmented_4_strided_k(nrows, nvecv, alpha, v%val(v%jmin,1), ldv, Mtmp, beta, w%val,Ntmp)
          else
            call dgemm_sB_augmented_4_k        (nrows, nvecv, alpha, v%val, Mtmp, beta, w%val,Ntmp)
          end if
          handled = .true.
        else if( nvecw .eq. 8 ) then
          if( strided_v ) then
            call dgemm_sB_augmented_8_strided_k(nrows, nvecv, alpha, v%val(v%jmin,1), ldv, Mtmp, beta, w%val,Ntmp)
          else
            call dgemm_sB_augmented_8_k        (nrows, nvecv, alpha, v%val, Mtmp, beta, w%val,Ntmp)
          end if
          handled = .true.
        end if
      end if
      deallocate(Mtmp)
    end if

    if( .not. handled ) then
      allocate(Mtmp(nvecv,nvecw))
      Mtmp = M%val(M%imin:M%imax,M%jmin:M%jmax)
      call dgemm_sB_augmented_generic(nrows,nvecw,nvecv,alpha,v%val(v%jmin,1),ldv, Mtmp, beta, w%val(w%jmin,1),ldw,Ntmp)
      deallocate(Mtmp)
    end if


    ! global reduction of result
    call MPI_Allreduce(Ntmp,NCtmp,nvecw*nvecw,MPI_DOUBLE_PRECISION,MPI_SUM,v%map%comm,iflag)

    N%val(N%imin:N%imax,N%jmin:N%jmax) = NCtmp
#ifdef PHIST_HIGH_PRECISION_KERNELS
    N%err(N%imin:N%imax,N%jmin:N%jmax) = 0.
#endif
    !--------------------------------------------------------------------------------
  end subroutine mvec_times_sdmat_augmented


  !==================================================================================
  ! special gemm routine for the augmented mvec times sdmat variant
  subroutine mvec_times_sdmat_add_mvec_times_sdMat(v,M,w,N,iflag)
    use mpi
    !--------------------------------------------------------------------------------
    type(MVec_t),  intent(in)    :: v
    type(SDMat_t), intent(in)    :: M
    type(Mvec_t),  intent(inout) :: w
    type(SDMat_t), intent(in) :: N
    integer,       intent(inout) :: iflag
    !--------------------------------------------------------------------------------
    integer :: nrows, nvecv, nvecw, ldv, ldw
    logical :: strided_v, strided_w, handled
    real(kind=8), allocatable :: Mtmp(:,:), MCtmp(:,:)
    real(kind=8), allocatable :: Ntmp(:,:), NCtmp(:,:)
    !--------------------------------------------------------------------------------

    ! determine data layout
    nrows = v%map%nlocal(v%map%me)
    nvecv = v%jmax-v%jmin+1
    nvecw = w%jmax-w%jmin+1
    ldv = size(v%val,1)
    ldw = size(w%val,1)
    if( .not. v%is_view .or. &
      & ( v%jmin .eq. lbound(v%val,1) .and. &
      &   v%jmax .eq. ubound(v%val,1)       ) ) then
      strided_v = .false.
    else
      strided_v = .true.
    end if

    if( .not. w%is_view .or. &
      & ( w%jmin .eq. lbound(w%val,1) .and. &
      &   w%jmax .eq. ubound(w%val,1)       ) ) then
      strided_w = .false.
    else
      strided_w = .true.
    end if


    if(  FORCE_HIGH_PRECISION .or.iand(iflag,PHIST_ROBUST_REDUCTIONS) .gt. 0 ) then
#ifndef PHIST_HIGH_PRECISION_KERNELS
      iflag = PHIST_NOT_IMPLEMENTED
      return
#else
      ! check if we can do it
      if( strided_w .or. (nvecw .ne. 1 .and. nvecw .ne. 2 .and. nvecw .ne. 4) ) then
        iflag = PHIST_NOT_IMPLEMENTED
        return
      end if

      ! copy M to buffer
      allocate(Mtmp(nvecv,nvecw),MCtmp(nvecv,nvecw))
      Mtmp = M%val(M%imin:M%imax,M%jmin:M%jmax)
      MCtmp = M%err(M%imin:M%imax,M%jmin:M%jmax)
      allocate(Ntmp(nvecw,nvecw),NCtmp(nvecw,nvecw))
      Ntmp = N%val(N%imin:N%imax,N%jmin:N%jmax)
      NCtmp = N%err(N%imin:N%imax,N%jmin:N%jmax)

      if( nvecw .eq. 4 ) then
        if( strided_v ) then
          call dgemm_sb_add_sd_prec_strided_k_4(v%paddedN,nvecv,v%val(v%jmin,1),ldv,Mtmp,MCtmp,w%val(w%jmin,1),Ntmp,NCtmp)
        else
          if( nvecv .eq. 1 ) then
            call dgemm_sb_add_sd_prec_1_4(v%paddedN,v%val(v%jmin,1),Mtmp,MCtmp,w%val(w%jmin,1),Ntmp,NCtmp)
          else if( nvecv .eq. 2 ) then
            call dgemm_sb_add_sd_prec_2_4(v%paddedN,v%val(v%jmin,1),Mtmp,MCtmp,w%val(w%jmin,1),Ntmp,NCtmp)
          else if( nvecv .eq. 4 ) then
            call dgemm_sb_add_sd_prec_4_4(v%paddedN,v%val(v%jmin,1),Mtmp,MCtmp,w%val(w%jmin,1),Ntmp,NCtmp)
          else
            call dgemm_sb_add_sd_prec_k_4(v%paddedN,nvecv,v%val(v%jmin,1),Mtmp,MCtmp,w%val(w%jmin,1),Ntmp,NCtmp)
          end if
        end if
      else if( nvecw .eq. 2 ) then
        if( strided_v ) then
          call dgemm_sb_add_sd_prec_strided_k_2(v%paddedN,nvecv,v%val(v%jmin,1),ldv,Mtmp,MCtmp,w%val(w%jmin,1),Ntmp,NCtmp)
        else
          if( nvecv .eq. 1 ) then
            call dgemm_sb_add_sd_prec_1_2(v%paddedN,v%val(v%jmin,1),Mtmp,MCtmp,w%val(w%jmin,1),Ntmp,NCtmp)
          else if( nvecv .eq. 2 ) then
            call dgemm_sb_add_sd_prec_2_2(v%paddedN,v%val(v%jmin,1),Mtmp,MCtmp,w%val(w%jmin,1),Ntmp,NCtmp)
          else if( nvecv .eq. 4 ) then
            call dgemm_sb_add_sd_prec_4_2(v%paddedN,v%val(v%jmin,1),Mtmp,MCtmp,w%val(w%jmin,1),Ntmp,NCtmp)
          else
            call dgemm_sb_add_sd_prec_k_2(v%paddedN,nvecv,v%val(v%jmin,1),Mtmp,MCtmp,w%val(w%jmin,1),Ntmp,NCtmp)
          end if
        end if
      else ! nvecw .eq. 1
        if( strided_v ) then
          call dgemm_sb_add_sd_prec_strided_k_1(v%paddedN,nvecv,v%val(v%jmin,1),ldv,Mtmp,MCtmp,w%val(w%jmin,1),Ntmp,NCtmp)
        else
          if( nvecv .eq. 1 ) then
            call dgemm_sb_add_sd_prec_1_1(v%paddedN,v%val(v%jmin,1),Mtmp,MCtmp,w%val(w%jmin,1),Ntmp,NCtmp)
          else if( nvecv .eq. 2 ) then
            call dgemm_sb_add_sd_prec_2_1(v%paddedN,v%val(v%jmin,1),Mtmp,MCtmp,w%val(w%jmin,1),Ntmp,NCtmp)
          else if( nvecv .eq. 4 ) then
            call dgemm_sb_add_sd_prec_4_1(v%paddedN,v%val(v%jmin,1),Mtmp,MCtmp,w%val(w%jmin,1),Ntmp,NCtmp)
          else
            call dgemm_sb_add_sd_prec_k_1(v%paddedN,nvecv,v%val(v%jmin,1),Mtmp,MCtmp,w%val(w%jmin,1),Ntmp,NCtmp)
          end if
        end if
      end if

      iflag = 0
      return
#endif
    end if


    handled = .false.
    allocate(Ntmp(nvecw,nvecw))
    Ntmp = transpose(N%val(N%imin:N%imax,N%jmin:N%jmax))
    if( .not. handled ) then
      allocate(Mtmp(nvecw,nvecv))
      Mtmp = transpose(M%val(M%imin:M%imax,M%jmin:M%jmax))
      ! recognize small block mvecs
      if( .not. strided_w ) then
        if( nvecw .eq. 1 ) then
          if( strided_v ) then
            call dgemm_sB_add_sd_1_strided_k(nrows, nvecv, v%val(v%jmin,1), ldv, Mtmp, w%val,Ntmp)
          else
            call dgemm_sB_add_sd_1_k        (nrows, nvecv, v%val, Mtmp, w%val,Ntmp)
          end if
          handled = .true.
        else if( nvecw .eq. 2 ) then
          if( strided_v ) then
            call dgemm_sB_add_sd_2_strided_k(nrows, nvecv, v%val(v%jmin,1), ldv, Mtmp, w%val,Ntmp)
          else
            call dgemm_sB_add_sd_2_k        (nrows, nvecv, v%val, Mtmp, w%val,Ntmp)
          end if
          handled = .true.
        else if( nvecw .eq. 4 ) then
          if( strided_v ) then
            call dgemm_sB_add_sd_4_strided_k(nrows, nvecv, v%val(v%jmin,1), ldv, Mtmp, w%val,Ntmp)
          else
            call dgemm_sB_add_sd_4_k        (nrows, nvecv, v%val, Mtmp, w%val,Ntmp)
          end if
          handled = .true.
        else if( nvecw .eq. 8 ) then
          if( strided_v ) then
            call dgemm_sB_add_sd_8_strided_k(nrows, nvecv, v%val(v%jmin,1), ldv, Mtmp, w%val,Ntmp)
          else
            call dgemm_sB_add_sd_8_k        (nrows, nvecv, v%val, Mtmp, w%val,Ntmp)
          end if
          handled = .true.
        end if
      end if
      deallocate(Mtmp)
    end if

    if( .not. handled ) then
      allocate(Mtmp(nvecv,nvecw))
      Mtmp = M%val(M%imin:M%imax,M%jmin:M%jmax)
      call dgemm_sB_add_sd_generic(nrows,nvecw,nvecv,v%val(v%jmin,1),ldv, Mtmp, w%val(w%jmin,1),ldw,Ntmp)
      deallocate(Mtmp)
    end if

    iflag = 0
    !--------------------------------------------------------------------------------
  end subroutine mvec_times_sdmat_add_mvec_times_sdmat


  !==================================================================================
  ! special gemm routine for mvec <- mvec*sdMat
  subroutine mvec_times_sdmat_inplace(v,M,iflag)
    !--------------------------------------------------------------------------------
    type(MVec_t),  intent(in)    :: v
    type(SDMat_t), intent(in)    :: M
    integer,       intent(inout) :: iflag
    !--------------------------------------------------------------------------------
    integer :: nrows, nvecv, nvecw, ldv
    logical :: strided_v
    real(kind=8), allocatable :: M_(:,:), MC_(:,:)
    !--------------------------------------------------------------------------------

    ! determine data layout
    nrows = v%map%nlocal(v%map%me)
    nvecv = v%jmax-v%jmin+1
    nvecw = M%jmax-M%jmin+1
    ldv = size(v%val,1)
    if( .not. v%is_view .or. &
      & ( v%jmin .eq. lbound(v%val,1) .and. &
      &   v%jmax .eq. ubound(v%val,1)       ) ) then
      strided_v = .false.
    else
      strided_v = .true.
    end if

    ! copy data to dense block
    allocate(M_(nvecv,nvecw),MC_(nvecv,nvecw))
    M_ = M%val(M%imin:M%imax,M%jmin:M%jmax)
#ifdef PHIST_HIGH_PRECISION_KERNELS
    MC_ = M%err(M%imin:M%imax,M%jmin:M%jmax)
#endif

    if(  FORCE_HIGH_PRECISION .or.iand(iflag,PHIST_ROBUST_REDUCTIONS) .gt. 0 ) then
#ifndef PHIST_HIGH_PRECISION_KERNELS
      iflag = PHIST_NOT_IMPLEMENTED
      return
#else
      ! check if we can do it
      if( strided_v ) then
        iflag = PHIST_NOT_IMPLEMENTED
        return
      end if

      if( nvecv .ne. nvecw ) then
        iflag = PHIST_NOT_IMPLEMENTED
        return
      end if

      if( nvecv .ne. 4 .and. nvecv .ne. 2 .and. nvecv .ne. 1 ) then
        iflag = PHIST_NOT_IMPLEMENTED
        return
      end if

      if( nvecv .eq. 4 ) then
        call dgemm_sB_inplace_prec_4(v%paddedN,v%val(v%jmin,1),M_(1,1),MC_(1,1))
      else if( nvecv .eq. 2 ) then
        call dgemm_sB_inplace_prec_2(v%paddedN,v%val(v%jmin,1),M_(1,1),MC_(1,1))
      else if( nvecv .eq. 1 ) then
        call dgemm_sB_inplace_prec_1(v%paddedN,v%val(v%jmin,1),M_(1,1),MC_(1,1))
      end if

      iflag = 0
      return
#endif
    end if

    ! special implementations
    if( .not. strided_v .and. nvecv .eq. nvecw ) then
      if( nvecv .eq. 8 ) then
        call dgemm_sB_8_8_inplace(nrows,v%val(v%jmin,1),M_(1,1))
        iflag = 0
        return
      else if( nvecv .eq. 4 ) then
        call dgemm_sB_4_4_inplace(nrows,v%val(v%jmin,1),M_(1,1))
        iflag = 0
        return
      else if( nvecv .eq. 2 ) then
        call dgemm_sB_2_2_inplace(nrows,v%val(v%jmin,1),M_(1,1))
        iflag = 0
        return
      else if( nvecv .eq. 1 ) then
        call dgemm_sB_1_1_inplace(nrows,v%val(v%jmin,1),M_(1,1))
        iflag = 0
        return
      end if
    end if

    ! generic case
    call dgemm_sB_generic_inplace(nrows,nvecw,nvecv,v%val(v%jmin,1),ldv,M_(1,1))
    iflag = 0

    !--------------------------------------------------------------------------------
  end subroutine mvec_times_sdmat_inplace


  !==================================================================================
  ! special gemm routine for mvecT_times_mvec
  subroutine mvecT_times_mvec(alpha,v,w,beta,m,iflag)
    use mpi
    !--------------------------------------------------------------------------------
    real(kind=8),  intent(in)    :: alpha
    type(MVec_t),  intent(in)    :: v
    type(MVec_t),  intent(in)    :: w
    real(kind=8),  intent(in)    :: beta
    type(SDMat_t), intent(inout) :: M
    integer,       intent(inout) :: iflag
    !--------------------------------------------------------------------------------
    integer :: nrows, nvecv, nvecw, ldv, ldw
    logical :: strided_v, strided_w
    logical :: handled, tmp_transposed
    real(kind=8), allocatable :: tmp(:,:), tmpC(:,:)
    real(kind=8), allocatable :: tmp_(:,:)
    real(kind=8), allocatable :: localBuff(:,:,:), globalBuff(:,:,:,:), globalBuff_(:,:,:,:)
    real(kind=8), allocatable :: mBuff(:,:), mBuffC(:,:)
    !--------------------------------------------------------------------------------

    ! check if we only need to scale
    if( alpha .eq. 0 ) then
      M%val(M%imin:M%imax,M%jmin:M%jmax) = beta*M%val(M%imin:M%imax,M%jmin:M%jmax)
#ifdef PHIST_HIGH_PRECISION_KERNELS
      M%err(M%imin:M%imax,M%jmin:M%jmax) = beta*M%err(M%imin:M%imax,M%jmin:M%jmax)
#endif
      return
    end if

    ! determine data layout
    nrows = v%map%nlocal(v%map%me)
    nvecv = v%jmax-v%jmin+1
    nvecw = w%jmax-w%jmin+1
    ldv = size(v%val,1)
    ldw = size(w%val,1)
    if( .not. v%is_view .or. &
      & ( v%jmin .eq. lbound(v%val,1) .and. &
      &   v%jmax .eq. ubound(v%val,1)       ) ) then
      strided_v = .false.
    else
      strided_v = .true.
    end if

    if( .not. w%is_view .or. &
      & ( w%jmin .eq. lbound(w%val,1) .and. &
      &   w%jmax .eq. ubound(w%val,1)       ) ) then
      strided_w = .false.
    else
      strided_w = .true.
    end if

    handled = .false.
    tmp_transposed = .false.

    ! check if we need higher precision
    if(  FORCE_HIGH_PRECISION .or.iand(iflag,PHIST_ROBUST_REDUCTIONS) .gt. 0 ) then
#ifndef PHIST_HIGH_PRECISION_KERNELS
      iflag = PHIST_NOT_IMPLEMENTED
      return
#else
      ! check if we can do it
      if( strided_w .or. strided_v ) then
        iflag = PHIST_NOT_IMPLEMENTED
        return
      end if

      if( nvecv .eq. 4 ) then
        handled = .true.
        tmp_transposed = .false.
        allocate(localBuff(4,nvecw,2))
        if( nvecw .eq. 4 ) then
          if( loc(v%val) .eq. loc(w%val) ) then
            call dgemm_sC_self_prec_4(v%paddedN, v%val(v%jmin,1), localBuff(1,1,1), localBuff(1,1,2))
          else
            call dgemm_sC_prec_4_4(v%paddedN, v%val(v%jmin,1), w%val(w%jmin,1), localBuff(1,1,1), localBuff(1,1,2))
          end if
        else if( nvecw .eq. 2 ) then
          call dgemm_sC_prec_4_2(v%paddedN, v%val(v%jmin,1), w%val(w%jmin,1), localBuff(1,1,1), localBuff(1,1,2))
        else if( nvecw .eq. 1 ) then
          call dgemm_sC_prec_4_1(v%paddedN, v%val(v%jmin,1), w%val(w%jmin,1), localBuff(1,1,1), localBuff(1,1,2))
        else
          call dgemm_sC_prec_4_k(v%paddedN, nvecw, v%val(v%jmin,1), w%val(w%jmin,1), localBuff(1,1,1), localBuff(1,1,2))
        end if
      end if
      if( .not. handled .and. nvecw .eq. 4 ) then
        handled = .true.
        tmp_transposed = .true.
        allocate(localBuff(4,nvecv,2))
        if( nvecv .eq. 2 ) then
          call dgemm_sC_prec_4_2(v%paddedN, w%val(w%jmin,1), v%val(v%jmin,1), localBuff(1,1,1), localBuff(1,1,2))
        else if( nvecv .eq. 1 ) then
          call dgemm_sC_prec_4_1(v%paddedN, w%val(w%jmin,1), v%val(v%jmin,1), localBuff(1,1,1), localBuff(1,1,2))
        else
          call dgemm_sC_prec_4_k(v%paddedN, nvecv, w%val(w%jmin,1), v%val(v%jmin,1), localBuff(1,1,1), localBuff(1,1,2))
        end if
      end if
      if( .not. handled .and. nvecv .eq. 2 ) then
        handled = .true.
        tmp_transposed = .false.
        allocate(localBuff(2,nvecw,2))
        if(nvecw .eq. 2 ) then
          if( loc(v%val) .eq. loc(w%val) ) then
            call dgemm_sC_self_prec_2(v%paddedN, v%val(v%jmin,1), localBuff(1,1,1), localBuff(1,1,2))
          else
            call dgemm_sC_prec_2_2(v%paddedN, v%val(v%jmin,1), w%val(w%jmin,1), localBuff(1,1,1), localBuff(1,1,2))
          end if
        else if( nvecw .eq. 1 ) then
          call dgemm_sC_prec_2_1(v%paddedN, v%val(v%jmin,1), w%val(w%jmin,1), localBuff(1,1,1), localBuff(1,1,2))
        else
          call dgemm_sC_prec_2_k(v%paddedN, nvecw, v%val(v%jmin,1), w%val(w%jmin,1), localBuff(1,1,1), localBuff(1,1,2))
        end if
      end if
      if( .not. handled .and. nvecw .eq. 2 ) then
        handled = .true.
        tmp_transposed = .true.
        allocate(localBuff(2,nvecv,2))
        if( nvecv .eq. 1 ) then
          call dgemm_sC_prec_2_1(v%paddedN, w%val(w%jmin,1), v%val(v%jmin,1), localBuff(1,1,1), localBuff(1,1,2))
        else
          call dgemm_sC_prec_2_k(v%paddedN, nvecv, w%val(w%jmin,1), v%val(v%jmin,1), localBuff(1,1,1), localBuff(1,1,2))
        end if
      end if
      if( .not. handled .and. nvecv .eq. 1 ) then
        handled = .true.
        tmp_transposed = .false.
        allocate(localBuff(1,nvecw,2))
        if(nvecw .eq. 1 ) then
          if( loc(v%val) .eq. loc(w%val) ) then
            call ddot_self_prec_1(v%paddedN, v%val(v%jmin,1), localBuff(1,1,1), localBuff(1,1,2))
          else
            call ddot_prec_1(v%paddedN, v%val(v%jmin,1), w%val(w%jmin,1), localBuff(1,1,1), localBuff(1,1,2))
          end if
        else
          call dgemm_sC_prec_1_k(v%paddedN, nvecw, v%val(v%jmin,1), w%val(w%jmin,1), localBuff(1,1,1), localBuff(1,1,2))
        end if
      end if
      if( .not. handled .and. nvecw .eq. 1 ) then
        handled = .true.
        tmp_transposed = .true.
        allocate(localBuff(1,nvecv,2))
        call dgemm_sC_prec_1_k(v%paddedN, nvecv, w%val(w%jmin,1), v%val(v%jmin,1), localBuff(1,1,1), localBuff(1,1,2))
      end if

      if( .not. handled ) then
        iflag = PHIST_NOT_IMPLEMENTED
        return
      end if
!write(*,*) 'here v', nvecv, v%val
!write(*,*) 'here w', nvecw, w%val
!write(*,*) 'here buff', localBuff(:,:,1)
!write(*,*) 'here buffC', localBuff(:,:,2)

      ! gather results
      if( tmp_transposed ) then
        allocate(globalBuff(nvecw,nvecv,2,v%map%nProcs))
        allocate(globalBuff_(nvecw,nvecv,v%map%nProcs,2))
      else
        allocate(globalBuff(nvecv,nvecw,2,v%map%nProcs))
        allocate(globalBuff_(nvecv,nvecw,v%map%nProcs,2))
      end if
      call MPI_Allgather(localBuff,2*nvecv*nvecw,MPI_DOUBLE_PRECISION, &
        &                globalBuff,2*nvecv*nvecw,MPI_DOUBLE_PRECISION, &
        &                v%map%comm, iflag)
      globalBuff_(:,:,:,1) = globalBuff(:,:,1,:)
      globalBuff_(:,:,:,2) = globalBuff(:,:,2,:)
!write(*,*) 'here globalBuff', globalBuff_(:,:,:,1)
!write(*,*) 'here globalBuffC', globalBuff_(:,:,:,2)
      ! MPI reduction
      if( nvecv*nvecw .eq. 1 ) then
        call prec_reduction_1(v%map%nProcs, globalBuff_(1,1,1,1), globalBuff_(1,1,1,2), &
          &                                 localBuff(1,1,1), localBuff(1,1,2)          )
      else if( nvecv*nvecw .eq. 2 ) then
        call prec_reduction_2(v%map%nProcs, globalBuff_(1,1,1,1), globalBuff_(1,1,1,2), &
          &                                 localBuff(1,1,1), localBuff(1,1,2)          )
      else if( nvecv*nvecw .eq. 4 ) then
        call prec_reduction_4(v%map%nProcs, globalBuff_(1,1,1,1), globalBuff_(1,1,1,2), &
          &                                 localBuff(1,1,1), localBuff(1,1,2)          )
      else if( mod(nvecv*nvecw,4) .eq. 0 ) then
        call prec_reduction_4k(v%map%nProcs, nvecv*nvecw, globalBuff_(1,1,1,1), globalBuff_(1,1,1,2), &
          &                                               localBuff(1,1,1), localBuff(1,1,2)          )
      else if( mod(nvecv*nvecw,2) .eq. 0 ) then
        call prec_reduction_2k(v%map%nProcs, nvecv*nvecw, globalBuff_(1,1,1,1), globalBuff_(1,1,1,2), &
          &                                               localBuff(1,1,1), localBuff(1,1,2)          )
      else
        call prec_reduction_k(v%map%nProcs, nvecv*nvecw, globalBuff_(1,1,1,1), globalBuff_(1,1,1,2), &
          &                                              localBuff(1,1,1), localBuff(1,1,2)          )
      end if
!write(*,*) 'here vTw', localBuff(:,:,1)
!write(*,*) 'here vTw _C', localBuff(:,:,2)

      ! copy it to a buffer again for the final calculation of m = alpha*res + beta*m
      allocate(mBuff(nvecv,nvecw), tmp(nvecv,nvecw), mBuffC(nvecv,nvecw), tmpC(nvecv,nvecw))
      mBuff = M%val(M%imin:M%imax,M%jmin:M%jmax)
      mBuffC = M%err(M%imin:M%imax,M%jmin:M%jmax)
      if( tmp_transposed ) then
        tmp = transpose(localBuff(:,:,1))
        tmpC = transpose(localBuff(:,:,2))
      else
        tmp = localBuff(:,:,1)
        tmpC = localBuff(:,:,2)
      end if
!write(*,*) 'here mBuff', mBuff
!write(*,*) 'here tmp', tmp
!write(*,*) 'here alpha beta', alpha, beta
      call daxpby_prec(nvecv*nvecw,alpha,tmp(1,1),tmpC(1,1),beta,mBuff(1,1),mBuffC(1,1))
!write(*,*) 'here mBuff', mBuff
!write(*,*) 'here mBuffC', mBuffC
      ! set result
      M%val(M%imin:M%imax,M%jmin:M%jmax) = mBuff
      M%err(M%imin:M%imax,M%jmin:M%jmax) = mBuffC
      return
#endif
    end if


    if( .not. strided_v ) then
      if( nvecv .eq. 1 ) then
        allocate(tmp(nvecv,nvecw))
        if( strided_w ) then
          call dgemm_sC_strided_1(nrows,nvecw,v%val,w%val(w%jmin,1),ldw,tmp)
        else if( nvecw .eq. 1 ) then
          call dgemm_sC_1_1(nrows,v%val,w%val,tmp)
        else
          call dgemm_sC_1(nrows,nvecw,v%val,w%val,tmp)
        end if
        handled = .true.
      else if( nvecv .eq. 2 ) then
        allocate(tmp(nvecv,nvecw))
        if( strided_w ) then
          call dgemm_sC_strided_2(nrows,nvecw,v%val,w%val(w%jmin,1),ldw,tmp)
        else if( nvecw .eq. 2 ) then
          call dgemm_sC_2_2(nrows,v%val,w%val,tmp)
        else
          call dgemm_sC_2(nrows,nvecw,v%val,w%val,tmp)
        end if
        handled = .true.
      else if( nvecv .eq. 4 ) then
        allocate(tmp(nvecv,nvecw))
        if( strided_w ) then
          call dgemm_sC_strided_4(nrows,nvecw,v%val,w%val(w%jmin,1),ldw,tmp)
        else if( nvecw .eq. 4 ) then
          call dgemm_sC_4_4(nrows,v%val,w%val,tmp)
        else
          call dgemm_sC_4(nrows,nvecw,v%val,w%val,tmp)
        end if
        handled = .true.
      else if( nvecv .eq. 8 ) then
        allocate(tmp(nvecv,nvecw))
        if( strided_w ) then
          call dgemm_sC_strided_8(nrows,nvecw,v%val,w%val(w%jmin,1),ldw,tmp)
        else if( nvecw .eq. 8 ) then
          call dgemm_sC_8_8(nrows,v%val,w%val,tmp)
        else
          call dgemm_sC_8(nrows,nvecw,v%val,w%val,tmp)
        end if
        handled = .true.
      end if
    end if

    if( .not. handled .and. .not. strided_w ) then
      if( nvecw .eq. 1 ) then
        allocate(tmp(nvecw,nvecv))
        if( strided_v ) then
          call dgemm_sC_strided_1(nrows,nvecv,w%val,v%val(v%jmin,1),ldv,tmp)
        else
          call dgemm_sC_1(nrows,nvecv,w%val,v%val,tmp)
        end if
        handled = .true.
        tmp_transposed = .true.
      else if( nvecw .eq. 2 ) then
        allocate(tmp(nvecw,nvecv))
        if( strided_v ) then
          call dgemm_sC_strided_2(nrows,nvecv,w%val,v%val(v%jmin,1),ldv,tmp)
        else
          call dgemm_sC_2(nrows,nvecv,w%val,v%val,tmp)
        end if
        handled = .true.
        tmp_transposed = .true.
      else if( nvecw .eq. 4 ) then
        allocate(tmp(nvecw,nvecv))
        if( strided_v ) then
          call dgemm_sC_strided_4(nrows,nvecv,w%val,v%val(v%jmin,1),ldv,tmp)
        else
          call dgemm_sC_4(nrows,nvecv,w%val,v%val,tmp)
        end if
        handled = .true.
        tmp_transposed = .true.
      else if( nvecw .eq. 8 ) then
        allocate(tmp(nvecw,nvecv))
        if( strided_v ) then
          call dgemm_sC_strided_8(nrows,nvecv,w%val,v%val(v%jmin,1),ldv,tmp)
        else
          call dgemm_sC_8(nrows,nvecv,w%val,v%val,tmp)
        end if
        handled = .true.
        tmp_transposed = .true.
      end if
    end if


    if( .not. handled ) then
      ! generic case
      allocate(tmp(nvecv,nvecw))

      call dgemm_sC_generic(nrows,nvecv,nvecw,v%val(v%jmin,1),ldv,w%val(w%jmin,1),ldw,tmp)

    end if

    allocate(tmp_(nvecv,nvecw))
    if( tmp_transposed ) then
      tmp_ = transpose(tmp)
      deallocate(tmp)
      allocate(tmp(nvecv,nvecw))
    else
      tmp_ = tmp
    end if
    call MPI_Allreduce(tmp_,tmp,nvecv*nvecw,MPI_DOUBLE_PRECISION,MPI_SUM,v%map%comm,iflag)

    M%val(M%imin:M%imax,M%jmin:M%jmax) = alpha*tmp+beta*M%val(M%imin:M%imax,M%jmin:M%jmax)
#ifdef PHIST_HIGH_PRECISION_KERNELS
    M%err(M%imin:M%imax,M%jmin:M%jmax) = 0.
#endif
    !--------------------------------------------------------------------------------
  end subroutine mvecT_times_mvec


  !==================================================================================
  ! special augmented double-gemm routine for mvecT_times_mvec_times_sdMat_inplace
  subroutine mvecT_times_mvec_times_sdMat_inplace(alpha,v,w,n,beta,m,iflag)
    use mpi
    !--------------------------------------------------------------------------------
    real(kind=8),  intent(in)    :: alpha
    type(MVec_t),  intent(in)    :: v
    type(MVec_t),  intent(inout) :: w
    type(SDMat_t), intent(in)    :: N
    real(kind=8),  intent(in)    :: beta
    type(SDMat_t), intent(inout) :: M
    integer,       intent(inout) :: iflag
    !--------------------------------------------------------------------------------
    integer :: nrows, nvecv, nvecw, ldv, ldw
    logical :: strided_v, strided_w
    logical :: handled, tmp_transposed
    real(kind=8), allocatable :: tmp(:,:), tmpC(:,:)
    real(kind=8), allocatable :: tmpN(:,:), tmpNC(:,:)
    real(kind=8), allocatable :: tmp_(:,:)
    real(kind=8), allocatable :: localBuff(:,:,:), globalBuff(:,:,:,:), globalBuff_(:,:,:,:)
    real(kind=8), allocatable :: mBuff(:,:), mBuffC(:,:)
    !--------------------------------------------------------------------------------

    ! check if we only need to scale
    if( alpha .eq. 0 ) then
      M%val(M%imin:M%imax,M%jmin:M%jmax) = beta*M%val(M%imin:M%imax,M%jmin:M%jmax)
#ifdef PHIST_HIGH_PRECISION_KERNELS
      M%err(M%imin:M%imax,M%jmin:M%jmax) = beta*M%err(M%imin:M%imax,M%jmin:M%jmax)
#endif
      call mvec_times_sdMat_inplace(w,n,iflag)
      return
    end if

    ! determine data layout
    nrows = v%map%nlocal(v%map%me)
    nvecv = v%jmax-v%jmin+1
    nvecw = w%jmax-w%jmin+1
    ldv = size(v%val,1)
    ldw = size(w%val,1)
    if( .not. v%is_view .or. &
      & ( v%jmin .eq. lbound(v%val,1) .and. &
      &   v%jmax .eq. ubound(v%val,1)       ) ) then
      strided_v = .false.
    else
      strided_v = .true.
    end if

    if( .not. w%is_view .or. &
      & ( w%jmin .eq. lbound(w%val,1) .and. &
      &   w%jmax .eq. ubound(w%val,1)       ) ) then
      strided_w = .false.
    else
      strided_w = .true.
    end if

    handled = .false.
    tmp_transposed = .false.

    ! allocate dense tmp mem for n
    allocate(tmpN (nvecw,nvecw)) ; tmpN  = N%val(N%imin:N%imax,N%jmin:N%jmax)
#ifdef PHIST_HIGH_PRECISION_KERNELS
    allocate(tmpNC(nvecw,nvecw)) ; tmpNC = N%err(N%imin:N%imax,N%jmin:N%jmax)
#endif

    ! check if we need higher precision
    if(  FORCE_HIGH_PRECISION .or.iand(iflag,PHIST_ROBUST_REDUCTIONS) .gt. 0 ) then
#ifndef PHIST_HIGH_PRECISION_KERNELS
      iflag = PHIST_NOT_IMPLEMENTED
      return
#else
      ! check if we can do it
      if( strided_w .or. strided_v ) then
        iflag = PHIST_NOT_IMPLEMENTED
        return
      end if

      if( nvecv .eq. 4 ) then
        handled = .true.
        tmp_transposed = .true.
        allocate(localBuff(nvecw,4,2))
        if( nvecw .eq. 4 ) then
          if( loc(v%val) .eq. loc(w%val) ) then
            call dgemm_fused_sCD_self_prec_4(v%paddedN, v%val(v%jmin,1), tmpN(1,1), tmpNC(1,1), localBuff(1,1,1), localBuff(1,1,2))
          else
            call dgemm_fused_sCD_prec_4_4(v%paddedN, v%val(v%jmin,1), w%val(w%jmin,1), tmpN(1,1), tmpNC(1,1), localBuff(1,1,1), localBuff(1,1,2))
          end if
        else if( nvecw .eq. 2 ) then
          call dgemm_fused_sCD_prec_4_2(v%paddedN, v%val(v%jmin,1), w%val(w%jmin,1), tmpN(1,1), tmpNC(1,1), localBuff(1,1,1), localBuff(1,1,2))
        else if( nvecw .eq. 1 ) then
          call dgemm_fused_sCD_prec_4_1(v%paddedN, v%val(v%jmin,1), w%val(w%jmin,1), tmpN(1,1), tmpNC(1,1), localBuff(1,1,1), localBuff(1,1,2))
        else
          handled = .false.
        end if
      end if
      if( .not. handled .and. nvecw .eq. 4 ) then
        handled = .true.
        tmp_transposed = .true.
        allocate(localBuff(4,nvecv,2))
        if( nvecv .eq. 2 ) then
          call dgemm_fused_sCD_prec_2_4(v%paddedN, v%val(v%jmin,1), w%val(w%jmin,1), tmpN(1,1), tmpNC(1,1), localBuff(1,1,1), localBuff(1,1,2))
        else if( nvecv .eq. 1 ) then
          call dgemm_fused_sCD_prec_1_4(v%paddedN, v%val(v%jmin,1), w%val(w%jmin,1), tmpN(1,1), tmpNC(1,1), localBuff(1,1,1), localBuff(1,1,2))
        else
          call dgemm_fused_sCD_prec_k_4(v%paddedN, nvecv, v%val(v%jmin,1), w%val(w%jmin,1), tmpN(1,1), tmpNC(1,1), localBuff(1,1,1), localBuff(1,1,2))
        end if
      end if
      if( .not. handled .and. nvecv .eq. 2 ) then
        handled = .true.
        tmp_transposed = .true.
        allocate(localBuff(nvecw,2,2))
        if(nvecw .eq. 2 ) then
          if( loc(v%val) .eq. loc(w%val) ) then
            call dgemm_fused_sCD_self_prec_2(v%paddedN, v%val(v%jmin,1), tmpN(1,1), tmpNC(1,1), localBuff(1,1,1), localBuff(1,1,2))
          else
            call dgemm_fused_sCD_prec_2_2(v%paddedN, v%val(v%jmin,1), w%val(w%jmin,1), tmpN(1,1), tmpNC(1,1), localBuff(1,1,1), localBuff(1,1,2))
          end if
        else if( nvecw .eq. 1 ) then
          call dgemm_fused_sCD_prec_2_1(v%paddedN, v%val(v%jmin,1), w%val(w%jmin,1), tmpN(1,1), tmpNC(1,1), localBuff(1,1,1), localBuff(1,1,2))
        else
          handled = .false.
        end if
      end if
      if( .not. handled .and. nvecw .eq. 2 ) then
        handled = .true.
        tmp_transposed = .true.
        allocate(localBuff(2,nvecv,2))
        if( nvecv .eq. 1 ) then
          call dgemm_fused_sCD_prec_1_2(v%paddedN, v%val(v%jmin,1), w%val(w%jmin,1), tmpN(1,1), tmpNC(1,1), localBuff(1,1,1), localBuff(1,1,2))
        else
          call dgemm_fused_sCD_prec_k_2(v%paddedN, nvecv, v%val(v%jmin,1), w%val(w%jmin,1), tmpN(1,1), tmpNC(1,1), localBuff(1,1,1), localBuff(1,1,2))
        end if
      end if
      if( .not. handled .and. nvecv .eq. 1 ) then
        handled = .true.
        tmp_transposed = .true.
        allocate(localBuff(nvecw,1,2))
        if(nvecw .eq. 1 ) then
          if( loc(v%val) .eq. loc(w%val) ) then
            call ddot_fused_scale_self_prec_1(v%paddedN, v%val(v%jmin,1), tmpN(1,1), tmpNC(1,1), localBuff(1,1,1), localBuff(1,1,2))
          else
            call ddot_fused_scale_prec_1(v%paddedN, v%val(v%jmin,1), w%val(w%jmin,1), tmpN(1,1), tmpNC(1,1), localBuff(1,1,1), localBuff(1,1,2))
          end if
        else
          handled = .false.
        end if
      end if
      if( .not. handled .and. nvecw .eq. 1 ) then
        handled = .true.
        tmp_transposed = .true.
        allocate(localBuff(1,nvecv,2))
        call dgemm_fused_sCD_prec_k_1(v%paddedN, nvecv, v%val(v%jmin,1), w%val(w%jmin,1), tmpN(1,1), tmpNC(1,1), localBuff(1,1,1), localBuff(1,1,2))
      end if

      if( .not. handled ) then
        iflag = PHIST_NOT_IMPLEMENTED
        return
      end if
!write(*,*) 'here v', nvecv, v%val
!write(*,*) 'here w', nvecw, w%val
!write(*,*) 'here buff', localBuff(:,:,1)
!write(*,*) 'here buffC', localBuff(:,:,2)

      ! gather results
      if( tmp_transposed ) then
        allocate(globalBuff(nvecw,nvecv,2,v%map%nProcs))
        allocate(globalBuff_(nvecw,nvecv,v%map%nProcs,2))
      else
        allocate(globalBuff(nvecv,nvecw,2,v%map%nProcs))
        allocate(globalBuff_(nvecv,nvecw,v%map%nProcs,2))
      end if
      call MPI_Allgather(localBuff,2*nvecv*nvecw,MPI_DOUBLE_PRECISION, &
        &                globalBuff,2*nvecv*nvecw,MPI_DOUBLE_PRECISION, &
        &                v%map%comm, iflag)
      globalBuff_(:,:,:,1) = globalBuff(:,:,1,:)
      globalBuff_(:,:,:,2) = globalBuff(:,:,2,:)
!write(*,*) 'here globalBuff', globalBuff_(:,:,:,1)
!write(*,*) 'here globalBuffC', globalBuff_(:,:,:,2)
      ! MPI reduction
      if( nvecv*nvecw .eq. 1 ) then
        call prec_reduction_1(v%map%nProcs, globalBuff_(1,1,1,1), globalBuff_(1,1,1,2), &
          &                                 localBuff(1,1,1), localBuff(1,1,2)          )
      else if( nvecv*nvecw .eq. 2 ) then
        call prec_reduction_2(v%map%nProcs, globalBuff_(1,1,1,1), globalBuff_(1,1,1,2), &
          &                                 localBuff(1,1,1), localBuff(1,1,2)          )
      else if( nvecv*nvecw .eq. 4 ) then
        call prec_reduction_4(v%map%nProcs, globalBuff_(1,1,1,1), globalBuff_(1,1,1,2), &
          &                                 localBuff(1,1,1), localBuff(1,1,2)          )
      else if( mod(nvecv*nvecw,4) .eq. 0 ) then
        call prec_reduction_4k(v%map%nProcs, nvecv*nvecw, globalBuff_(1,1,1,1), globalBuff_(1,1,1,2), &
          &                                               localBuff(1,1,1), localBuff(1,1,2)          )
      else if( mod(nvecv*nvecw,2) .eq. 0 ) then
        call prec_reduction_2k(v%map%nProcs, nvecv*nvecw, globalBuff_(1,1,1,1), globalBuff_(1,1,1,2), &
          &                                               localBuff(1,1,1), localBuff(1,1,2)          )
      else
        call prec_reduction_k(v%map%nProcs, nvecv*nvecw, globalBuff_(1,1,1,1), globalBuff_(1,1,1,2), &
          &                                              localBuff(1,1,1), localBuff(1,1,2)          )
      end if
!write(*,*) 'here vTw', localBuff(:,:,1)
!write(*,*) 'here vTw _C', localBuff(:,:,2)

      ! copy it to a buffer again for the final calculation of m = alpha*res + beta*m
      allocate(mBuff(nvecv,nvecw), tmp(nvecv,nvecw), mBuffC(nvecv,nvecw), tmpC(nvecv,nvecw))
      mBuff = M%val(M%imin:M%imax,M%jmin:M%jmax)
      mBuffC = M%err(M%imin:M%imax,M%jmin:M%jmax)
      if( tmp_transposed ) then
        tmp = transpose(localBuff(:,:,1))
        tmpC = transpose(localBuff(:,:,2))
      else
        tmp = localBuff(:,:,1)
        tmpC = localBuff(:,:,2)
      end if
!write(*,*) 'here mBuff', mBuff
!write(*,*) 'here tmp', tmp
!write(*,*) 'here alpha beta', alpha, beta
      call daxpby_prec(nvecv*nvecw,alpha,tmp(1,1),tmpC(1,1),beta,mBuff(1,1),mBuffC(1,1))
!write(*,*) 'here mBuff', mBuff
!write(*,*) 'here mBuffC', mBuffC
      ! set result
      M%val(M%imin:M%imax,M%jmin:M%jmax) = mBuff
      M%err(M%imin:M%imax,M%jmin:M%jmax) = mBuffC
      return
#endif
    end if


    if( .not. strided_v ) then
      if( nvecv .eq. 1 ) then
        allocate(tmp(nvecv,nvecw))
        if( strided_w ) then
          call dgemm_fused_sCD_1_strided_k(nrows,nvecw,v%val,w%val(w%jmin,1),ldw,tmpN,tmp)
        else
          if( loc(v%val) .eq. loc(w%val) ) then
            call dgemm_fused_sCD_1_self(nrows,v%val,tmpN,tmp)
          else
            call dgemm_fused_sCD_1_k(nrows,nvecw,v%val,w%val,tmpN,tmp)
          end if
        end if
        handled = .true.
      else if( nvecv .eq. 2 ) then
        allocate(tmp(nvecv,nvecw))
        if( strided_w ) then
          call dgemm_fused_sCD_2_strided_k(nrows,nvecw,v%val,w%val(w%jmin,1),ldw,tmpN,tmp)
        else
          if( loc(v%val) .eq. loc(w%val) ) then
            call dgemm_fused_sCD_2_self(nrows,v%val,tmpN,tmp)
          else
            call dgemm_fused_sCD_2_k(nrows,nvecw,v%val,w%val,tmpN,tmp)
          end if
        end if
        handled = .true.
      else if( nvecv .eq. 4 ) then
        allocate(tmp(nvecv,nvecw))
        if( strided_w ) then
          call dgemm_fused_sCD_4_strided_k(nrows,nvecw,v%val,w%val(w%jmin,1),ldw,tmpN,tmp)
        else
          if( loc(v%val) .eq. loc(w%val) ) then
            call dgemm_fused_sCD_4_self(nrows,v%val,tmpN,tmp)
          else
            call dgemm_fused_sCD_4_k(nrows,nvecw,v%val,w%val,tmpN,tmp)
          end if
        end if
        handled = .true.
      else if( nvecv .eq. 8 ) then
        allocate(tmp(nvecv,nvecw))
        if( strided_w ) then
          call dgemm_fused_sCD_8_strided_k(nrows,nvecw,v%val,w%val(w%jmin,1),ldw,tmpN,tmp)
        else
          if( loc(v%val) .eq. loc(w%val) ) then
            call dgemm_fused_sCD_8_self(nrows,v%val,tmpN,tmp)
          else
            call dgemm_fused_sCD_8_k(nrows,nvecw,v%val,w%val,tmpN,tmp)
          end if
        end if
        handled = .true.
      end if
    end if

    if( .not. handled .and. .not. strided_w ) then
      if( nvecw .eq. 1 ) then
        allocate(tmp(nvecv,nvecw))
        if( strided_v ) then
          call dgemm_fused_sCD_strided_k_1(nrows,nvecv,v%val(v%jmin,1),ldv,w%val,tmpN,tmp)
        else
          call dgemm_fused_sCD_k_1(nrows,nvecv,v%val,w%val,tmpN,tmp)
        end if
        handled = .true.
      else if( nvecw .eq. 2 ) then
        allocate(tmp(nvecv,nvecw))
        if( strided_v ) then
          call dgemm_fused_sCD_strided_k_2(nrows,nvecv,v%val(v%jmin,1),ldv,w%val,tmpN,tmp)
        else
          call dgemm_fused_sCD_k_2(nrows,nvecv,v%val,w%val,tmpN,tmp)
        end if
        handled = .true.
      else if( nvecw .eq. 4 ) then
        allocate(tmp(nvecv,nvecw))
        if( strided_v ) then
          call dgemm_fused_sCD_strided_k_4(nrows,nvecv,v%val(v%jmin,1),ldv,w%val,tmpN,tmp)
        else
          call dgemm_fused_sCD_k_4(nrows,nvecv,v%val,w%val,tmpN,tmp)
        end if
        handled = .true.
      else if( nvecw .eq. 8 ) then
        allocate(tmp(nvecv,nvecw))
        if( strided_v ) then
          call dgemm_fused_sCD_strided_k_8(nrows,nvecv,v%val(v%jmin,1),ldv,w%val,tmpN,tmp)
        else
          call dgemm_fused_sCD_k_8(nrows,nvecv,v%val,w%val,tmpN,tmp)
        end if
        handled = .true.
      end if
    end if


    if( .not. handled ) then
      ! generic case
      allocate(tmp(nvecv,nvecw))

      if( nvecv .eq. nvecw .and. loc(v%val(v%jmin,1)) .eq. loc(w%val(w%jmin,1)) ) then
        call dgemm_fused_sCD_generic_self(nrows,nvecw,w%val(w%jmin,1),ldw,tmpN,tmp)
      else
        call dgemm_fused_sCD_generic(nrows,nvecv,nvecw,v%val(v%jmin,1),ldv,w%val(w%jmin,1),ldw,tmpN,tmp)
      end if

    end if

    allocate(tmp_(nvecv,nvecw))
    if( tmp_transposed ) then
      tmp_ = transpose(tmp)
      deallocate(tmp)
      allocate(tmp(nvecv,nvecw))
    else
      tmp_ = tmp
    end if
    call MPI_Allreduce(tmp_,tmp,nvecv*nvecw,MPI_DOUBLE_PRECISION,MPI_SUM,v%map%comm,iflag)

    M%val(M%imin:M%imax,M%jmin:M%jmax) = alpha*tmp+beta*M%val(M%imin:M%imax,M%jmin:M%jmax)
#ifdef PHIST_HIGH_PRECISION_KERNELS
    M%err(M%imin:M%imax,M%jmin:M%jmax) = 0.
#endif
    !--------------------------------------------------------------------------------
  end subroutine mvecT_times_mvec_times_sdMat_inplace


  !==================================================================================
  ! orthogonalize v with v = QR, fill with orthogonal random vectors if not full rank
  subroutine mvec_QR(v,R,nullSpaceDim)
    !--------------------------------------------------------------------------------
    type(MVec_t),   intent(inout) :: v
    type(SDMat_t),  intent(inout) :: R
    integer,        intent(out)   :: nullSpaceDim
    !--------------------------------------------------------------------------------
    real(kind=8), parameter :: eps = 1.e-10
    integer :: i, i_, nvec, rank, k, iflag
    type(MVec_t) :: vi, vipn
    type(SDMat_t) :: Ripn
    real(kind=8) :: rii(1:1)
    integer :: idum
    !--------------------------------------------------------------------------------
#if defined(TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
    write(*,*) 'entering mvec_QR'
#endif

    nvec = v%jmax-v%jmin+1
    idum = 0

    ! setup views vi, vipn, Ripn
    vipn%jmax = v%jmax
    vipn%is_view = .true.
    vi%is_view = .true.
    vi%val=>v%val
    vipn%val=>v%val
    vi%map = v%map
    vipn%map = v%map
    Ripn%jmax = R%jmax
    Ripn%is_view = .true.
    Ripn%val => R%val
#ifdef PHIST_HIGH_PRECISION_KERNELS
    Ripn%err => R%err
#endif
    Ripn%comm = R%comm

    ! simple modified gram schmidt... probably SLOW
    ! R = 0
    nullSpaceDim = 0
    R%val(R%imin:R%imax,R%jmin:R%jmax) = 0.
    i_ = 0
    do i = 0, nvec-1, 1
      ! create view of column i
      vi%jmin = v%jmin+i
      vi%jmax = v%jmin+i
      call mvec_norm2(vi,rii,idum)
      R%val(R%imin+i_,R%jmin+i) = rii(1)
      if( rii(1) .lt. eps ) then
        nullSpaceDim = nullSpaceDim + 1
        rii(1) = 0
        R%val(R%imin+i_,R%jmin+i) = 0
      else
        rii(1) = 1._8/rii(1)
        call mvec_scale(vi,rii(1))
        if( i .lt. nvec-1 ) then
          vipn%jmin = v%jmin+i+1
          Ripn%imin = R%imin+i_
          Ripn%imax = R%imin+i_
          Ripn%jmin = R%jmin+i+1
          call mvecT_times_mvec(1._8,vi,vipn,0._8,Ripn,idum)
          call mvec_times_sdmat(-1._8,vi,Ripn,1._8,vipn,idum)
        end if
        ! copy vi to column i_, because previous vector was not linearly independent
        if( i .ne. i_ ) then
!$omp parallel do schedule(static)
          do k = 1, v%map%nlocal(v%map%me), 1
            v%val(v%jmin+i_,k) = v%val(v%jmin+i,k)
          end do
        end if
        i_ = i_ + 1
      end if
    end do

    ! try to generate orthogonal random vectors
    if( nullSpaceDim .gt. 0 ) then
      rank = nvec - nullSpaceDim
      ! reuse vipn for view
      vipn%jmin = v%jmin+rank
      call mvec_random(vipn)

      ! reuse Ripn for temporary storage
      Ripn%val=>null()
      Ripn%is_view = .false.
      Ripn%imin = 1
      Ripn%imax = max(1,rank)
      Ripn%jmin = 1
      Ripn%jmax = nullSpaceDim
      allocate(Ripn%val(Ripn%imax,Ripn%jmax))
      Ripn%val = 0._8
#ifdef PHIST_HIGH_PRECISION_KERNELS
      allocate(Ripn%err(Ripn%imax,Ripn%jmax))
      Ripn%err = 0._8
#endif

      if( rank .gt. 0 ) then
        ! orthog. random vectors wrt. previous vectors
        vi%jmin = v%jmin
        vi%jmax = v%jmin+rank-1
        vipn%jmin = v%jmin+rank
        call mvecT_times_mvec(1._8,vi,vipn,0._8,Ripn,idum)
        call mvec_times_sdmat(-1._8,vi,Ripn,1._8,vipn,idum)
      end if

      ! orthog. random vectors wrt. each other
      Ripn%imin = 1
      Ripn%imax = 1
      do i = rank, nvec-1, 1
        ! create view of column i
        vi%jmin = v%jmin+i
        vi%jmax = v%jmin+i
        call mvec_norm2(vi,rii,idum)
        if( rii(1) .lt. eps ) then
          write(*,*) 'error during orthogonalization'
          flush(6)
          nullSpaceDim = -1
          deallocate(Ripn%val)
#ifdef PHIST_HIGH_PRECISION_KERNELS
          deallocate(Ripn%err)
#endif
          return
        end if
        rii(1) = 1._8/rii(1)
        call mvec_scale(vi,rii(1))
        if( i .lt. nvec-1 ) then
          vipn%jmin = v%jmin+i+1
          Ripn%jmin = i-rank+2
          call mvecT_times_mvec(1._8,vi,vipn,0._8,Ripn,idum)
          call mvec_times_sdmat(-1._8,vi,Ripn,1._8,vipn,idum)
        end if
      end do

      deallocate(Ripn%val)
#ifdef PHIST_HIGH_PRECISION_KERNELS
      deallocate(Ripn%err)
#endif
    end if

    !--------------------------------------------------------------------------------
  end subroutine mvec_QR
 

  !==================================================================================
  ! wrapper routines 

  subroutine phist_Dmvec_create(mvec_ptr, map_ptr, nvec, ierr) &
    & bind(C,name='phist_Dmvec_create_f') ! circumvent bug in opari (openmp instrumentalization)
    use, intrinsic :: iso_c_binding
    use :: omp_lib
    !--------------------------------------------------------------------------------
    type(C_PTR),        intent(out) :: mvec_ptr
    type(C_PTR),        value       :: map_ptr
    integer(C_INT32_T), value       :: nvec
    integer(C_INT),     intent(out) :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: mvec
    type(Map_t), pointer :: map
    integer :: i, nt, padding
    type(C_PTR) :: rawMem
#if defined(TESTING)
    integer(C_INTPTR_T) :: dummy
#endif
    !--------------------------------------------------------------------------------

    if( nvec .le. 0 ) then
      ierr=-44
      return
    end if

    call c_f_pointer(map_ptr, map)
    allocate(mvec)
    mvec%is_view = .false.
    mvec%jmin = 1
    mvec%jmax = nvec
    mvec%map = map
    nt = omp_get_max_threads()
    padding = nt*4
    mvec%paddedN = (map%nlocal(map%me)/padding+1)*padding
#ifdef F_DEBUG
    write(*,*) 'creating new mvec with dimensions:', nvec, map%nlocal(map%me), 'address', transfer(c_loc(mvec),dummy)
    flush(6)
#endif
    ! allocate memory
    call allocate_aligned( int((/nvec,mvec%paddedN/),kind=C_SIZE_T), mvec%val, ierr)
    if ( ierr /= 0 ) return

    ! that should hopefully help in cases of NUMA
    call dset_1(size(mvec%val), mvec%val, 0._8)

    mvec_ptr = c_loc(mvec)
    ierr = 0

  end subroutine phist_Dmvec_create


  subroutine phist_Dmvec_delete(mvec_ptr, ierr) bind(C,name='phist_Dmvec_delete_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value       :: mvec_ptr
    integer(C_INT),     intent(out) :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: mvec
#ifdef F_DEBUG
    integer(C_INTPTR_T) :: dummy
#endif
    !--------------------------------------------------------------------------------

#ifdef F_DEBUG
    write(*,*) 'deleting mvec at address', transfer(mvec_ptr,dummy)
    flush(6)
#endif
    if( c_associated(mvec_ptr) ) then
      call c_f_pointer(mvec_ptr, mvec)
      if( .not. mvec%is_view) then
        call deallocate_aligned(mvec%val)
      end if
      deallocate(mvec)
    end if
    ierr = 0

  end subroutine phist_Dmvec_delete


  subroutine phist_Dmvec_extract_view(mvec_ptr, raw_ptr, lda, ierr) bind(C,name='phist_Dmvec_extract_view_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value       :: mvec_ptr
    type(C_PTR),        intent(out) :: raw_ptr
    integer(C_INT32_T), intent(out) :: lda!, stride
    integer(C_INT),     intent(out) :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: mvec
#ifdef F_DEBUG
    integer(C_INTPTR_T) :: dummy
#endif
    !--------------------------------------------------------------------------------

#ifdef F_DEBUG
    write(*,*) 'extract view of mvec at address', transfer(mvec_ptr,dummy)
    flush(6)
#endif
    if( c_associated(mvec_ptr) ) then
      call c_f_pointer(mvec_ptr, mvec)
      raw_ptr = c_loc(mvec%val(mvec%jmin,1))
      lda = size(mvec%val,1)
      !stride = size(vec%val,2)
      ierr = 0
    else
      ierr = -88
    end if

  end subroutine phist_Dmvec_extract_view


  subroutine phist_Dmvec_my_length(mvec_ptr, mylen, ierr) bind(C,name='phist_Dmvec_my_length_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value       :: mvec_ptr
    integer(C_INT32_T), intent(out) :: mylen
    integer(C_INT),     intent(out) :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: mvec
    !--------------------------------------------------------------------------------

    if( c_associated(mvec_ptr) ) then
      call c_f_pointer(mvec_ptr, mvec)
      mylen = mvec%map%nlocal(mvec%map%me)
      ierr = 0
    else
      ierr = -88
    end if

  end subroutine phist_Dmvec_my_length


  subroutine phist_Dmvec_get_map(mvec_ptr, map_ptr, ierr) bind(C,name='phist_Dmvec_get_map_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value       :: mvec_ptr
    type(C_PTR),        intent(out) :: map_ptr
    integer(C_INT),     intent(out) :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: mvec
    !--------------------------------------------------------------------------------

    if( c_associated(mvec_ptr) ) then
      call c_f_pointer(mvec_ptr, mvec)
      map_ptr = c_loc(mvec%map)
      ierr = 0
    else
      ierr = -88
    end if

  end subroutine phist_Dmvec_get_map


  subroutine phist_Dmvec_num_vectors(mvec_ptr, nvec, ierr) bind(C,name='phist_Dmvec_num_vectors_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value       :: mvec_ptr
    integer(C_INT),     intent(out) :: nvec
    integer(C_INT),     intent(out) :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: mvec
    !--------------------------------------------------------------------------------

    if( c_associated(mvec_ptr) ) then
      call c_f_pointer(mvec_ptr, mvec)
      nvec = mvec%jmax-mvec%jmin+1
      ierr = 0
    else
      ierr = -88
    end if

  end subroutine phist_Dmvec_num_vectors


  subroutine phist_Dmvec_view_block(mvec_ptr, view_ptr, jmin, jmax, ierr) bind(C,name='phist_Dmvec_view_block_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: mvec_ptr
    type(C_PTR),        intent(inout) :: view_ptr
    integer(C_INT),     value         :: jmin, jmax
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: mvec, view
!#ifdef F_DEBUG
!    integer(C_INTPTR_T) :: dummy
!#endif
    !--------------------------------------------------------------------------------

    if( jmin .lt. 0 .or. jmax .lt. jmin ) then
      write(*,*) 'Trying to create a view with jmin', jmin, 'and jmax', jmax
      flush(6)
      ierr = -88
      return
    end if

!#ifdef F_DEBUG
!    write(*,*) 'create view of mvec at address', transfer(mvec_ptr,dummy)
!    flush(6)
!#endif
    if( c_associated(mvec_ptr) ) then
      call c_f_pointer(mvec_ptr, mvec)
      if( c_associated(view_ptr) ) then
        call c_f_pointer(view_ptr,view)
        if( .not. view%is_view ) then
          write(*,*) 'mvec not a view!'
          flush(6)
          ierr = -88
          return
        end if
!#ifdef F_DEBUG
!      write(*,*) 'reusing view at address', transfer(view_ptr,dummy)
!      flush(6)
!#endif
      else
        allocate(view)
        view_ptr = c_loc(view)
!#ifdef F_DEBUG
!      write(*,*) 'created new view at address', transfer(view_ptr,dummy)
!      flush(6)
!#endif
      end if
      view%is_view = .true.
      view%jmin = mvec%jmin+jmin
      view%jmax = mvec%jmin+jmax
      view%paddedN = mvec%paddedN
      view%map=mvec%map
      view%val=>mvec%val
      ierr = 0
    else
      ierr = -88
    end if

  end subroutine phist_Dmvec_view_block


  subroutine phist_Dmvec_get_block(mvec_ptr, block_ptr, jmin, jmax, ierr) bind(C,name='phist_Dmvec_get_block_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: mvec_ptr, block_ptr
    integer(C_INT),     value         :: jmin, jmax
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: mvec, block
    type(MVec_t) :: view
    !--------------------------------------------------------------------------------

    if( .not. c_associated(mvec_ptr) .or. .not. c_associated(block_ptr) ) then
      ierr = -88
      return
    end if

    call c_f_pointer(mvec_ptr, mvec)
    call c_f_pointer(block_ptr,block)

#ifdef F_DEBUG
    if( .not. map_compatible_map(mvec%map, block%map) ) then
      ierr = -1
      return
    end if
#endif

    if( jmax-jmin .ne. block%jmax-block%jmin .or. &
      & jmax-jmin .gt. mvec%jmax-mvec%jmin        ) then
      ierr = -1
      return
    end if

    ! create view and let mvec_add_mvec handle the rest!
    view%is_view = .true.
    view%map = mvec%map
    view%jmin = mvec%jmin+jmin
    view%jmax = mvec%jmin+jmax
    view%val=>mvec%val
    call mvec_add_mvec(1._8,view,0._8,block)

    ierr = 0

  end subroutine phist_Dmvec_get_block


  subroutine phist_Dmvec_set_block(mvec_ptr, block_ptr, jmin, jmax, ierr) bind(C,name='phist_Dmvec_set_block_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: mvec_ptr, block_ptr
    integer(C_INT),     value         :: jmin, jmax
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: mvec, block
    type(MVec_t) :: view
    !--------------------------------------------------------------------------------

    if( .not. c_associated(mvec_ptr) .or. .not. c_associated(block_ptr) ) then
      ierr = -88
      return
    end if

    call c_f_pointer(mvec_ptr, mvec)
    call c_f_pointer(block_ptr,block)

#ifdef TESTING
    if( .not. map_compatible_map(mvec%map, block%map) ) then
      ierr = -1
      return
    end if
#endif

    if( jmax-jmin .ne. block%jmax-block%jmin .or. &
      & jmax-jmin .gt. mvec%jmax-mvec%jmin        ) then
      ierr = -1
      return
    end if

    ! create view and let mvec_add_mvec handle the rest!
    view%is_view = .true.
    view%map = mvec%map
    view%jmin = mvec%jmin+jmin
    view%jmax = mvec%jmin+jmax
    view%val=>mvec%val
    view%paddedN = mvec%paddedN
    call mvec_add_mvec(1._8,block,0._8,view)

    ierr = 0

  end subroutine phist_Dmvec_set_block


  subroutine phist_Dmvec_gather_mvecs(mvec_ptr, block_ptr_list, nblocks, ierr) bind(C,name='phist_Dmvec_gather_mvecs_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: mvec_ptr
    type(C_PTR),        intent(in)    :: block_ptr_list(*)
    integer(C_INT),     value         :: nblocks
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: mvec, tmp
    type(MVec_t), allocatable :: block_list(:)
    integer :: i
    !--------------------------------------------------------------------------------

    if( .not. c_associated(mvec_ptr) ) then
      ierr = -88
      return
    end if
    call c_f_pointer(mvec_ptr, mvec)

    allocate(block_list(nblocks))
    do i = 1, nblocks, 1
      if( .not. c_associated(block_ptr_list(i)) ) then
        ierr = -88
        return
      end if
      call c_f_pointer(block_ptr_list(i),tmp)

#ifdef TESTING
    if( .not. map_compatible_map(mvec%map, tmp%map) ) then
      ierr = -1
      return
    end if
#endif

      block_list(i)%jmin = tmp%jmin
      block_list(i)%jmax = tmp%jmax
      block_list(i)%is_view = .true.
      block_list(i)%map = tmp%map
      block_list(i)%val => tmp%val
      block_list(i)%paddedN = tmp%paddedN
    end do

    call mvec_gather_mvecs(mvec,block_list)

    ierr = 0

  end subroutine phist_Dmvec_gather_mvecs


  subroutine phist_Dmvec_scatter_mvecs(mvec_ptr, block_ptr_list, nblocks, ierr) bind(C,name='phist_Dmvec_scatter_mvecs_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: mvec_ptr
    type(C_PTR),        intent(in)    :: block_ptr_list(*)
    integer(C_INT),     value         :: nblocks
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: mvec, tmp
    type(MVec_t), allocatable :: block_list(:)
    integer :: i
    !--------------------------------------------------------------------------------

    if( .not. c_associated(mvec_ptr) ) then
      ierr = -88
      return
    end if
    call c_f_pointer(mvec_ptr, mvec)

    allocate(block_list(nblocks))
    do i = 1, nblocks, 1
      if( .not. c_associated(block_ptr_list(i)) ) then
        ierr = -88
        return
      end if
      call c_f_pointer(block_ptr_list(i),tmp)

#ifdef TESTING
    if( .not. map_compatible_map(mvec%map, tmp%map) ) then
      ierr = -1
      return
    end if
#endif

      block_list(i)%jmin = tmp%jmin
      block_list(i)%jmax = tmp%jmax
      block_list(i)%is_view = .true.
      block_list(i)%map = tmp%map
      block_list(i)%val => tmp%val
      block_list(i)%paddedN = tmp%paddedN
    end do

    call mvec_scatter_mvecs(mvec,block_list)

    ierr = 0

  end subroutine phist_Dmvec_scatter_mvecs


  subroutine phist_Dmvec_put_value(mvec_ptr, val, ierr) &
    & bind(C,name='phist_Dmvec_put_value_f') ! circumvent bug in opari (openmp instrumentalization)
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: mvec_ptr
    real(C_DOUBLE),     value         :: val
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: mvec
    !--------------------------------------------------------------------------------

    if( .not. c_associated(mvec_ptr) ) then
      ierr = -88
      return
    end if

    call c_f_pointer(mvec_ptr, mvec)

    if( .not. mvec%is_view .or. &
      & ( mvec%jmin .eq. lbound(mvec%val,1) .and. &
      &   mvec%jmax .eq. ubound(mvec%val,1)       ) ) then

      call dset_1(size(mvec%val,1)*mvec%map%nlocal(mvec%map%me), mvec%val(1,1), val)
    else
      call dset_general(mvec%jmax-mvec%jmin+1, mvec%map%nlocal(mvec%map%me), mvec%val(mvec%jmin,1), size(mvec%val,1), val)
    end if

    ierr = 0

  end subroutine phist_Dmvec_put_value


  subroutine phist_Dmvec_put_func(V_ptr, elemFunc_ptr, last_arg, ierr) &
    & bind(C,name='phist_Dmvec_put_func_f') ! circumvent bug in opari (openmp instrumentalization)
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value :: V_ptr
    type(C_FUNPTR),     value       :: elemFunc_ptr
    type(C_PTR),        value       :: last_arg
    integer(C_INT),     intent(out) :: ierr
    !--------------------------------------------------------------------------------
    type(mvec_t), pointer :: V
    procedure(mvecElemFunc), pointer :: elemFunc
    !--------------------------------------------------------------------------------
    integer(kind=8) :: i
    integer(kind=4) :: j
    integer(kind=G_GIDX_T) :: ii
    integer(kind=G_LIDX_T) :: jj
    integer(C_INT) :: thread_ierr
    !--------------------------------------------------------------------------------

    ierr=0

    if( .not. c_associated(V_ptr) ) then
      ierr = -88
      return
    end if

    call c_f_pointer(V_ptr, V)

    ! get procedure pointer
    call c_f_procpointer(elemFunc_ptr, elemFunc)

! note: element function assumes 0-based indexing

!$omp parallel do schedule(static) private(ii,j,jj,thread_ierr)
    do i = 1, V%map%nlocal(V%map%me), 1
      ii = V%map%distrib(V%map%me)+i-2
      do j=V%jmin,V%jmax
        jj=j-V%jmin
        thread_ierr=elemFunc(ii,jj,C_LOC(V%val(j,i)),last_arg)
        !write(*,*) ii,jj,V%val(j,i)
        if (thread_ierr/=0) then
!$omp critical
          ierr=thread_ierr
!$omp end critical
        end if
      end do
    end do

end subroutine phist_Dmvec_put_func

  subroutine phist_Dmvec_random(mvec_ptr, ierr) bind(C,name='phist_Dmvec_random_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: mvec_ptr
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: mvec
    integer :: nvec, nrows, ldx
    integer(kind=8) :: pre_skip, post_skip
    !--------------------------------------------------------------------------------

    if( .not. c_associated(mvec_ptr) ) then
      ierr = -88
      return
    end if

    call c_f_pointer(mvec_ptr, mvec)
    call mvec_random(mvec)

    ierr = 0

  end subroutine phist_Dmvec_random


  subroutine phist_Dmvec_print(mvec_ptr, ierr) bind(C,name='phist_Dmvec_print_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: mvec_ptr
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: mvec
    integer :: i,j, iproc
    !--------------------------------------------------------------------------------

    if( .not. c_associated(mvec_ptr) ) then
      ierr = -88
      return
    end if

    call c_f_pointer(mvec_ptr, mvec)

    do iproc = 0, mvec%map%nProcs-1, 1

      if( iproc .eq. mvec%map%me ) then
        do i = 1, mvec%map%nlocal(mvec%map%me), 1
          do j=mvec%jmin,mvec%jmax
            write(*,'(G25.16,A2)',advance='no') mvec%val(j,i),'  '
          end do
          write(*,*)
        end do
        flush(6)
      end if

      call MPI_Barrier(mvec%map%comm, ierr)
    end do

    ierr = 0

  end subroutine phist_Dmvec_print


  subroutine phist_Dmvec_norm2(mvec_ptr, vnrm, iflag) bind(C,name='phist_Dmvec_norm2_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: mvec_ptr
    real(C_DOUBLE),     intent(out)   :: vnrm(*)
    integer(C_INT),     intent(inout) :: iflag
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: mvec
    integer :: nvec
    !--------------------------------------------------------------------------------

    if( .not. c_associated(mvec_ptr) ) then
      iflag = -88
      return
    end if

    call c_f_pointer(mvec_ptr, mvec)

    nvec = mvec%jmax-mvec%jmin+1
    call mvec_norm2(mvec, vnrm(1:nvec), iflag)
    
  end subroutine phist_Dmvec_norm2


  subroutine phist_Dmvec_scale(x_ptr, alpha, ierr) bind(C,name='phist_Dmvec_scale_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: x_ptr
    real(C_DOUBLE),     value         :: alpha
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: x
    !--------------------------------------------------------------------------------

    if( .not. c_associated(x_ptr) ) then
      ierr = -88
      return
    end if

    call c_f_pointer(x_ptr, x)

    call mvec_scale(x,alpha)

    ierr = 0

  end subroutine phist_Dmvec_scale


  subroutine phist_Dmvec_vscale(x_ptr, alpha, ierr) bind(C,name='phist_Dmvec_vscale_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: x_ptr
    real(C_DOUBLE),     intent(in)    :: alpha(*)
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: x
    !--------------------------------------------------------------------------------

    if( .not. c_associated(x_ptr) ) then
      ierr = -88
      return
    end if

    call c_f_pointer(x_ptr, x)

    call mvec_vscale(x,alpha)

    ierr = 0

  end subroutine phist_Dmvec_vscale



  subroutine phist_Dmvec_add_mvec(alpha, x_ptr, beta, y_ptr, ierr) bind(C,name='phist_Dmvec_add_mvec_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    real(C_DOUBLE),     value         :: alpha, beta
    type(C_PTR),        value         :: x_ptr, y_ptr
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: x, y
    !--------------------------------------------------------------------------------

    if( .not. c_associated(y_ptr) ) then
      ierr = -88
      return
    end if
    if( .not. c_associated(x_ptr) ) then
      ierr = -88
      return
    end if

    call c_f_pointer(x_ptr, x)
    call c_f_pointer(y_ptr, y)

    if( x%jmin .lt. 1 ) then
      ierr = -1
      return
    end if
    if( y%jmin .lt. 1 ) then
      ierr = -1
      return
    end if
    if( y%jmax-y%jmin .ne. x%jmax-x%jmin ) then
      ierr = -1
      return
    end if

#ifdef TESTING
    if( .not. map_compatible_map(x%map, y%map) ) then
      ierr = -1
      return
    end if
#endif

    call mvec_add_mvec(alpha, x, beta, y)
    ierr = 0

  end subroutine phist_Dmvec_add_mvec


  subroutine phist_Dmvec_vadd_mvec(alpha, x_ptr, beta, y_ptr, ierr) bind(C,name='phist_Dmvec_vadd_mvec_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    real(C_DOUBLE),     intent(in)    :: alpha(*)
    type(C_PTR),        value         :: x_ptr, y_ptr
    real(C_DOUBLE),     value         :: beta
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: x, y
    integer :: nvec
    !--------------------------------------------------------------------------------

    !if( beta .ne. 0 .and. .not. c_associated(y_ptr) ) then
    if( .not. c_associated(y_ptr) ) then
      ierr = -88
      return
    end if

    !if( any(alpha(1:nvec) .ne. 0) .and. .not. c_associated(x_ptr) ) then
    if( .not. c_associated(x_ptr) ) then
      ierr = -88
      return
    end if

    call c_f_pointer(x_ptr, x)
    call c_f_pointer(y_ptr, y)

    if( y%jmin .lt. 1 ) then
      ierr = -1
      return
    end if
    if( x%jmin .lt. 1 ) then
      ierr = -1
      return
    end if

    nvec = y%jmax-y%jmin+1
    if( nvec .ne. x%jmax-x%jmin+1 ) then
      ierr = -1
      return
    end if

#ifdef TESTING
    if( .not. map_compatible_map(x%map, y%map) ) then
      ierr = -1
      return
    end if
#endif


    call mvec_vadd_mvec(alpha, x, beta, y)
    ierr = 0

  end subroutine phist_Dmvec_vadd_mvec


  subroutine phist_Dmvec_dot_mvec(x_ptr, y_ptr, dot, iflag) bind(C,name='phist_Dmvec_dot_mvec_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: x_ptr, y_ptr
    real(C_DOUBLE),     intent(out)   :: dot(*)
    integer(C_INT),     intent(inout) :: iflag
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: x, y
    !--------------------------------------------------------------------------------

    if( .not. c_associated(x_ptr) .or. .not. c_associated(y_ptr) ) then
      iflag = -88
      return
    end if
    call c_f_pointer(x_ptr, x)
    call c_f_pointer(y_ptr, y)

    if( y%jmax-y%jmin .ne. x%jmax-x%jmin ) then
      iflag = -1
      return
    end if

#ifdef TESTING
    if( .not. map_compatible_map(x%map, y%map) ) then
      iflag = -1
      return
    end if
#endif

    call mvec_dot_mvec(x,y,dot,iflag)

  end subroutine phist_Dmvec_dot_mvec


  subroutine phist_Dmvec_times_sdMat(alpha, v_ptr, M_ptr, beta, w_ptr, iflag) bind(C,name='phist_Dmvec_times_sdMat_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: v_ptr, w_ptr
    real(C_DOUBLE),     value         :: alpha, beta
    type(C_PTR),        value         :: M_ptr
    integer(C_INT),     intent(inout) :: iflag
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: v, w
    type(SDMat_t), pointer :: M
    !--------------------------------------------------------------------------------

    if( .not. c_associated(v_ptr) .or. &
      & .not. c_associated(w_ptr) .or. &
      & .not. c_associated(M_ptr)      ) then
      iflag = -88
      return
    end if

    call c_f_pointer(v_ptr,v)
    call c_f_pointer(w_ptr,w)
    call c_f_pointer(M_ptr,M)

    if( v%jmax-v%jmin .ne. M%imax-M%imin .or. &
      & M%jmax-M%jmin .ne. w%jmax-w%jmin      ) then
      iflag = -1
      return
    end if

#ifdef TESTING
    if( .not. map_compatible_map(v%map, w%map) ) then
      iflag = -1
      return
    end if
#endif

    call mvec_times_sdmat(alpha,v,M,beta,w,iflag)

  end subroutine phist_Dmvec_times_sdMat


  subroutine phist_Dmvec_times_sdMat_augmented(alpha, v_ptr, M_ptr, beta, w_ptr, N_ptr, iflag) bind(C,name='phist_Dmvec_times_sdMat_augmented_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: v_ptr, w_ptr
    real(C_DOUBLE),     value         :: alpha, beta
    type(C_PTR),        value         :: M_ptr, N_ptr
    integer(C_INT),     intent(inout) :: iflag
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: v, w
    type(SDMat_t), pointer :: M, N
    !--------------------------------------------------------------------------------

    if( .not. c_associated(v_ptr) .or. &
      & .not. c_associated(w_ptr) .or. &
      & .not. c_associated(M_ptr) .or. &
      & .not. c_associated(N_ptr)      ) then
      iflag = -88
      return
    end if

    call c_f_pointer(v_ptr,v)
    call c_f_pointer(w_ptr,w)
    call c_f_pointer(M_ptr,M)
    call c_f_pointer(N_ptr,N)

    if( v%jmax-v%jmin .ne. M%imax-M%imin .or. &
      & M%jmax-M%jmin .ne. w%jmax-w%jmin .or. &
      & N%imax-N%imin .ne. N%jmax-N%jmin .or. &
      & w%jmax-w%jmin .ne. N%imax-N%imin      ) then
      iflag = -1
      return
    end if

#ifdef TESTING
    if( .not. map_compatible_map(v%map, w%map) ) then
      iflag = -1
      return
    end if
#endif

    call mvec_times_sdmat_augmented(alpha,v,M,beta,w,N,iflag)

  end subroutine phist_Dmvec_times_sdMat_augmented


  subroutine phist_Dmvec_times_sdMat_add_mvec_times_sdMat(v_ptr, M_ptr, w_ptr, N_ptr, iflag) bind(C,name='phist_Dmvec_times_sdMat_add_mvec_times_sdMat_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: v_ptr, w_ptr
    type(C_PTR),        value         :: M_ptr, N_ptr
    integer(C_INT),     intent(inout) :: iflag
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: v, w
    type(SDMat_t), pointer :: M, N
    !--------------------------------------------------------------------------------

    if( .not. c_associated(v_ptr) .or. &
      & .not. c_associated(w_ptr) .or. &
      & .not. c_associated(M_ptr) .or. &
      & .not. c_associated(N_ptr)      ) then
      iflag = -88
      return
    end if

    call c_f_pointer(v_ptr,v)
    call c_f_pointer(w_ptr,w)
    call c_f_pointer(M_ptr,M)
    call c_f_pointer(N_ptr,N)

    if( v%jmax-v%jmin .ne. M%imax-M%imin .or. &
      & M%jmax-M%jmin .ne. w%jmax-w%jmin .or. &
      & N%imax-N%imin .ne. N%jmax-N%jmin .or. &
      & w%jmax-w%jmin .ne. N%imax-N%imin      ) then
      iflag = -1
      return
    end if

#ifdef TESTING
    if( .not. map_compatible_map(v%map, w%map) ) then
      iflag = -1
      return
    end if
#endif

    call mvec_times_sdMat_add_mvec_times_sdMat(v,M,w,N,iflag)

  end subroutine phist_Dmvec_times_sdMat_add_mvec_times_sdMat



  subroutine phist_Dmvec_times_sdMat_inplace(v_ptr, M_ptr, iflag) bind(C,name='phist_Dmvec_times_sdMat_inplace_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: v_ptr
    type(C_PTR),        value         :: M_ptr
    integer(C_INT),     intent(inout) :: iflag
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: v
    type(SDMat_t), pointer :: M
    !--------------------------------------------------------------------------------

    if( .not. c_associated(v_ptr) .or. &
      & .not. c_associated(M_ptr)      ) then
      iflag = -88
      return
    end if

    call c_f_pointer(v_ptr,v)
    call c_f_pointer(M_ptr,M)

    if( v%jmax-v%jmin .ne. M%imax-M%imin .or. &
      & M%jmax-M%jmin .gt. v%jmax-v%jmin      ) then
      iflag = -1
      return
    end if

    call mvec_times_sdmat_inplace(v,M,iflag)

  end subroutine phist_Dmvec_times_sdMat_inplace


  subroutine phist_DmvecT_times_mvec(alpha, v_ptr, w_ptr, beta, M_ptr, iflag) bind(C,name='phist_DmvecT_times_mvec_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: v_ptr, w_ptr, M_ptr
    real(C_DOUBLE),     value         :: alpha, beta
    integer(C_INT),     intent(inout) :: iflag
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: v, w
    type(SDMat_t), pointer :: M
    !--------------------------------------------------------------------------------

    if( .not. c_associated(v_ptr) .or. &
      & .not. c_associated(w_ptr) .or. &
      & .not. c_associated(M_ptr)      ) then
      iflag = -88
      return
    end if

    call c_f_pointer(v_ptr,v)
    call c_f_pointer(w_ptr,w)
    call c_f_pointer(M_ptr,M)

    if( v%jmax-v%jmin .ne. M%imax-M%imin .or. &
      & M%jmax-M%jmin .ne. w%jmax-w%jmin      ) then
      iflag = -1
      return
    end if

#ifdef TESTING
    if( .not. map_compatible_map(v%map, w%map) ) then
      iflag = -1
      return
    end if
#endif

    call mvecT_times_mvec(alpha,v,w,beta,M,iflag)

  end subroutine phist_DmvecT_times_mvec


  subroutine phist_DmvecT_times_mvec_times_sdMat_inplace(alpha, v_ptr, w_ptr, N_ptr, beta, M_ptr, iflag) bind(C,name='phist_DmvecT_times_mvec_times_sdMat_inplace_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: v_ptr, w_ptr, M_ptr, N_ptr
    real(C_DOUBLE),     value         :: alpha, beta
    integer(C_INT),     intent(inout) :: iflag
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: v, w
    type(SDMat_t), pointer :: M, N
    !--------------------------------------------------------------------------------

    if( .not. c_associated(v_ptr) .or. &
      & .not. c_associated(w_ptr) .or. &
      & .not. c_associated(N_ptr) .or. &
      & .not. c_associated(M_ptr)      ) then
      iflag = -88
      return
    end if

    call c_f_pointer(v_ptr,v)
    call c_f_pointer(w_ptr,w)
    call c_f_pointer(M_ptr,M)
    call c_f_pointer(N_ptr,N)

    if( v%jmax-v%jmin .ne. M%imax-M%imin .or. &
      & M%jmax-M%jmin .ne. w%jmax-w%jmin .or. &
      & N%imax-N%imin .ne. N%jmax-N%jmin .or. &
      & w%jmax-w%jmin .ne. N%imax-N%imin      ) then
      iflag = -1
      return
    end if

#ifdef TESTING
    if( .not. map_compatible_map(v%map, w%map) ) then
      iflag = -1
      return
    end if
#endif

    call mvecT_times_mvec_times_sdMat_inplace(alpha,v,w,N,beta,M,iflag)

  end subroutine phist_DmvecT_times_mvec_times_sdMat_inplace



  subroutine phist_Dmvec_to_mvec(v_ptr, w_ptr, ierr) bind(C,name='phist_Dmvec_to_mvec_f')
    use, intrinsic :: iso_c_binding
    use mpi
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: v_ptr, w_ptr
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: v, w
    !--------------------------------------------------------------------------------
    integer :: nvec, i, nglobal
    integer(kind=8), allocatable :: globalIdx(:), invGlobalIdx(:)
    integer, allocatable :: counts(:), offsets(:)
    real(kind=8), allocatable :: globalVec(:,:), sendBuff(:,:)
    !--------------------------------------------------------------------------------

    if( .not. c_associated(v_ptr) .or. &
      & .not. c_associated(w_ptr)      ) then
      ierr = -88
      return
    end if

    call c_f_pointer(v_ptr,v)
    call c_f_pointer(w_ptr,w)

    if( v%jmax-v%jmin .ne. w%jmax-w%jmin ) then
      ierr = -1
      return
    end if

    ! we assume v and w have different maps...
    if( .not. map_compatible_map(v%map, w%map, reorder=.true.) ) then
      ierr = -1
      return
    end if

    ! we do this in two steps, first put all elements in ascending order, then redistribute them
    ! (this is not optimal, but we shouldn't need this subroutine very often)
    nvec = v%jmax-v%jmin+1
    nglobal = v%map%distrib(v%map%nProcs)-1
    if( allocated(v%map%global_idx) ) then

      allocate(globalIdx(nglobal))
      ! get all data globally! doesn't work for too large matrices
      call MPI_Allgatherv(v%map%global_idx, v%map%nlocal(v%map%me), MPI_INTEGER8, &
        &                 globalIdx, v%map%nlocal, int(v%map%distrib-1), MPI_INTEGER8, &
        &                 v%map%comm, ierr)

      ! get the inverse global index
      allocate(invGlobalIdx(nglobal))
      !write(*,*) 'distributed globalIdx', globalIdx
      do i = 1, nglobal
        invGlobalIdx(globalIdx(i)) = i
      end do
      !write(*,*) 'calculated invGlobalIdx', invGlobalIdx

    else

      allocate(invGlobalIdx(nglobal))
      do i = 1, nglobal
        invGlobalIdx(i) = i
      end do
      !write(*,*) 'generated invGlobalIdx', invGlobalIdx

    end if

    ! get all data globally! doesn't work for too large matrices
    allocate(globalVec(nvec,v%map%distrib(v%map%nProcs)-1))
    allocate(counts(0:size(v%map%nlocal)-1))
    counts = v%map%nlocal*nvec
    allocate(offsets(0:size(v%map%distrib)-1))
    offsets = int(v%map%distrib-1)*nvec
    allocate(sendBuff(nvec,v%map%nlocal(v%map%me)))
    sendBuff = v%val(v%jmin:v%jmax,1:v%map%nlocal(v%map%me))
    call MPI_Allgatherv(sendBuff, counts(v%map%me), MPI_DOUBLE_PRECISION, &
      &                 globalVec, counts, offsets, MPI_DOUBLE_PRECISION, &
      &                 v%map%comm, ierr)


    ! copy data to target vector w
    if( allocated(w%map%global_idx) ) then
      do i = 1, w%map%nlocal(w%map%me)
        w%val(w%jmin:w%jmax,i) = globalVec(:,invGlobalIdx(w%map%global_idx(i)))
      end do
    else
      do i = 1, w%map%nlocal(w%map%me)
        w%val(w%jmin:w%jmax,i) = globalVec(:,invGlobalIdx(i+w%map%distrib(w%map%me)-1))
      end do
    end if

    ierr = 0

  end subroutine phist_Dmvec_to_mvec


  subroutine phist_Dmvec_QR(v_ptr, R_ptr, ierr) bind(C,name='phist_Dmvec_QR_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: v_ptr, R_ptr
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: v
    type(SDMat_t), pointer :: R
    !--------------------------------------------------------------------------------

    if( .not. c_associated(v_ptr) .or. &
      & .not. c_associated(R_ptr)      ) then
      ierr = -88
      return
    end if

    call c_f_pointer(v_ptr,v)
    call c_f_pointer(R_ptr,R)

    if( v%jmax-v%jmin .ne. R%imax-R%imin .or. &
      & v%jmax-v%jmin .ne. R%jmax-R%jmin      ) then
      ierr = -1
      return
    end if

    call mvec_QR(v,R,ierr)

  end subroutine phist_Dmvec_QR

end module mvec_module

