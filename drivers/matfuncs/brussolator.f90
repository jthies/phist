#include "phist_config.h"
#ifdef PHIST_HAVE_GHOST
#include "ghost/config.h"
#endif

#ifndef PHIST_HAVE_GHOST
#define G_LIDX_T C_INT32_T
#define G_GIDX_T C_INT64_T
#else
#ifdef GHOST_HAVE_LONGIDX_LOCAL
#define G_LIDX_T C_INT64_T
#else
#define G_LIDX_T C_INT32_T
#endif
#ifdef GHOST_HAVE_LONGIDX_GLOBAL
#define G_GIDX_T C_INT64_T
#else
#define G_GIDX_T C_INT32_T
#endif
#endif

! simple symmetric positive definite tridiagonal Toeplitz matrix (e.g. 1D Poisson discretization)
! From the NEP collection: MVMBWM, a matrix generator for a brusselator wave model in chemical reaction
module brussolator_module
  implicit none
  private

  public :: Brussolator_rowFunc
  public :: Brussolator_initDimensions

  integer(kind=8) :: m
  real(kind=8) :: tau1, tau2

  real(kind=8), parameter :: delt1 = 0.008
  real(kind=8), parameter :: delt2 = 0.5*delt1
  real(kind=8), parameter :: alpha = 2.
  real(kind=8), parameter :: beta  = 5.45
  real(kind=8), parameter :: L     = 0.51302

contains

  subroutine Brussolator_initDimensions(newM, nrows, maxnne_per_row) bind(C, name='Brussolator_initDimensions')
    use, intrinsic :: iso_c_binding
    integer(kind=C_INT), value :: newM
    integer(kind=G_GIDX_T), intent(out) :: nrows
    integer(kind=G_LIDX_T), intent(out) :: maxnne_per_row
    real(kind=8) :: h

    ! dimensions
    m = newM
    nrows = 2*m
    maxnne_per_row = 4

    ! precalculate some values
    h = 1._8/(m+1)
    tau1 = delt1/(h**2*L**2)
    tau2 = delt2/(h**2*L**2)

  end subroutine Brussolator_initDimensions


  function Brussolator_rowFunc(row, nnz, cols, vals) result(the_result) bind(C, name='Brussolator_rowFunc')
    use, intrinsic :: iso_c_binding
    integer(G_GIDX_T), value :: row
    integer(G_LIDX_T), intent(inout) :: nnz
    integer(G_GIDX_T), intent(inout) :: cols(*)
    real(C_DOUBLE),    intent(inout) :: vals(*)
    integer(C_INT) :: the_result

    nnz = 0

    if( row < m ) then

      ! lower diagonal
      if( row .gt. 0 ) then
        nnz = nnz + 1
        vals(nnz) = tau1
        cols(nnz) = row-1
      end if

      ! diagonal
      nnz = nnz + 1
      vals(nnz) = -2._8*tau1+beta-1
      cols(nnz) = row

      ! upper diagonal
      nnz = nnz + 1
      vals(nnz) = tau1
      cols(nnz) = row+1

      ! upper block
      nnz = nnz + 1
      vals(nnz) = alpha**2
      cols(nnz) = row+m

    else ! row >= m

      ! lower block
      nnz = nnz + 1
      vals(nnz) = -beta
      cols(nnz) = row-m

      ! lower diagonal
      nnz = nnz + 1
      vals(nnz) = tau2
      cols(nnz) = row-1

      ! diagonal
      nnz = nnz + 1
      vals(nnz) = -2._8*tau2-alpha**2
      cols(nnz) = row

      ! upper diagonal
      if( row+1 .lt. 2*m ) then
        nnz = nnz + 1
        vals(nnz) = tau2
        cols(nnz) = row+1
      end if

    end if


    the_result = 0
    
  end function Brussolator_rowFunc

end module brussolator_module
