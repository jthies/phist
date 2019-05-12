/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/

#include "phist_config_fortran.h"

! simple symmetric positive definite tridiagonal Toeplitz matrix (e.g. 1D Poisson discretization)
module tritoeplitz_module
  implicit none
  private

  public :: TriToeplitz_rowFunc
  public :: TriToeplitz_initDimensions

  integer(kind=8) :: n

  real, parameter :: diagonalDominance = 1.e-8

contains

  subroutine TriToeplitz_initDimensions(newN, nrows, maxnne_per_row) bind(C, name='TriToeplitz_initDimensions')
    use, intrinsic :: iso_c_binding
    integer(kind=C_INT), value :: newN
    integer(kind=G_GIDX_T), intent(out) :: nrows
    integer(kind=G_LIDX_T), intent(out) :: maxnne_per_row

    n = 2_8**newN
    nrows = n
    maxnne_per_row = 3

  end subroutine TriToeplitz_initDimensions


  function TriToeplitz_rowFunc(row, nnz, cols, vals) result(the_result) bind(C, name='TriToeplitz_rowFunc')
    use, intrinsic :: iso_c_binding
    integer(G_GIDX_T), value :: row
    integer(G_LIDX_T), intent(inout) :: nnz
    integer(G_GIDX_T), intent(inout) :: cols(*)
    real(C_DOUBLE),    intent(inout) :: vals(*)
    integer(C_INT) :: the_result

    nnz = 0

    ! lower diagonal
    if( row .gt. 0 ) then
      nnz = nnz + 1
      vals(nnz) = -1./(1+diagonalDominance)
      cols(nnz) = row-1
    end if

    ! diagonal
    nnz = nnz + 1
    vals(nnz) = 2.
    cols(nnz) = row

    ! upper diagonal
    if( row+1 .lt. n ) then
      nnz = nnz + 1
      vals(nnz) = -1./(1+diagonalDominance)
      cols(nnz) = row+1
    end if

    the_result = 0
    
  end function TriToeplitz_rowFunc

end module tritoeplitz_module
