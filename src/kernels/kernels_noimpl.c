#include "essex_kernels.h"

// these implementations can be used if a kernel package
// does not implement all four data types (cf. epetra/ for
// an example, which supports only double (D))

// \name Matrix input from a file

//@{

//! read a matrix from a MatrixMarket (ASCII) file
_SUBROUTINE_(crsMat_read_mm)(_TYPE_(crsMat_ptr)* A, const char* filename,int* ierr)
  {
  *ierr=-99;
  }

//! read a matrix from a Ghost CRS (binary) file.
_SUBROUTINE_(crsMat_read_bin)(_TYPE_(crsMat_ptr)* A, const char* filename,int* ierr)
  {
  *ierr=-99;
  }

//! read a matrix from a Harwell-Boeing (HB) file
_SUBROUTINE_(crsMat_read_hb)(_TYPE_(crsMat_ptr)* A, const char* filename,int* ierr)
  {
  *ierr=-99;
  }
//!@}

//! \name get information about the data distribution in a matrix (maps)

//!@{
//! get the row distribution of the matrix
_SUBROUTINE_(crsMat_get_row_map)(_TYPE_(const_crsMat_ptr) A, const_map_ptr_t* map, int* ierr)
  {
  *ierr=-99;
  }

//! get column distribution of a matrix
_SUBROUTINE_(crsMat_get_col_map)(_TYPE_(const_crsMat_ptr) A, const_map_ptr_t* map, int* ierr)
  {
  *ierr=-99;
  }

//! get the map for vectors x in y=A*x
_SUBROUTINE_(crsMat_get_domain_map)(_TYPE_(const_crsMat_ptr) A, const_map_ptr_t* map, int* ierr)
  {
  *ierr=-99;
  }

//! get the map for vectors y in y=A*x
_SUBROUTINE_(crsMat_get_range_map)(_TYPE_(const_crsMat_ptr) A, const_map_ptr_t* map, int* ierr)
  {
  *ierr=-99;
  }
//@}

//! \name constructors

//@{
//! create a block-vector. The entries are stored contiguously
//! at val in column major ordering.
_SUBROUTINE_(mvec_create)(const_map_ptr_t map, int nvec, 
        _TYPE_(mvec_ptr)* V, int* ierr)
  {
  *ierr=-99;
  }

//! create a serial dense n x m matrix on all procs, with column major
//! ordering.
_SUBROUTINE_(sdMat_create)(int n, int m, _TYPE_(sdMat_ptr)* M, int* ierr)
  {
  *ierr=-99;
  }

//@}

//! retrieve local length of the vectors in V
_SUBROUTINE_(mvec_my_length)(_TYPE_(const_mvec_ptr) V, int* len, int* ierr)
  {
  *ierr=-99;
  }

//!
_SUBROUTINE_(mvec_extract_view)(_TYPE_(mvec_ptr) V, _ST_** val, int vec, int* ierr)
  {
  *ierr=-99;
  }

//!
_SUBROUTINE_(sdMat_extract_view)(_TYPE_(mvec_ptr) V, _ST_** val, int* ierr)
  {
  *ierr=-99;
  }

//! \name destructors

//@{

//!
_SUBROUTINE_(crsMat_delete)(_TYPE_(crsMat_ptr) A, int* ierr)
  {
  *ierr=-99;
  }

//!
_SUBROUTINE_(mvec_delete)(_TYPE_(mvec_ptr) V, int* ierr)
  {
  *ierr=-99;
  }

//!
_SUBROUTINE_(sdMat_delete)(_TYPE_(sdMat_ptr) M, int* ierr)
  {
  *ierr=-99;
  }

//@}

//! put scalar value into all elements of a multi-vector
_SUBROUTINE_(mvec_put_value)(_TYPE_(mvec_ptr) V, _ST_ value, int* ierr)
  {
  *ierr=-99;
  }

//! put scalar value into all elements of a serial dense matrix
_SUBROUTINE_(sdMat_put_value)(_TYPE_(mvec_ptr) V, _ST_ value, int* ierr)
  {
  *ierr=-99;
  }



//! \name Numerical functions
//!@{

//! y=alpha*A*x+beta*y.
_SUBROUTINE_(crsMat_X_mvec)(_ST_ alpha, _TYPE_(const_crsMat_ptr) A, 
        _TYPE_(const_mvec_ptr) x, _ST_ beta, _TYPE_(mvec_ptr) y, int* ierr)
  {
  *ierr=-99;
  }

//! dot product of vectors v_i and w_i, i=1..numvecs
_SUBROUTINE_(mvec_dot_mvec)(_TYPE_(const_mvec_ptr) v, 
                            _TYPE_(const_mvec_ptr) w, 
                            _ST_* s, int* ierr)
  {
  *ierr=-99;
  }
  
//! n x m multi-vector times m x m dense matrix gives n x m multi-vector,
//! W=alpha*V*C + beta*W
_SUBROUTINE_(mvec_X_sdMat)(_ST_ alpha, _TYPE_(const_mvec_ptr) V,
                                       _TYPE_(const_sdMat_ptr) C,
                                       _ST_ beta,
                                       _TYPE_(mvec_ptr) W, int* ierr)
  {
  *ierr=-99;
  }

//! 'tall skinny' QR decomposition, V=Q*R, Q'Q=I, R upper triangular.
_SUBROUTINE_(mvec_QR)(_TYPE_(const_mvec_ptr) V, _TYPE_(mvec_ptr) Q, _TYPE_(sdMat_ptr) R, int* ierr)
  {
  *ierr=-99;
  }

//!@}

