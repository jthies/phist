#include "phist_kernels.h"
#include "phist_macros.h"

// these implementations can be used if a kernel package
// does not implement all four data types (cf. epetra/ for
// an example, which supports only double (D))

void _SUBR_(type_avail)(int *ierr)
  {
  *ierr=-99;
  }

// \name Matrix input from a file

//@{

//! read a matrix from a MatrixMarket (ASCII) file
void _SUBR_(crsMat_read_mm)(_TYPE_(crsMat_ptr)* A, const char* filename,int* ierr)
  {
  *ierr=-99;
  }

//! read a matrix from a Ghost CRS (binary) file.
void _SUBR_(crsMat_read_bin)(_TYPE_(crsMat_ptr)* A, const char* filename,int* ierr)
  {
  *ierr=-99;
  }

//! read a matrix from a Harwell-Boeing (HB) file
void _SUBR_(crsMat_read_hb)(_TYPE_(crsMat_ptr)* A, const char* filename,int* ierr)
  {
  *ierr=-99;
  }
//!@}

//! \name get information about the data distribution in a matrix (maps)

//!@{
//! get the row distribution of the matrix
void _SUBR_(crsMat_get_row_map)(_TYPE_(const_crsMat_ptr) A, const_map_ptr_t* map, int* ierr)
  {
  *ierr=-99;
  }

//! get column distribution of a matrix
void _SUBR_(crsMat_get_col_map)(_TYPE_(const_crsMat_ptr) A, const_map_ptr_t* map, int* ierr)
  {
  *ierr=-99;
  }

//! get the map for vectors x in y=A*x
void _SUBR_(crsMat_get_domain_map)(_TYPE_(const_crsMat_ptr) A, const_map_ptr_t* map, int* ierr)
  {
  *ierr=-99;
  }

//! get the map for vectors y in y=A*x
void _SUBR_(crsMat_get_range_map)(_TYPE_(const_crsMat_ptr) A, const_map_ptr_t* map, int* ierr)
  {
  *ierr=-99;
  }
//@}

//! \name constructors

//@{
//! create a block-vector. The entries are stored contiguously
//! at val in column major ordering.
void _SUBR_(mvec_create)(_TYPE_(mvec_ptr)* V, 
        const_map_ptr_t map, int nvec, int* ierr)
  {
  *ierr=-99;
  }

//! create a serial dense n x m matrix on all procs, with column major
//! ordering.
void _SUBR_(sdMat_create)(_TYPE_(sdMat_ptr)* M, 
int n, int m, int* ierr)
  {
  *ierr=-99;
  }

//@}

//! retrieve local length of the vectors in V
void _SUBR_(mvec_my_length)(_TYPE_(const_mvec_ptr) V, int* len, int* ierr)
  {
  *ierr=-99;
  }

//! retrieve number of vectors/columns in V
void _SUBR_(mvec_num_vectors)(_TYPE_(const_mvec_ptr) V, int* nvec, int* ierr)
  {
  *ierr=-99;
  }

//!
void _SUBR_(mvec_extract_view)(_TYPE_(mvec_ptr) V, _ST_** val, int* lda, int* ierr)
  {
  *ierr=-99;
  }

//!
void _SUBR_(sdMat_extract_view)(_TYPE_(sdMat_ptr) V, _ST_** val, int* lda, int* ierr)
  {
  *ierr=-99;
  }

//! \name destructors

//@{

//!
void _SUBR_(crsMat_delete)(_TYPE_(crsMat_ptr) A, int* ierr)
  {
  *ierr=-99;
  }

//!
void _SUBR_(mvec_delete)(_TYPE_(mvec_ptr) V, int* ierr)
  {
  *ierr=-99;
  }

//!
void _SUBR_(sdMat_delete)(_TYPE_(sdMat_ptr) M, int* ierr)
  {
  *ierr=-99;
  }

//@}

//! put scalar value into all elements of a multi-vector
void _SUBR_(mvec_put_value)(_TYPE_(mvec_ptr) V, _ST_ value, int* ierr)
  {
  *ierr=-99;
  }

//! put scalar value into all elements of a serial dense matrix
void _SUBR_(sdMat_put_value)(_TYPE_(mvec_ptr) V, _ST_ value, int* ierr)
  {
  *ierr=-99;
  }

//! put random numbers into all elements of a multi-vector
void _SUBR_(mvec_random)(_TYPE_(mvec_ptr) V, int* ierr)
  {
  *ierr=-99;
  }

//! put random numbers into all elements of a serial dense matrix
void _SUBR_(sdMat_random)(_TYPE_(sdMat_ptr) M, int* ierr)
  {
  *ierr=-99;
  }


//! \name Numerical functions
//!@{

//! y=alpha*x+beta*y
void _SUBR_(mvec_add_mvec)(_ST_ alpha, _TYPE_(const_mvec_ptr) X,
                            _ST_ beta,  _TYPE_(mvec_ptr)       Y, 
                            int* ierr)
  {
  *ierr=-99;
  }

//! y=alpha*A*x+beta*y.
void _SUBR_(crsMat_times_mvec)(_ST_ alpha, _TYPE_(const_crsMat_ptr) A, 
        _TYPE_(const_mvec_ptr) x, _ST_ beta, _TYPE_(mvec_ptr) y, int* ierr)
  {
  *ierr=-99;
  }

//! dot product of vectors v_i and w_i, i=1..numvecs
void _SUBR_(mvec_dot_mvec)(_TYPE_(const_mvec_ptr) v, 
                            _TYPE_(const_mvec_ptr) w, 
                            _ST_* s, int* ierr)
  {
  *ierr=-99;
  }
  
//! n x m multi-vector times m x m dense matrix gives n x m multi-vector,
//! W=alpha*V*C + beta*W
void _SUBR_(mvec_times_sdMat)(_ST_ alpha, _TYPE_(const_mvec_ptr) V,
                                       _TYPE_(const_sdMat_ptr) C,
                                       _ST_ beta,
                                       _TYPE_(mvec_ptr) W, int* ierr)
  {
  *ierr=-99;
  }

//! 'tall skinny' QR decomposition, V=Q*R, Q'Q=I, R upper triangular.
void _SUBR_(mvec_QR)(_TYPE_(const_mvec_ptr) V, _TYPE_(mvec_ptr) Q, _TYPE_(sdMat_ptr) R, int* ierr)
  {
  *ierr=-99;
  }

//!@}

