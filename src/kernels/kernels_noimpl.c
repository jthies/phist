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
        const_map_ptr_t map, lidx_t nvec, int* ierr)
  {
  *ierr=-99;
  }

//! create a serial dense n x m matrix on all procs, with column major
//! ordering.
void _SUBR_(sdMat_create)(_TYPE_(sdMat_ptr)* M, 
int nrows, int ncols, int* ierr)
  {
  *ierr=-99;
  }

//@}

//! retrieve local length of the vectors in V
void _SUBR_(mvec_my_length)(_TYPE_(const_mvec_ptr) V, lidx_t* len, int* ierr)
  {
  *ierr=-99;
  }

//! retrieve the map of the vectors in V
void _SUBR_(mvec_get_map)(_TYPE_(const_mvec_ptr) V, const_map_ptr_t* map, int* ierr)
  {
  *ierr=-99;
  }

//! retrieve number of vectors/columns in V
void _SUBR_(mvec_num_vectors)(_TYPE_(const_mvec_ptr) V, int* nvec, int* ierr)
  {
  *ierr=-99;
  }

//!
void _SUBR_(mvec_extract_view)(_TYPE_(mvec_ptr) V, _ST_** val, lidx_t* lda, int* ierr)
  {
  *ierr=-99;
  }

//!
void _SUBR_(sdMat_extract_view)(_TYPE_(sdMat_ptr) V, _ST_** val, lidx_t* lda, int* ierr)
  {
  *ierr=-99;
  }

//! get a new vector that is a view of some columns of the original one,
//! Vblock = V(:,jmin:jmax). The new object Vblock is created but does not
//! allocate memory for the vector entries, instead using the entries from V
//! directly. When mvec_delete(Vblock) is called, the library has to take care
//! that the value array is not deleted 
void _SUBR_(mvec_view_block)(_TYPE_(mvec_ptr) V,
                             _TYPE_(mvec_ptr)* Vblock,
                             int jmin, int jmax, int* ierr)
  {
  *ierr=-99;
  }

//! get a new vector that is a copy of some columns of the original one,  
//! Vblock = V(:,jmin:jmax). The object Vblock must be created beforehand 
//! and the corresponding columns of V are copied into the value array    
//! of Vblock. V is not modified.
void _SUBR_(mvec_get_block)(_TYPE_(const_mvec_ptr) V,
                             _TYPE_(mvec_ptr) Vblock,
                             int jmin, int jmax, int* ierr)
  {
  *ierr=-99;
  }

//! given a multi-vector Vblock, set V(:,jmin:jmax)=Vblock by copying the corresponding
//! vectors. Vblock is not modified.
void _SUBR_(mvec_set_block)(_TYPE_(mvec_ptr) V,
                             _TYPE_(const_mvec_ptr) Vblock,
                             int jmin, int jmax, int* ierr)
  {
  *ierr=-99;
  }

//! get a new matrix that is a copy of some rows and columns of the original one,  
//! Mblock = M(imin:imax,jmin:jmax). The object Mblock must be created beforehand 
//! and the corresponding columns of M are copied into the value array    
//! of Mblock. M is not modified.
void _SUBR_(sdMat_get_block)(_TYPE_(const_mvec_ptr) M, 
                             _TYPE_(mvec_ptr) Mblock,
                             int imin, int imax, int jmin, int jmax, int* ierr)
  {
  *ierr=-99;
  }

//! given a serial dense matrix Mblock, set M(imin:imax,jmin:jmax)=Mblock by 
//! copying the corresponding elements. Mblock is not modified.
void _SUBR_(sdMat_set_block)(_TYPE_(sdMat_ptr) M, 
                             _TYPE_(const_sdMat_ptr) Mblock,
                             int imin, int imax, int jmin, int jmax, int* ierr)
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

//! normalize (in the 2-norm) each column of v and return ||v||_2
//! for each vector i in vnrm[i] (must be pre-allocated by caller)
void _SUBR_(mvec_normalize)(_TYPE_(mvec_ptr) V,
                            _MT_* vnrm, int* ierr)
  {
  *ierr=-99;
  }


//! scale each column i of v and by scalar[i]
void _SUBR_(mvec_scale)(_TYPE_(mvec_ptr) V, 
                        _ST_* scalar, int* ierr)
  {
  *ierr=-99;
  }

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
  
//! dense tall skinny matrix-matrix product yielding a serial dense matrix
//! C=alpha*V'*W+beta*C. C is replicated on all MPI processes sharing V and W.
void _SUBR_(mvecT_times_mvec)(_ST_ alpha, _TYPE_(const_mvec_ptr) V, 
                                       _TYPE_(const_mvec_ptr) W, 
                                       _ST_ beta, _TYPE_(sdMat_ptr) C, int* ierr)

  {
  *ierr=-99;
  }

//! 'tall skinny' QR decomposition, V=Q*R, Q'Q=I, R upper triangular.
//! Q is computed in place of V.
void _SUBR_(mvec_QR)(_TYPE_(mvec_ptr) V, _TYPE_(sdMat_ptr) R, int* ierr)
  {
  *ierr=-99;
  }

//!@}

