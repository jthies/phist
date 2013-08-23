#include "phist_kernels.h"
#include "phist_macros.h"

// these implementations can be used if a kernel package
// does not implement all four data types (cf. epetra/ for
// an example, which supports only double (D))

void SUBR(type_avail)(int *ierr)
  {
  *ierr=-99;
  }

// \name Matrix input from a file

//@{

//! read a matrix from a MatrixMarket (ASCII) file
void SUBR(crsMat_read_mm)(TYPE(crsMat_ptr)* A, const char* filename,int* ierr)
  {
  *ierr=-99;
  }

//! read a matrix from a Ghost CRS (binary) file.
void SUBR(crsMat_read_bin)(TYPE(crsMat_ptr)* A, const char* filename,int* ierr)
  {
  *ierr=-99;
  }

//! read a matrix from a Harwell-Boeing (HB) file
void SUBR(crsMat_read_hb)(TYPE(crsMat_ptr)* A, const char* filename,int* ierr)
  {
  *ierr=-99;
  }
//!@}

//! \name get information about the data distribution in a matrix (maps)

//!@{
//! get the row distribution of the matrix
void SUBR(crsMat_get_row_map)(TYPE(const_crsMat_ptr) A, const_map_ptr_t* map, int* ierr)
  {
  *ierr=-99;
  }

//! get column distribution of a matrix
void SUBR(crsMat_get_col_map)(TYPE(const_crsMat_ptr) A, const_map_ptr_t* map, int* ierr)
  {
  *ierr=-99;
  }

//! get the map for vectors x in y=A*x
void SUBR(crsMat_get_domain_map)(TYPE(const_crsMat_ptr) A, const_map_ptr_t* map, int* ierr)
  {
  *ierr=-99;
  }

//! get the map for vectors y in y=A*x
void SUBR(crsMat_get_range_map)(TYPE(const_crsMat_ptr) A, const_map_ptr_t* map, int* ierr)
  {
  *ierr=-99;
  }
//@}

//! \name constructors

//@{
//! create a block-vector. The entries are stored contiguously
//! at val in column major ordering.
void SUBR(mvec_create)(TYPE(mvec_ptr)* V, 
        const_map_ptr_t map, lidx_t nvec, int* ierr)
  {
  *ierr=-99;
  }

//! create a block-vector as view of raw data. The map tells the object
//! how many rows it should 'see' in the data (at most lda, the leading
//! dimension of the 2D array values).
void SUBR(mvec_create_view)(TYPE(mvec_ptr)* V, const_map_ptr_t map, 
        _ST_* values, lidx_t lda, int nvec,
        int* ierr)
  {
  *ierr=-99;
  }

//! create a serial dense n x m matrix on all procs in comm, with column major
//! ordering and the capability to communicate.
void SUBR(sdMat_create)(TYPE(sdMat_ptr)* M, 
int nrows, int ncols, const_comm_ptr_t comm, int* ierr)
  {
  *ierr=-99;
  }

//@}

//! retrieve local length of the vectors in V
void SUBR(mvec_my_length)(TYPE(const_mvec_ptr) V, lidx_t* len, int* ierr)
  {
  *ierr=-99;
  }

//! retrieve the map of the vectors in V
void SUBR(mvec_get_map)(TYPE(const_mvec_ptr) V, const_map_ptr_t* map, int* ierr)
  {
  *ierr=-99;
  }

//! retrieve the comm used for MPI communication in V
void SUBR(mvec_get_comm)(TYPE(const_mvec_ptr) V, const_comm_ptr_t* comm, int* ierr)
  {
  *ierr=-99;
  }
//! retrieve number of vectors/columns in V
void SUBR(mvec_num_vectors)(TYPE(const_mvec_ptr) V, int* nvec, int* ierr)
  {
  *ierr=-99;
  }

//! get number of cols in local dense matrix
void SUBR(sdMat_get_nrows)(TYPE(const_sdMat_ptr) M, int* nrows, int* ierr)
  {
  *ierr=-99;
  }

//! get number of cols in local dense matrix
void SUBR(sdMat_get_ncols)(TYPE(const_sdMat_ptr) M, int* ncols, int* ierr)
  {
  *ierr=-99;
  }


//!
void SUBR(mvec_extract_view)(TYPE(mvec_ptr) V, _ST_** val, lidx_t* lda, int* ierr)
  {
  *ierr=-99;
  }

//!
void SUBR(sdMat_extract_view)(TYPE(sdMat_ptr) V, _ST_** val, lidx_t* lda, int* ierr)
  {
  *ierr=-99;
  }

//! get a new vector that is a view of some columns of the original one,
//! Vblock = V(:,jmin:jmax). The new object Vblock is created but does not
//! allocate memory for the vector entries, instead using the entries from V
//! directly. When mvec_delete(Vblock) is called, the library has to take care
//! that the value array is not deleted 
void SUBR(mvec_view_block)(TYPE(mvec_ptr) V,
                             TYPE(mvec_ptr)* Vblock,
                             int jmin, int jmax, int* ierr)
  {
  *ierr=-99;
  }

//! get a new vector that is a copy of some columns of the original one,  
//! Vblock = V(:,jmin:jmax). The object Vblock must be created beforehand 
//! and the corresponding columns of V are copied into the value array    
//! of Vblock. V is not modified.
void SUBR(mvec_get_block)(TYPE(const_mvec_ptr) V,
                             TYPE(mvec_ptr) Vblock,
                             int jmin, int jmax, int* ierr)
  {
  *ierr=-99;
  }

//! given a multi-vector Vblock, set V(:,jmin:jmax)=Vblock by copying the corresponding
//! vectors. Vblock is not modified.
void SUBR(mvec_set_block)(TYPE(mvec_ptr) V,
                             TYPE(const_mvec_ptr) Vblock,
                             int jmin, int jmax, int* ierr)
  {
  *ierr=-99;
  }

//! get a new matrix that is a copy of some rows and columns of the original one,  
//! Mblock = M(imin:imax,jmin:jmax). The object Mblock must be created beforehand 
//! and the corresponding columns of M are copied into the value array    
//! of Mblock. M is not modified.
void SUBR(sdMat_get_block)(TYPE(const_mvec_ptr) M, 
                             TYPE(mvec_ptr) Mblock,
                             int imin, int imax, int jmin, int jmax, int* ierr)
  {
  *ierr=-99;
  }

//! given a serial dense matrix Mblock, set M(imin:imax,jmin:jmax)=Mblock by 
//! copying the corresponding elements. Mblock is not modified.
void SUBR(sdMat_set_block)(TYPE(sdMat_ptr) M, 
                             TYPE(const_sdMat_ptr) Mblock,
                             int imin, int imax, int jmin, int jmax, int* ierr)
  {
  *ierr=-99;
  }

//! \name destructors

//@{

//!
void SUBR(crsMat_delete)(TYPE(crsMat_ptr) A, int* ierr)
  {
  *ierr=-99;
  }

//!
void SUBR(mvec_delete)(TYPE(mvec_ptr) V, int* ierr)
  {
  *ierr=-99;
  }

//!
void SUBR(sdMat_delete)(TYPE(sdMat_ptr) M, int* ierr)
  {
  *ierr=-99;
  }

//@}

//! put scalar value into all elements of a multi-vector
void SUBR(mvec_put_value)(TYPE(mvec_ptr) V, _ST_ value, int* ierr)
  {
  *ierr=-99;
  }

//! put scalar value into all elements of a serial dense matrix
void SUBR(sdMat_put_value)(TYPE(mvec_ptr) V, _ST_ value, int* ierr)
  {
  *ierr=-99;
  }

//! put random numbers into all elements of a multi-vector
void SUBR(mvec_random)(TYPE(mvec_ptr) V, int* ierr)
  {
  *ierr=-99;
  }

//! put random numbers into all elements of a serial dense matrix
void SUBR(sdMat_random)(TYPE(sdMat_ptr) M, int* ierr)
  {
  *ierr=-99;
  }


//! \name Numerical functions
//!@{

//! compute the 2-norm) of each column of v                   
//! (vnrm[i] must be pre-allocated by caller)
void SUBR(mvec_norm2)(TYPE(const_mvec_ptr) V,
                            _MT_* vnrm, int* ierr)
{
*ierr=-99;
}

//! normalize (in the 2-norm) each column of v and return ||v||_2
//! for each vector i in vnrm[i] (must be pre-allocated by caller)
void SUBR(mvec_normalize)(TYPE(mvec_ptr) V,
                            _MT_* vnrm, int* ierr)
  {
  *ierr=-99;
  }


//! scale each column i of v and by scalar
void SUBR(mvec_scale)(TYPE(mvec_ptr) V, 
                        _ST_ scalar, int* ierr)
  {
  *ierr=-99;
  }

//! scale each column i of v and by scalar[i]
void SUBR(mvec_vscale)(TYPE(mvec_ptr) V, 
                        _ST_* scalar, int* ierr)
  {
  *ierr=-99;
  }

//! y=alpha*x+beta*y
void SUBR(mvec_add_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) X,
                            _ST_ beta,  TYPE(mvec_ptr)       Y, 
                            int* ierr)
  {
  *ierr=-99;
  }

//! B=alpha*A+beta*B
void SUBR(sdMat_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) A,
                            _ST_ beta,  TYPE(sdMat_ptr)       B, 
                            int* ierr)
  {
  *ierr=-99;
  }

//! y=alpha*A*x+beta*y.
void SUBR(crsMat_times_mvec)(_ST_ alpha, TYPE(const_crsMat_ptr) A, 
        TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* ierr)
  {
  *ierr=-99;
  }

//! dot product of vectors v_i and w_i, i=1..numvecs
void SUBR(mvec_dot_mvec)(TYPE(const_mvec_ptr) v, 
                            TYPE(const_mvec_ptr) w, 
                            _ST_* s, int* ierr)
  {
  *ierr=-99;
  }
  
//! W=V*C
void SUBR(mvec_times_sdMat)(_ST_ alpha, TYPE(const_mvec_ptr) V, 
                                       TYPE(const_sdMat_ptr) C, 
                                       _ST_ beta, TYPE(mvec_ptr) W, int* ierr)
  {
  *ierr=-99;
  }

//! C=V*W
void SUBR(sdMat_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) V, 
                                       TYPE(const_sdMat_ptr) V, 
                                       _ST_ beta, TYPE(sdMat_ptr) W, int* ierr)
  {
  *ierr=-99;
  }

//! dense tall skinny matrix-matrix product yielding a serial dense matrix
//! C=alpha*V'*W+beta*C. C is replicated on all MPI processes sharing V and W.
void SUBR(mvecT_times_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) V, 
                                       TYPE(const_mvec_ptr) W, 
                                       _ST_ beta, TYPE(sdMat_ptr) C, int* ierr)
  {
  *ierr=-99;
  }

//! 'tall skinny' QR decomposition, V=Q*R, Q'Q=I, R upper triangular.
//! Q is computed in place of V.
void SUBR(mvec_QR)(TYPE(mvec_ptr) V, TYPE(sdMat_ptr) R, int* ierr)
  {
  *ierr=-99;
  }

//!@}

