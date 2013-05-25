#include "phist_macros.h"

/*! these functions need to be provided by a kernel module  
    in order to be used by the iterative solvers implemen-  
    ted in this library. Matrices and vectors are passed    
    around as void pointers, which is not nice but it is    
    portable and interoperable with C and Fortran. 
    
    The numerical functions do not allocate memory, the user has
    to call the appropriate create functions beforehand (an 
    exception are the read functions for sparse matrix input and
    the functions to retrieve pointers to the map, see below).
    
    We currently request four types of objects. A 'map', defining
    the data distribution of a vector (this can be obtained from 
    a sparse matrix using the get_*_map functions). A sparse matrix 
    (crsMat), denoted by A below, typically distributed over several 
    MPI processes. A serial dense matrix (sdMat), which is replicated
    on all processes (not associated with an MPI_Comm). And
    a row-distributed vector, which may have several columns (mvec).

    Each subroutine has as its last argument the integer flag ierr, which
    should be 0 for success, positive for warnings and negative for errors.
    These standard error codes should be used:
    
    -99: not implemented.
    -88: cast of input args from void to required type failed.
    -77: caught an exception.
    -1: main mode of failure, for instance "memory allocation failed" in a "constructor".

*/


//! opaque pointer to multi-vector objects
typedef void* _TYPE_(mvec_ptr);

//! opaque pointer to const multi-vector objects
typedef const void* _TYPE_(const_mvec_ptr);

//! opaque pointer to serial dense matrix objects
typedef void* _TYPE_(sdMat_ptr);

//! opaque pointer to const serial dense matrix objects
typedef const void* _TYPE_(const_sdMat_ptr);

//! opaque pointer to CRS matrix objects
typedef void* _TYPE_(crsMat_ptr);

//! opaque pointer to const CRS matrix objects
typedef const void* _TYPE_(const_crsMat_ptr);

#ifdef __cplusplus
extern "C" {
#endif

//! returns 0 if the library implements the data type, -99 otherwise.
void _SUBR_(type_avail)(int* ierr);

//! \name Matrix input from a file

//@{

//! read a matrix from a MatrixMarket (ASCII) file
void _SUBR_(crsMat_read_mm)(_TYPE_(crsMat_ptr)* A, 
        const char* filename,int* ierr);

//! read a matrix from a Harwell-Boeing file
void _SUBR_(crsMat_read_hb)(_TYPE_(crsMat_ptr)* A, const char* filename,int* ierr);

//! read a matrix from a Ghost CRS (binary) file.
void _SUBR_(crsMat_read_bin)(_TYPE_(crsMat_ptr)* A, const char* filename,int* ierr);

//!@}

//! \name get information about the data distribution in a matrix (maps)

//!@{

//! get the row distribution of the matrix
void _SUBR_(crsMat_get_row_map)(_TYPE_(const_crsMat_ptr) A, 
        const_map_ptr_t* map, int* ierr);

//! get column distribution of a matrix
void _SUBR_(crsMat_get_col_map)(_TYPE_(const_crsMat_ptr) A, 
        const_map_ptr_t* map, int* ierr);

//! get the map for vectors x in y=A*x
void _SUBR_(crsMat_get_domain_map)(_TYPE_(const_crsMat_ptr) A, 
        const_map_ptr_t* map, int* ierr);

//! get the map for vectors y in y=A*x
void _SUBR_(crsMat_get_range_map)(_TYPE_(const_crsMat_ptr) A, 
        const_map_ptr_t* map, int* ierr);
//@}

//! \name constructors

//@{
//! create a block-vector.
void _SUBR_(mvec_create)(_TYPE_(mvec_ptr)* V, const_map_ptr_t map, int nvec, 
        int* ierr);

//! create a serial dense n x m matrix on all procs, with column major
//! ordering.
void _SUBR_(sdMat_create)(_TYPE_(sdMat_ptr)* M, 
        int nrows, int ncols, int* ierr);

//@}

//! \name destructors

//@{

//!
void _SUBR_(crsMat_delete)(_TYPE_(crsMat_ptr) A, int* ierr);

//!
void _SUBR_(mvec_delete)(_TYPE_(mvec_ptr) V, int* ierr);

//!
void _SUBR_(sdMat_delete)(_TYPE_(sdMat_ptr) M, int* ierr);

//@}

//! \name getting data from objects
//@{

//! retrieve local length of the vectors in V
void _SUBR_(mvec_my_length)(_TYPE_(const_mvec_ptr) V, int* len, int* ierr);

//! retrieve number of vectors/columns in V
void _SUBR_(mvec_num_vectors)(_TYPE_(const_mvec_ptr) V, int* nvec, int* ierr);

//! extract view from multi-vector
void _SUBR_(mvec_extract_view)(_TYPE_(mvec_ptr) V, _ST_** val,
        int* lda, int* ierr);

//! extract view from serial dense matrix
void _SUBR_(sdMat_extract_view)(_TYPE_(sdMat_ptr) M, _ST_** val,
        int* lda, int* ierr);

//@}

//! put scalar value into all elements of a multi-vector
void _SUBR_(mvec_put_value)(_TYPE_(mvec_ptr) V, _ST_ value, int* ierr);

//! put scalar value into all elements of a serial dense matrix
void _SUBR_(sdMat_put_value)(_TYPE_(mvec_ptr) V, _ST_ value, int* ierr);

//! put random numbers into all elements of a multi-vector
void _SUBR_(mvec_random)(_TYPE_(mvec_ptr) V, int* ierr);

//! put random numbers into all elements of a serial dense matrix
void _SUBR_(sdMat_random)(_TYPE_(sdMat_ptr) V, int* ierr);

//! \name Numerical functions
//!@{

//! y=alpha*x+beta*y
void _SUBR_(mvec_add_mvec)(_ST_ alpha, _TYPE_(const_mvec_ptr) X,
                            _ST_ beta,  _TYPE_(mvec_ptr)       Y,     
                            int* ierr);

//! y=alpha*A*x+beta*y. The scalars alpha and beta are expected to be of the
//! same type as the entries in the vectors and matrix. Mixing of types is
//! not allowed.
void _SUBR_(crsMat_times_mvec)(_ST_ alpha, _TYPE_(const_crsMat_ptr) A, 
        _TYPE_(const_mvec_ptr) x, _ST_ beta, _TYPE_(mvec_ptr) y, int* ierr);

//! dot product of vectors v_i and w_i, i=1..numvecs
void _SUBR_(mvec_dot_mvec)(_TYPE_(const_mvec_ptr) V, 
                            _TYPE_(const_mvec_ptr) W, 
                            _ST_* vw, int* ierr);

//! dense tall skinny matrix-matrix product yielding a serial dense matrix
//! C=alpha*V'*W+beta*C. C is replicated on all MPI processes sharing V and W.
void _SUBR_(mvecT_times_mvec)(_ST_ alpha, _TYPE_(const_mvec_ptr) V, 
                                       _TYPE_(const_mvec_ptr) W, 
                                       _ST_ beta, _TYPE_(sdMat_ptr) C, int* ierr);

//! n x m multi-vector times m x m dense matrix gives n x m multi-vector,
//! W=alpha*V*C + beta*W
void _SUBR_(mvec_times_sdMat)(_ST_ alpha, _TYPE_(const_mvec_ptr) V, 
                                       _TYPE_(const_sdMat_ptr) C,
                           _ST_ beta,  _TYPE_(mvec_ptr) W, 
                                       int* ierr);

//! 'tall skinny' QR decomposition, V=Q*R, Q'Q=I, R upper triangular.
void _SUBR_(mvec_QR)(_TYPE_(const_mvec_ptr) V, 
                            _TYPE_(mvec_ptr) Q, 
                            _TYPE_(sdMat_ptr) R, int* ierr);

//!@}


#ifdef __cplusplus
} //extern "C"
#endif
