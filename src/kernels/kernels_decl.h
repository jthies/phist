
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

//! \name Matrix input from a file

//@{

//! read a matrix from a MatrixMarket (ASCII) file
_SUBROUTINE_(read_crsMat_mm)(void** A, char* filename,int* ierr);

//! read a matrix from a Harwell-Boeing file
_SUBROUTINE_(read_crsMat_hb)(void** A, char* filename,int* ierr);

//! read a matrix from a Ghost CRS (binary) file.
_SUBROUTINE_(read_crsMat_bin)(void** A, char* filename,int* ierr);

//!@}

//! \name get information about the data distribution in a matrix (maps)

//!@{
//! get the row distribution of the matrix
_SUBROUTINE_(get_row_map)(void* A, void** map, int* ierr);

//! get column distribution of a matrix
_SUBROUTINE_(get_col_map)(void* A, void** map, int* ierr);

//! get the map for vectors x in y=A*x
_SUBROUTINE_(get_domain_map)(void* A, void** map, int* ierr);

//! get the map for vectors y in y=A*x
_SUBROUTINE_(get_range_map)(void* A, void** map, int* ierr);
//@}

//! \name constructors

//@{
//! create a block-vector. The entries are stored contiguously
//! at val in column major ordering.
_SUBROUTINE_(create_mvec)(void* map, int nvec, void** V, _ST_** val, int* ierr);

//! create a serial dense n x m matrix on all procs, with column major
//! ordering.
_SUBROUTINE_(create_sdMat)(int nrows, int ncols, void** M, _ST_** val, int* ierr);

//@}

//! \name destructors

//@{

//!
_SUBROUTINE_(delete_crsMat)(void* A, int* ierr);

//!
_SUBROUTINE_(delete_mvec)(void* v, int* ierr);

//!
_SUBROUTINE_(delete_sdMat)(void* M, int* ierr);

//@}

//! \name Numerical functions
//!@{

//! y=alpha*A*x+beta*y. The scalars alpha and beta are expected to be of the
//! same type as the entries in the vectors and matrix. Mixing of types is
//! not allowed.
_SUBROUTINE_(crsMat_X_mvec)(_ST_ alpha, void* A, void* x, _ST_ beta, void* y, int* ierr);

//! dot product of vectors v_i and w_i, i=1..numvecs
_SUBROUTINE_(mvec_dot_mvec)(void* v, void* w, _ST_* vw, int* ierr);

//! dense tall skinny matrix-matrix product yielding a serial dense matrix
//! C=alpha*V'*W+beta*C. C is replicated on all MPI processes sharing V and W.
_SUBROUTINE_(mvecT_X_mvec)(_ST_ alpha, void* V, void* W, _ST_ beta, void* C, int* ierr);

//! 'tall skinny' QR decomposition, V=Q*R, Q'Q=I, R upper triangular.
_SUBROUTINE_(mvec_QR)(void* V, void* Q, void* R, int* ierr);

//!@}

