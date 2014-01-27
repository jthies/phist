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

#ifdef __cplusplus
extern "C" {
#endif

//! returns 0 if the library implements the data type, -99 otherwise.
void SUBR(type_avail)(int* ierr);

//! \name Matrix input from a file

//@{

//! read a matrix from a MatrixMarket (ASCII) file
void SUBR(crsMat_read_mm)(TYPE(crsMat_ptr)* A, 
        const char* filename,int* ierr);

//! read a matrix from a Harwell-Boeing file
void SUBR(crsMat_read_hb)(TYPE(crsMat_ptr)* A, const char* filename,int* ierr);

//! read a matrix from a Ghost CRS (binary) file.
void SUBR(crsMat_read_bin)(TYPE(crsMat_ptr)* A, const char* filename,int* ierr);

//!@}

//! \name get information about the data distribution in a matrix (maps)

//!@{

//! get the row distribution of the matrix
void SUBR(crsMat_get_row_map)(TYPE(const_crsMat_ptr) A, 
        const_map_ptr_t* map, int* ierr);

//! get column distribution of a matrix
void SUBR(crsMat_get_col_map)(TYPE(const_crsMat_ptr) A, 
        const_map_ptr_t* map, int* ierr);

//! get the map for vectors x in y=A*x
void SUBR(crsMat_get_domain_map)(TYPE(const_crsMat_ptr) A, 
        const_map_ptr_t* map, int* ierr);

//! get the map for vectors y in y=A*x
void SUBR(crsMat_get_range_map)(TYPE(const_crsMat_ptr) A,
        const_map_ptr_t* map, int* ierr);
//@}

//! \name constructors

//@{
//! create a block-vector.
void SUBR(mvec_create)(TYPE(mvec_ptr)* V, const_map_ptr_t map, int nvec, 
        int* ierr);

//! create a block-vector as view of raw data. The map tells the object
//! how many rows it should 'see' in the data (at most lda, the leading
//! dimension of the 2D array values).
void SUBR(mvec_create_view)(TYPE(mvec_ptr)* V, const_map_ptr_t map, 
        _ST_* values, lidx_t lda, int nvec, 
        int* ierr);

//! create a serial dense n x m matrix on all procs in comm, 
//! with column major ordering and the capability to communicate.
//! TODO - ghost can use int64_t as lidx_t, should the sdMat have
//! long ints as well? It should be possible to view the local
//! part of a vector as sdMat, so at least nrows should become
//! lidx_t, I think.
void SUBR(sdMat_create)(TYPE(sdMat_ptr)* M, 
        int nrows, int ncols, const_comm_ptr_t comm, int* ierr);

//@}

//! \name destructors

//@{

//!
void SUBR(crsMat_delete)(TYPE(crsMat_ptr) A, int* ierr);

//!
void SUBR(mvec_delete)(TYPE(mvec_ptr) V, int* ierr);

//!
void SUBR(sdMat_delete)(TYPE(sdMat_ptr) M, int* ierr);

//@}

//! \name getting data from objects
//@{

//! retrieve local length of the vectors in V
void SUBR(mvec_my_length)(TYPE(const_mvec_ptr) V, lidx_t* len, int* ierr);

//! retrieve the map of the vectors in V
void SUBR(mvec_get_map)(TYPE(const_mvec_ptr) V, const_map_ptr_t* map, int* ierr);

//! retrieve the comm used for MPI communication in V
void SUBR(mvec_get_comm)(TYPE(const_mvec_ptr) V, const_comm_ptr_t* comm, int* ierr);

//! retrieve number of vectors/columns in V
void SUBR(mvec_num_vectors)(TYPE(const_mvec_ptr) V, int* nvec, int* ierr);

//! get number of cols in local dense matrix
void SUBR(sdMat_get_nrows)(TYPE(const_sdMat_ptr) M, int* nrows, int* ierr);

//! get number of cols in local dense matrix
void SUBR(sdMat_get_ncols)(TYPE(const_sdMat_ptr) M, int* ncols, int* ierr);

//! extract view from multi-vector. Sets the user-provided val pointer to point to the
//! beginning of the first vector, and puts the leading dimension of the array into lda,
//! such that the first element of vector i is val[i*lda]. Note that lda may be larger
//! than the actual local vector length as obtained by mvec_my_length. This function is
//! dangerous in the sense that it would force the underlying kernel lib to implement the
//! data layout in this given fashion. A library that does not guarantee this should return
//! -99 here ("not implemented")
void SUBR(mvec_extract_view)(TYPE(mvec_ptr) V, _ST_** val,
        lidx_t* lda, int* ierr);

//! extract view from serial dense matrix. See comment for mvec_extract_view for details.
void SUBR(sdMat_extract_view)(TYPE(sdMat_ptr) M, _ST_** val,
        lidx_t* lda, int* ierr);

//@}

//! get a new vector that is a view of some columns of the original one,
//! Vblock = V(:,jmin:jmax). The new object Vblock is created but does not
//! allocate memory for the vector entries, instead using the entries from V
//! directly. When mvec_delete(Vblock) is called, the library has to take care
//! that the value array is not deleted.
//! If on entry, *Vblock!=NULL, this function should delete *Vblock so that
//! repeated use of view_block does not lead to memory holes. It is crucial
//! that you pass in a NULL pointer if you want a new object, otherwise you
//! may get a segfault.
//! 
//! The user is responsible for deleting the objects in the correct order: First
//! the view of V, then V itself. It is allowed to create a view of a view, but 
//! then again, the order of deletion has to be observed by the user. For instance,
//! this code may run into trouble:
//!
//!     Dmvec_ptr_t A, Av, Avv;
//!     phist_Dmvec_create(&A,...);
//!     phist_Dmvec_view_block(A,&Av,...);
//!     phist_Dmvec_view_block(Av,&Avv,...);
//!     phist_Dmvec_delete(Av,...);
//! (do something with Avv)
void SUBR(mvec_view_block)(TYPE(mvec_ptr) V, 
                             TYPE(mvec_ptr)* Vblock,
                             int jmin, int jmax, int* ierr);

//! get a new vector that is a copy of some columns of the original one,  
//! Vblock = V(:,jmin:jmax). The object Vblock must be created beforehand 
//! and the corresponding columns of V are copied into the value array    
//! of Vblock. V is not modified.
void SUBR(mvec_get_block)(TYPE(const_mvec_ptr) V, 
                             TYPE(mvec_ptr) Vblock,
                             int jmin, int jmax, int* ierr);

//! given a multi-vector Vblock, set V(:,jmin:jmax)=Vblock by copying the corresponding
//! vectors. Vblock is not modified.
void SUBR(mvec_set_block)(TYPE(mvec_ptr) V, 
                             TYPE(const_mvec_ptr) Vblock,
                             int jmin, int jmax, int* ierr);

//! get a new matrix that is a view of some rows and columns of the original one, 
//! Mblock = M(imin:imax,jmin:jmax). The behavior is analogous to mvec_view_block.
//! If on entry, *Mblock!=NULL, this function should delete *Mblock so that
//! repeated use of view_block does not lead to memory holes. It is crucial
//! that you pass in a NULL pointer if you want a new object, otherwise you
//! may get a segfault.
void SUBR(sdMat_view_block)(TYPE(sdMat_ptr) M, 
                             TYPE(sdMat_ptr)* Mblock,
                             int imin, int imax, int jmin, int jmax, int* ierr);

//! get a new matrix that is a copy of some rows and columns of the original one,  
//! Mblock = M(imin:imax,jmin:jmax). The object Mblock must be created beforehand 
//! and the corresponding columns of M are copied into the value array    
//! of Mblock. M is not modified.
void SUBR(sdMat_get_block)(TYPE(const_sdMat_ptr) M, 
                             TYPE(sdMat_ptr) Mblock,
                             int imin, int imax, int jmin, int jmax, int* ierr);

//! given a serial dense matrix Mblock, set M(imin:imax,jmin:jmax)=Mblock by 
//! copying the corresponding elements. Mblock is not modified.
void SUBR(sdMat_set_block)(TYPE(sdMat_ptr) M, 
                             TYPE(const_sdMat_ptr) Mblock,
                             int imin, int imax, int jmin, int jmax, int* ierr);

//! put scalar value into all elements of a multi-vector
void SUBR(mvec_put_value)(TYPE(mvec_ptr) V, _ST_ value, int* ierr);

//! put scalar value into all elements of a serial dense matrix
void SUBR(sdMat_put_value)(TYPE(sdMat_ptr) V, _ST_ value, int* ierr);

//! put random numbers into all elements of a multi-vector
void SUBR(mvec_random)(TYPE(mvec_ptr) V, int* ierr);

//! put random numbers into all elements of a serial dense matrix
void SUBR(sdMat_random)(TYPE(sdMat_ptr) V, int* ierr);

//! print a vector to the screen (for debugging)
void SUBR(mvec_print)(TYPE(const_mvec_ptr) V, int* ierr);

//! print an sdMat to the screen (for debugging)
void SUBR(sdMat_print)(TYPE(const_sdMat_ptr) M, int* ierr);

//! \name Numerical functions
//!@{

//! compute the 2-norm) of each column of v
//! (vnrm[i] must be pre-allocated by caller)
void SUBR(mvec_norm2)(TYPE(const_mvec_ptr) V, 
                            _MT_* vnrm, int* ierr);

//! normalize (in the 2-norm) each column of v and return ||v||_2
//! for each vector i in vnrm[i] (must be pre-allocated by caller)
void SUBR(mvec_normalize)(TYPE(mvec_ptr) V, 
                            _MT_* vnrm, int* ierr);

//! scale each column i of v and by scalar
void SUBR(mvec_scale)(TYPE(mvec_ptr) V, 
                            _ST_ scalar, int* ierr);

//! scale each column i of v and by scalar[i]
void SUBR(mvec_vscale)(TYPE(mvec_ptr) V, 
                            const _ST_* scalar, int* ierr);

//! y=alpha*x+beta*y.
//! This function can also be used for special cases such as
//! alpha=0 => scale y
//! alpha=1, beta=0 => copy y=x
//! alpha!=0, beta=1: 'axpy' operation
void SUBR(mvec_add_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) X,
                            _ST_ beta,  TYPE(mvec_ptr)       Y,     
                            int* ierr);

//! y[i]=alpha[i]*x[i]+beta*y[i]
void SUBR(mvec_vadd_mvec)(const _ST_ alpha[], TYPE(const_mvec_ptr) X,
                          const _ST_ beta,  TYPE(mvec_ptr)       Y,     
                            int* ierr);

//! B=alpha*A+beta*B
void SUBR(sdMat_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) A,
                            _ST_ beta,  TYPE(sdMat_ptr)       B,     
                            int* ierr);

//! y=alpha*A*x+beta*y. The scalars alpha and beta are expected to be of the
//! same type as the entries in the vectors and matrix. Mixing of types is
//! not allowed.
void SUBR(crsMat_times_mvec)(_ST_ alpha, TYPE(const_crsMat_ptr) A, 
        TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* ierr);

//! y=alpha*A^H*x+beta*y. The scalars alpha and beta are expected to be of the
//! same type as the entries in the vectors and matrix. Mixing of types is
//! not allowed. In the complex case, the conjugate transpose is used.
void SUBR(crsMatT_times_mvec)(_ST_ alpha, TYPE(const_crsMat_ptr) A, 
        TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* ierr);

//! y[i]=alpha*(A*x[i]+shifts[i]*x[i]) + beta*y[i]
void SUBR(crsMat_times_mvec_vadd_mvec)(_ST_ alpha, TYPE(const_crsMat_ptr) A,
        const _ST_ shifts[], TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* ierr);

//! dot product of vectors v_i and w_i, i=1..numvecs
void SUBR(mvec_dot_mvec)(TYPE(const_mvec_ptr) V, 
                            TYPE(const_mvec_ptr) W, 
                            _ST_* vw, int* ierr);

//! dense tall skinny matrix-matrix product yielding a serial dense matrix
//! C=alpha*V'*W+beta*C. C is replicated on all MPI processes sharing V and W.
void SUBR(mvecT_times_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) V, 
                                       TYPE(const_mvec_ptr) W, 
                                       _ST_ beta, TYPE(sdMat_ptr) C, int* ierr);

//! n x m multi-vector times m x k dense matrix gives n x k multi-vector,
//! W=alpha*V*C + beta*W
void SUBR(mvec_times_sdMat)(_ST_ alpha, TYPE(const_mvec_ptr) V, 
                                       TYPE(const_sdMat_ptr) C,
                           _ST_ beta,  TYPE(mvec_ptr) W, 
                                       int* ierr);

//! n x m serial dense matrix times m x k serial dense matrix gives n x k serial dense matrix,
//! C=alpha*V*W + beta*C
void SUBR(sdMat_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) V, 
                                           TYPE(const_sdMat_ptr) W, 
                               _ST_ beta, TYPE(sdMat_ptr) C,
                                       int* ierr);

//! n x m conj. transposed serial dense matrix times m x k serial dense matrix gives m x k serial dense matrix,
//! C=alpha*V*W + beta*C
void SUBR(sdMatT_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) V, 
                                           TYPE(const_sdMat_ptr) W, 
                               _ST_ beta, TYPE(sdMat_ptr) C,
                                       int* ierr);


//! 'tall skinny' QR decomposition, V=Q*R, Q'Q=I, R upper triangular.   
//! Q is computed in place of V. If V does not have full rank, ierr>0   
//! indicates the dimension of the null-space of V. The first m-ierr    
//! columns of Q are an orthogonal basis of the column space of V, the  
//! remaining columns form a basis for the null space.                  
//!                                                                     
//! TODO: our current three implementations of this function (ghost,    
//!     tpetra and epetra) are identical, it would be nicer to move this
//!     function into core/ so we have a single implementation.         
void SUBR(mvec_QR)(TYPE(mvec_ptr) V, 
                     TYPE(sdMat_ptr) R, int* ierr);


#ifdef PHIST_KERNEL_LIB_FORTRAN
void SUBR(crsMat_create_fromRowFunc)(TYPE(crsMat_ptr) *A, int nrows, int ncols, int maxnne, void (*rowFunPtr)(int32_t,int32_t*,int32_t*,void*), int *ierr);
void SUBR(mvec_gather_mvecs)(TYPE(mvec_ptr) V, TYPE(const_mvec_ptr) W[], int nblocks, int *ierr);
void SUBR(mvec_scatter_mvecs)(TYPE(const_mvec_ptr) V, TYPE(mvec_ptr) W[], int nblocks, int *ierr);
void SUBR(mvec_times_sdMat_inplace)(TYPE(mvec_ptr) V, TYPE(const_sdMat_ptr) M, int *ierr);
#endif


//!@}



#ifdef __cplusplus
} //extern "C"
#endif
