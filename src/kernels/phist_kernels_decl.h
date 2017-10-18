/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_macros.h"

//! \addtogroup kernels
//@{

/*! @file phist_kernels_decl.h 
    @brief basic operations involving matrices and vectors 

    these functions need to be provided by a kernel module  
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
    (sparseMat), denoted by A below, typically distributed over several 
    MPI processes. A small dense matrix (sdMat), which is replicated
    on all processes (not associated with an MPI_Comm). And
    a row-distributed vector, which may have several columns (mvec).

    STATUS FLAG
    -----------

    Each subroutine has as its last argument the integer flag iflag, which
    should be 0 for success, positive for warnings and negative for errors.
    These standard error codes should be used:
    
    -99: not implemented.
    -88: cast of input args from void to required type failed.
    -77: caught an exception.

    etc., for a full list of defined return codes see phist_defs.h

    -1..-9: function specific modes of failure, which will have to be looked up
            in the source code (documentation) right now.

    OPTIONAL ARGUMENTS
    ------------------
    
    Some of the kernel functions allow passing in flags via *iflag to provide the kernel 
    library with information about the situation in which the function is called, see the 
    documentation of the individual functions for usage details.
    
    NOTE: The kernel library is *not* required to honour the flags given. For instance, if 
          you tell a function to read a matrix from a file and repartition it for 
          minimizing communication, it may or may not do so, depending on wether TPLs are 
          available or the kernel lib implements the feature at all.
*/

#ifdef __cplusplus
extern "C" {
#endif


//! returns 0 if the library implements the data type, -99 otherwise.
void SUBR(type_avail)(int* iflag);

//!   \defgroup crsmat Sparse matrix functions (sparseMat_t) 
//@{
//! \name Matrix input from a file
//!
//! optional flags:
//! * PHIST_SPARSEMAT_PERM_GLOBAL
//! * PHIST_SPARSEMAT_DIST2_COLOR (feature required for CARP-CG)
//!
//! If a predefined map is used, the permutation and partitioning is defined solely by the map.
//!
///@{

//! read a matrix from a MatrixMarket (ASCII) file \ingroup(crsmat)
void SUBR(sparseMat_read_mm)(TYPE(sparseMat_ptr)* A, phist_const_comm_ptr comm,
        const char* filename,int* iflag);

//! read a matrix from a Harwell-Boeing file \ingroup(crsmat)
void SUBR(sparseMat_read_hb)(TYPE(sparseMat_ptr)* A, phist_const_comm_ptr comm,
        const char* filename,int* iflag);

//! read a matrix from a Ghost CRS (binary) file \ingroup(crsmat)
void SUBR(sparseMat_read_bin)(TYPE(sparseMat_ptr)* A, phist_const_comm_ptr comm,
const char* filename,int* iflag);

//! read a matrix from a MatrixMarket (ASCII) file with predefined context \ingroup(crsmat)
void SUBR(sparseMat_read_mm_with_context)(TYPE(sparseMat_ptr)* A, phist_const_context_ptr ctx,
        const char* filename,int* iflag);

//! read a matrix from a Harwell-Boeing file with predefined context \ingroup(crsmat)
void SUBR(sparseMat_read_hb_with_context)(TYPE(sparseMat_ptr)* A, phist_const_context_ptr ctx,
        const char* filename,int* iflag);

//! read a matrix from a Ghost CRS (binary) file with predefined context \ingroup(crsmat)
void SUBR(sparseMat_read_bin_with_context)(TYPE(sparseMat_ptr)* A, phist_const_context_ptr ctx,
const char* filename,int* iflag);

///@}

//! \name get information about the data distribution in a matrix (maps)
///@{

//! get the row distribution of the matrix
void SUBR(sparseMat_get_row_map)(TYPE(const_sparseMat_ptr) A, 
        phist_const_map_ptr* map, int* iflag);

//! get column distribution of a matrix
void SUBR(sparseMat_get_col_map)(TYPE(const_sparseMat_ptr) A, 
        phist_const_map_ptr* map, int* iflag);

//! get the map for vectors x in y=A*x
void SUBR(sparseMat_get_domain_map)(TYPE(const_sparseMat_ptr) A, 
        phist_const_map_ptr* map, int* iflag);

//! get the map for vectors y in y=A*x
void SUBR(sparseMat_get_range_map)(TYPE(const_sparseMat_ptr) A,
        phist_const_map_ptr* map, int* iflag);

//! get the complete context for creating a similar matrix
void SUBR(sparseMat_get_context)(TYPE(const_sparseMat_ptr) A,
        phist_const_context_ptr* ctx, int* iflag);
///@}
//@}

//! \name constructors
//@{

//! create a block-vector. \ingroup mvec
void SUBR(mvec_create)(TYPE(mvec_ptr)* V, phist_const_map_ptr map, int nvec, 
        int* iflag);

//! construct small dense matrix \ingroup sdmat

//! create a small dense n x m matrix on all procs in comm,   
//! with column major ordering.
//! If comm!=NULL, the object has the capability to communicate
//! and can be used in functions like mvecT_times_mvec if the  
//! map of the mvecs uses the same comm. Otherwise, it is a lo-
//! cal object (MPI_COMM_SELF is assumed).
void SUBR(sdMat_create)(TYPE(sdMat_ptr)* M, 
        int nrows, int ncols, phist_const_comm_ptr comm, int* iflag);

//@}

//! \name destructors
//@{

//! delete sparseMat \ingroup crsmat
void SUBR(sparseMat_delete)(TYPE(sparseMat_ptr) A, int* iflag);

//! delete mvec \ingroup mvec
void SUBR(mvec_delete)(TYPE(mvec_ptr) V, int* iflag);

//! delete sdMat \ingroup sdmat
void SUBR(sdMat_delete)(TYPE(sdMat_ptr) M, int* iflag);

//@}

//! \name getting data from objects
//@{

//! retrieve the map of the vectors in V \ingroup mvec
void SUBR(mvec_get_map)(TYPE(const_mvec_ptr) V, phist_const_map_ptr* map, int* iflag);

//! retrieve number of vectors/columns in V \ingroup mvec
void SUBR(mvec_num_vectors)(TYPE(const_mvec_ptr) V, int* nvec, int* iflag);

//! extract view from multi-vector. \ingroup mvec

//! Sets the user-provided val pointer to point to the
//! beginning of the first vector, and puts the leading dimension of the array into lda,
//! such that the element j of vector i is val[i*lda+j]. If PHIST_MVECS_ROW_MAJOR is
//! #defined, the storage is transposed and element j of vector i is found at val[j*lda+i].
//! Note that lda may be larger than the actual local vector length as obtained by 
//! mvec_my_length. This function is dangerous in the sense that it would force the 
//! underlying kernel lib to implement the data layout in either of these formats.
//! library that does not guarantee this should return -99 here ("not implemented").
//!
//! Note on accelerators: the pointer obtained lives on the CPU, in order to make sure
//!     the data is consistent with the copy on a GPU, always call mvec_from_device before
//!     accessing the raw data and mvec_upload after manipulating it. Make sure that you
//!     know which copy is up-to-date at what point.
//!
void SUBR(mvec_extract_view)(TYPE(mvec_ptr) V, _ST_** V_raw,
        phist_lidx* lda, int* iflag);

//! get number of cols in local dense matrix \ingroup sdmat
void SUBR(sdMat_get_nrows)(TYPE(const_sdMat_ptr) M, int* nrows, int* iflag);

//! get number of cols in local dense matrix. \ingroup sdmat
void SUBR(sdMat_get_ncols)(TYPE(const_sdMat_ptr) M, int* ncols, int* iflag);

//! extract view from small dense matrix. \ingroup sdmat

//! See comment for mvec_extract_view for details,
//! the macro indicating row-major storage layout is PHIST_SDMATS_ROW_MAJOR.
void SUBR(sdMat_extract_view)(TYPE(sdMat_ptr) M, _ST_** M_raw,
        phist_lidx* lda, int* iflag);

#ifdef PHIST_HIGH_PRECISION_KERNELS
//! extract pointer to least significant bits of double-double matrix elements,
//! this obviously makes some assumptions on how the high precision arithmetic is
//! implemented by the kernel lib, but as we will only support these features with
//! builtin and ghost kernels for now, it's the easiest way of providing some
//! common implemetations for e.g. lapack-like routines in higher precision.
//!
//! NOTE: the array is assumed to have the same storage layout as the one obtained by
//! sdMat_extract_view and should be synchronized with accelerator data by sdMat_from/to_device
void SUBR(sdMat_extract_error)(TYPE(sdMat_ptr) M, _ST_** MC_raw, int* iflag);
#endif

//@}

//! \name data transfer between host and device (for kernel libs that support GPUs).
//!
//! These functions should return 0 if a GPU transfer is not needed, for instance because
//! the kernel lib does not support GPUs or the object was not created to live on a GPU.
//@{

//! copy multi-vector (mvec) data from the host CPU to the GPU. \ingroup mvec
void SUBR(mvec_to_device)(TYPE(mvec_ptr) V, int* iflag);

//! copy multi-vector (mvec) data from the GPU to the host CPU \ingroup mvec
void SUBR(mvec_from_device)(TYPE(mvec_ptr) V, int* iflag);

//! copy small dense matrix (sdmat) data from the host CPU to the GPU. \ingroup sdmat
void SUBR(sdMat_to_device)(TYPE(sdMat_ptr) M, int* iflag);

//! copy small dense matrix (sdmat) data from the GPU to the host CPU \ingroup sdmat
void SUBR(sdMat_from_device)(TYPE(sdMat_ptr) M, int* iflag);

//@}

//!   \defgroup mvec and sdMat I/O
/*! note: these functions only read and write the data, not the map, comm etc.
If you want to read an mvec from file, you have to create it first with a consitent map object.
 */
//@{

//! write mvec data to file
void SUBR(mvec_write_bin)(TYPE(const_mvec_ptr) V, const char* filename, int* iflag);

//! read mvec data from file
void SUBR(mvec_read_bin)(TYPE(mvec_ptr) V, const char* filename, int* iflag);

//! write sdMat data to file
void SUBR(sdMat_write_bin)(TYPE(const_sdMat_ptr) M, const char* filename, int* iflag);

//! read sdMat data from file
void SUBR(sdMat_read_bin)(TYPE(sdMat_ptr) M, const char* filename, int* iflag);

//@}

//!   \defgroup mvec Multi-vector functions (mvec_t) 
//@{

//! this function can be e.g. used to permute or redistribute vectors, the vector
//! entries of v_in will be copied into v_out, which may be based on a different map.
void SUBR(mvec_to_mvec)(TYPE(const_mvec_ptr) v_in, TYPE(mvec_ptr) v_out, int* iflag);

//! get a new vector that is a view of some columns of the original one.

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
//!     Dmvec_ptr A, Av, Avv;
//!     phist_Dmvec_create(&A,...);
//!     phist_Dmvec_view_block(A,&Av,...);
//!     phist_Dmvec_view_block(Av,&Avv,...);
//!     phist_Dmvec_delete(Av,...);
//! (do something with Avv)
void SUBR(mvec_view_block)(TYPE(mvec_ptr) V, 
                             TYPE(mvec_ptr)* Vblock,
                             int jmin, int jmax, int* iflag);

//! get a new vector that is a copy of some columns of the original one,  
//! Vblock = V(:,jmin:jmax). The object Vblock must be created beforehand 
//! and the corresponding columns of V are copied into the value array    
//! of Vblock. V is not modified.
void SUBR(mvec_get_block)(TYPE(const_mvec_ptr) V, 
                             TYPE(mvec_ptr) Vblock,
                             int jmin, int jmax, int* iflag);

//! given a multi-vector Vblock, set V(:,jmin:jmax)=Vblock by copying the corresponding
//! vectors. Vblock is not modified.
void SUBR(mvec_set_block)(TYPE(mvec_ptr) V, 
                             TYPE(const_mvec_ptr) Vblock,
                             int jmin, int jmax, int* iflag);

//@}

//!   \defgroup sdmat Small dense matrix functions (sdMat_t) 
//@{

//! get a new matrix that is a view of some rows and columns of the original one, 
//! Mblock = M(imin:imax,jmin:jmax). The behavior is analogous to mvec_view_block.
//! If on entry, *Mblock!=NULL, this function should delete *Mblock so that
//! repeated use of view_block does not lead to memory holes. It is crucial
//! that you pass in a NULL pointer if you want a new object, otherwise you
//! may get a segfault.
void SUBR(sdMat_view_block)(TYPE(sdMat_ptr) M, 
                             TYPE(sdMat_ptr)* Mblock,
                             int imin, int imax, int jmin, int jmax, int* iflag);

//! get a new matrix that is a copy of some rows and columns of the original one,  
//! Mblock = M(imin:imax,jmin:jmax). The object Mblock must be created beforehand 
//! and the corresponding columns of M are copied into the value array    
//! of Mblock. M is not modified.
void SUBR(sdMat_get_block)(TYPE(const_sdMat_ptr) M, 
                             TYPE(sdMat_ptr) Mblock,
                             int imin, int imax, int jmin, int jmax, int* iflag);

//! given a small dense matrix Mblock, set M(imin:imax,jmin:jmax)=Mblock by 
//! copying the corresponding elements. Mblock is not modified.
void SUBR(sdMat_set_block)(TYPE(sdMat_ptr) M, 
                             TYPE(const_sdMat_ptr) Mblock,
                             int imin, int imax, int jmin, int jmax, int* iflag);

//@}

//! \name initialize/fill mvecs and sdMats
//@{

//! put scalar value into all elements of a multi-vector \ingroup mvec
void SUBR(mvec_put_value)(TYPE(mvec_ptr) V, _ST_ value, int* iflag);

//! put scalar value into all elements of a small dense matrix \ingroup mvec
void SUBR(sdMat_put_value)(TYPE(sdMat_ptr) V, _ST_ value, int* iflag);

//! put random numbers into all elements of a multi-vector \ingroup mvec
void SUBR(mvec_random)(TYPE(mvec_ptr) V, int* iflag);

//! put random numbers into all elements of a small dense matrix \ingroup sdmat
void SUBR(sdMat_random)(TYPE(sdMat_ptr) V, int* iflag);

//! set all mvec elements V(i,j) by calling a function for each element
void SUBR(mvec_put_func)(TYPE(mvec_ptr) V,
        phist_mvec_elemFunc elemFunPtr, void* last_arg, int *iflag);

//! put identity matrix into a small dense matrix \ingroup sdmat

//! If M is not square and k=min(nrows,ncols), M(1:k,1:k)=I and all other entries are set to 0.
//! If *iflag|PHIST_SDMAT_RUN_ON_HOST, only the host memory is  updated. Otherwise (by default), 
//! both host and device mem  are updated (if applicable).
void SUBR(sdMat_identity)(TYPE(sdMat_ptr) M, int* iflag);

//@}

//! print a vector to the screen (for debugging) \ingroup mvec
void SUBR(mvec_print)(TYPE(const_mvec_ptr) V, int* iflag);

//! print an sdMat to the screen (for debugging) \ingroup sdmat
void SUBR(sdMat_print)(TYPE(const_sdMat_ptr) M, int* iflag);

//! \name Numerical functions
//@{

//! column-wise 2-norm \ingroup mvec

//! compute the 2-norm) of each column of v
//! (vnrm[i] must be pre-allocated by caller)
void SUBR(mvec_norm2)(TYPE(const_mvec_ptr) V, 
                        _MT_* vnrm, int *iflag);

//! normalize each column. \ingroup mvec

//! normalize (in the 2-norm) each column of v and return ||v||_2
//! for each vector i in vnrm[i] (must be pre-allocated by caller)
void SUBR(mvec_normalize)(TYPE(mvec_ptr) V, 
                            _MT_* vnrm, int* iflag);

//! scale each column i of v and by scalar. \ingroup mvec
void SUBR(mvec_scale)(TYPE(mvec_ptr) V, 
                            _ST_ scalar, int* iflag);

//! scale each column i of v and by scalar[i]. \ingroup mvec
void SUBR(mvec_vscale)(TYPE(mvec_ptr) V, 
                            const _ST_* scalar, int* iflag);

//! y=alpha*x+beta*y. \ingroup mvec

//! This function can also be used for special cases such as
//! alpha=0 => scale y
//! alpha=1, beta=0 => copy y=x
//! alpha!=0, beta=1: 'axpy' operation
void SUBR(mvec_add_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) X,
                            _ST_ beta,  TYPE(mvec_ptr)       Y,     
                            int* iflag);

//! y[i]=alpha[i]*x[i]+beta*y[i]. \ingroup mvec
void SUBR(mvec_vadd_mvec)(const _ST_ alpha[], TYPE(const_mvec_ptr) X,
                                _ST_ beta,    TYPE(mvec_ptr)       Y,     
                          int* iflag);


//! dot product of vectors v_i and w_i, i=1..numvecs. \ingroup mvec
void SUBR(mvec_dot_mvec)(TYPE(const_mvec_ptr) V, 
                            TYPE(const_mvec_ptr) W, 
                            _ST_* vw, int* iflag);

//! inner product of two multi-vectors. \ingroup mvec

//! dense tall skinny matrix-matrix product yielding a small dense matrix
//! C=alpha*V'*W+beta*C. C is replicated on all MPI processes sharing V and W.
void SUBR(mvecT_times_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) V, 
                                       TYPE(const_mvec_ptr) W, 
                                       _ST_ beta, TYPE(sdMat_ptr) C, int* iflag);

//! W=alpha*V*C + beta*W \ingroup mvec

//! n x m multi-vector times m x k dense matrix gives n x k multi-vector
void SUBR(mvec_times_sdMat)(_ST_ alpha, TYPE(const_mvec_ptr) V, 
                                       TYPE(const_sdMat_ptr) C,
                           _ST_ beta,  TYPE(mvec_ptr) W, 
                                       int* iflag);


//! V <- V*M \ingroup mvec

//! Note that M may be rectangular with ncols<nrows, in which case only the first
//! ncols columns of V are overwritten.
//! A naive default implementation for this rather uncommon kernel is available in common/kernels_no_inplace_VC.cpp
//! so that we can easily support kernel libraries that don't have it.
void SUBR(mvec_times_sdMat_inplace)(TYPE(mvec_ptr) V, TYPE(const_sdMat_ptr) M, int *iflag);

//! W <- = V*C + W*D \ingroup mvec

//! augmented kernel with two multi-vectors and two sdMats.
//! A naive default implementation for this rather uncommon kernel is available in common/kernels_no_inplace_VC.cpp
//! so that we can easily support kernel libraries that don't have it.
void SUBR(mvec_times_sdMat_add_mvec_times_sdMat)(TYPE(const_mvec_ptr) V,
                                                 TYPE(const_sdMat_ptr) C,
                                                 TYPE(mvec_ptr) W,
                                                 TYPE(const_sdMat_ptr) D,
                                                 int* iflag);


//! B=alpha*A+beta*B. \ingroup sdmat
void SUBR(sdMat_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) A,
                            _ST_ beta,  TYPE(sdMat_ptr)       B,     
                            int* iflag);

//! B=alpha*A'+beta*B. \ingroup sdmat
void SUBR(sdMatT_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) A,
                            _ST_ beta,  TYPE(sdMat_ptr)       B,     
                            int* iflag);


//! C=beta*C+alpha*A*B. \ingroup sdmat

//! n x m small dense matrix times m x k small dense matrix gives n x k small dense matrix,
//! C=alpha*V*W + beta*C
void SUBR(sdMat_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) V, 
                                         TYPE(const_sdMat_ptr) W, 
                              _ST_ beta,       TYPE(sdMat_ptr) C,
                              int* iflag);

//! C=beta*C+alpha*V'*W. \ingroup sdmat

//! m x n conj. transposed small dense matrix times m x k small dense matrix gives n x k small dense matrix,
//! C=alpha*V'*W + beta*C
void SUBR(sdMatT_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) V, 
                                          TYPE(const_sdMat_ptr) W, 
                              _ST_ beta,        TYPE(sdMat_ptr) C,
                              int* iflag);

//! C=beta*C+alpha*V*W'. \ingroup sdmat

//! n x m small dense matrix times conj. transposed k x m small dense matrix gives n x k small dense matrix,
//! C=alpha*V*W' + beta*C
void SUBR(sdMat_times_sdMatT)(_ST_ alpha, TYPE(const_sdMat_ptr) V, 
                                          TYPE(const_sdMat_ptr) W, 
                              _ST_ beta,        TYPE(sdMat_ptr) C,
                              int* iflag);

//! \addtogroup crsmat
//@{

//! Exchange elements of x between different processes, used to overlap spMVM communication with other operations.
//! Set the flag PHIST_SPMVM_ONLY_LOCAL in the call to sparseMat_times_mvec* to indicate all data is already there!
void SUBR(sparseMat_times_mvec_communicate)(TYPE(const_sparseMat_ptr) A, TYPE(const_mvec_ptr) x, int* iflag);

//! y=alpha*A*x+beta*y.

//! The scalars alpha and beta are expected to be of the
//! same type as the entries in the vectors and matrix. Mixing of types is
//! not allowed.
void SUBR(sparseMat_times_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, 
        TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag);

//! y=alpha*A^H*x+beta*y.

//! The scalars alpha and beta are expected to be of the
//! same type as the entries in the vectors and matrix. Mixing of types is
//! not allowed. In the complex case, the conjugate transpose is used.
void SUBR(sparseMatT_times_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, 
        TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag);

//! y[i]=alpha*(A*x+shift*x) + beta*y
void SUBR(sparseMat_times_mvec_add_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A,
        _ST_ shift, TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag);

//! y[i]=alpha*(A*x[i]+shifts[i]*x[i]) + beta*y[i]
void SUBR(sparseMat_times_mvec_vadd_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A,
        const _ST_ shifts[], TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag);

//! 'tall skinny' QR decomposition, V=Q*R, Q'Q=I, R upper triangular. \ingroup mvec

//! Q is computed in place of V. If V does not have full rank, iflag>0  
//! indicates the dimension of the null-space of V. The first m-iflag   
//! columns of Q are an orthogonal basis of the column space of V, the  
//! remaining columns form a basis for the null space.                  
//!                                                                     
//! As it is quite a high demand from the kernel lib to supply a rank-  
//! revealing QR, it is not strictly necessary to provide this function.
//! The orthog routine (in core/phist_orthog.h) for instance checks for 
//! a return value of -99 (not implemented) and uses a PHIST-based      
//! implementation of CholQR instead, which only requires simpler       
//! kernels.                                                            
//! For kernel libs supporting accelerators, we assume that R is sync'd 
//! after the call, that is, no sdMat_to/from_device needs to be called.
void SUBR(mvec_QR)(TYPE(mvec_ptr) V, 
                     TYPE(sdMat_ptr) R, int* iflag);

//! create matrix from a function that returns entries row-wise

//! this is the same form in which matrices are created from functions
//! in ghost and how the test problems in essex/physics are defined.
//!
//! optional flags:
//! * PHIST_SPARSEMAT_PERM_LOCAL/GLOBAL
//! * PHIST_SPARSEMAT_PERM_OPT_SINGLE/BLOCKSPMVM
//! * PHIST_SPARSEMAT_OPT_CARP
//! * PHIST_SPARSEMAT_DIST2_COLOR (if the kernel lib supports it)
//!
void SUBR(sparseMat_create_fromRowFunc)(TYPE(sparseMat_ptr) *A, phist_const_comm_ptr comm,
        phist_gidx nrows, phist_gidx ncols, phist_lidx maxnne,
        phist_sparseMat_rowFunc rowFunPtr, void* last_arg, int *iflag);

//! create a sparse matrix from a row func and use a distribution prescribed by a given context
//! (that is, assume the same shape as another matrix)
//!
//! optional flags: as fromRowFunc
void SUBR(sparseMat_create_fromRowFuncAndContext)(TYPE(sparseMat_ptr) *vA, phist_const_context_ptr ctx,
        phist_lidx maxnne,phist_sparseMat_rowFunc rowFunPtr,void* last_arg,
        int *iflag);
                

// These are not used or tested, perhaps useful in the future?
#ifdef PHIST_KERNEL_LIB_BUILTIN
void SUBR(mvec_gather_mvecs)(TYPE(mvec_ptr) V, TYPE(const_mvec_ptr) W[], int nblocks, int *iflag);
void SUBR(mvec_scatter_mvecs)(TYPE(const_mvec_ptr) V, TYPE(mvec_ptr) W[], int nblocks, int *iflag);
#endif

//@}

#ifdef IS_COMPLEX
# ifdef IS_DOUBLE
//! mixed real/complex operation: split mvec into real and imag part.
//! if either reV or imV are NULL, it is not touched.
void SUBR(mvec_split)(TYPE(const_mvec_ptr) V, phist_Dmvec* reV, phist_Dmvec* imV, int *iflag);

//! mixed real/complex operation: copy separate real and imaginary part into complex vector
//! if either reV or imV are NULL, it is not touched.
void SUBR(mvec_combine)(TYPE(mvec_ptr) V, phist_Dconst_mvec_ptr reV, phist_Dconst_mvec_ptr imV, int *iflag);
# else
//! mixed real/complex operation: split mvec into real and imag part.
//! if either reV or imV are NULL, it is not touched.
void SUBR(mvec_split)(TYPE(const_mvec_ptr) V, phist_Smvec* reV, phist_Smvec* imV, int *iflag);

//! mixed real/complex operation: copy separate real and imaginary part into complex vector
//! if either reV or imV are NULL, it is not touched.
void SUBR(mvec_combine)(TYPE(mvec_ptr) V, phist_Sconst_mvec_ptr reV, phist_Sconst_mvec_ptr imV, int *iflag);
# endif
#endif

#ifdef __cplusplus
} //extern "C"
#endif
//@}
