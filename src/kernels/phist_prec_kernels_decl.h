//! \defgroup prec mixed-precision kernels implemented directly in PHIST.
//! these kernels are highly accurate but not necessarily efficient, they
//! are intended for use with sdMats, that is the matrix should be small 
//! enough s.t. serial execution won't hurt even on a supercomputer.
//! 
//! The way to use these kernels is by
//! compiling with a kernel lib that allows you to enable PHIST_HIGH_RECISION_KERNELS
//! (e.g. the builtin kernels on Haswell)
//! Note: you should always make sure that on hybrid CPU/GPU systems the sdMat values on the
//! host and device are synchronized by calling sdMat_from(to)_device before (after) using 
//! these functions. They all assume that the data obtained from sdMat_extract_view/error is
//! up-to-date with the device, and do not call to_device afterwards.

#ifdef __cplusplus
extern "C" {
#endif

//! cholesky decomposition. \ingroup prec

//! stable cholesky factorization with pivoting and rank-recognition for hpd. matrix
//! returns permuted lower triangular cholesky factor M for M <- M*M'
void SUBR(sdMat_cholesky)(TYPE(sdMat_ptr) M, int* perm, int* rank, int* iflag);

//! backward substitution. \ingroup prec

//! backward substitution for pivoted upper triangular cholesky factor
void SUBR(sdMat_backwardSubst_sdMat)(const TYPE(sdMat_ptr) R, int* perm, int rank, TYPE(sdMat_ptr) X, int* iflag);

//! forward substitution. \ingroup prec

//! forward substitution for pivoted conj. transposed upper triangular cholesky factor
void SUBR(sdMat_forwardSubst_sdMat)(const TYPE(sdMat_ptr) R, int* perm, int rank, TYPE(sdMat_ptr) X, int* iflag);

//! given B=V'V, compute (in place) B^ s.t. V*B^ is orthonormal. 
//! The inverse of B^ is returned in B_1 (if it is not NULL), the 
//! rank of V in *rank.
//! If V does not have full rank, the last n-*rank columns of B and B_1
//! will be zero.
void SUBR(sdMat_qb)(TYPE(sdMat_ptr) B, TYPE(sdMat_ptr) B_1,int* rank, int* iflag);

#ifdef __cplusplus
} // extern "C"
#endif
