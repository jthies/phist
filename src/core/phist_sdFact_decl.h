//! \defgroup sdfact sdMat kernels for factoring small dense matrices.
//!
//! These kernels are highly accurate but not necessarily efficient, If the
//! kernel lib supports PHIST_HIGH_PRECISION_KERNELS (e.g. the builtin kernels
//! with AVX2), the input *iflag=PHIST_ROBUST_REDUCTIONS enables the high precision
//! variants. Otherwise, standard precision is used.
//! 
//! Note: you should always make sure that on hybrid CPU/GPU systems the sdMat values on the
//! host and device are synchronized by calling sdMat_from(to)_device before (after) using 
//! these functions. They all assume that the data obtained from sdMat_extract_view/error is
//! up-to-date with the device, and do not call to_device afterwards.

//! cholesky decomposition. \ingroup sdFact

//! stable cholesky factorization with pivoting and rank-recognition for hp(s)d. matrix 
//! returns column-permuted upper triangular cholesky factor R in M s.t. M=R'*R.        
//! If M is singular with rank(M)=r<m (M being m x m), the last m-r rows of R will be 0.
//! On exit, *rank=r is returned.
void SUBR(sdMat_cholesky)(TYPE(sdMat_ptr) M, int* perm, int* rank, int* iflag);

//! backward substitution. \ingroup sdFact

//! backward substitution for pivoted upper triangular cholesky factor
void SUBR(sdMat_backwardSubst_sdMat)(const TYPE(sdMat_ptr) R, int* perm, int rank, TYPE(sdMat_ptr) X, int* iflag);

//! forward substitution. \ingroup sdFact

//! forward substitution for pivoted conj. transposed upper triangular cholesky factor
void SUBR(sdMat_forwardSubst_sdMat)(const TYPE(sdMat_ptr) R, int* perm, int rank, TYPE(sdMat_ptr) X, int* iflag);

//! given B=V'V, compute (in place) B^ s.t. V*B^ is orthonormal. 
//! The inverse of B^ is returned in B_1 (if it is not NULL), the 
//! rank of V in *rank.
//! If V does not have full rank, the last n-*rank columns of B and B_1
//! will be zero.
void SUBR(sdMat_qb)(TYPE(sdMat_ptr) B, TYPE(sdMat_ptr) B_1,int* rank, int* iflag);

//! computes in-place the inverse of a Hermitian and positive semi-definite matrix using Cholesky factorization.
//! If A is singular (actually semi-definite, that is), the pseudo-inverse is computed using rank-revealing Cholesky.
//! The rank of A on input is returned as *rank.
void SUBR(sdMat_inv)(TYPE(sdMat_ptr) A_hpd, int* rank, int* iflag);

//! computes in-place the pseudo-inverse of an arbitrary square matrix. The rank of A on input is returned as *rank.
void SUBR(sdMat_pinv)(TYPE(sdMat_ptr) A_gen, int* rank, int* iflag);

//! singular value decomposition, A = U*Sigma*Vt

//! A, Sigma are  m x n, U m x m, V n x n.
//! The function just calls the lapack routine XGESVD with JOBU=JOBV='A',
//! so all the left and right singular vectors are computed and A is filled
//! with garbage on output. The transpose of V is returned, Vt = transpose(V),
//! and the singular values are sorted on the diagonal of Sigma by decreading
//! magnitude.
//!
//! We also allow Sigma to have dimension (min(m,n),1), in that case we return the
//! diagonal entries of Sigma only (the actual singular values).
void SUBR(sdMat_svd)(TYPE(sdMat_ptr) A, TYPE(sdMat_ptr) U, TYPE(sdMat_ptr) Sigma, TYPE(sdMat_ptr) Vt, int* iflag);
