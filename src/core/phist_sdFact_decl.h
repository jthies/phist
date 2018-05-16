/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/

//! \addtogroup sdFact
//!@{

//! \brief cholesky decomposition
//!
//! stable cholesky factorization with pivoting and rank-recognition for hp(s)d. matrix 
//! returns column-permuted upper triangular cholesky factor R in M s.t. M=R'*R.        
//! If M is singular with rank(M)=r<m (M being m x m), the last m-r rows of R will be 0.
//! On exit, *rank=r is returned.
void SUBR(sdMat_cholesky)(TYPE(sdMat_ptr) M, int* perm, int* rank, _MT_ rankTol, int* iflag);

//! \brief backward substitution
//!
//! backward substitution for pivoted upper triangular cholesky factor
void SUBR(sdMat_backwardSubst_sdMat)(const TYPE(sdMat_ptr) R, int* perm, int rank, TYPE(sdMat_ptr) X, int* iflag);

//! \brief forward substitution
//!
//! forward substitution for pivoted conj. transposed upper triangular cholesky factor
void SUBR(sdMat_forwardSubst_sdMat)(const TYPE(sdMat_ptr) R, int* perm, int rank, TYPE(sdMat_ptr) X, int* iflag);

//! \brief given B=V'V, compute (in place) B^ s.t. V*B^ is orthonormal.
//!
//! The inverse of B^ is returned in B_1 (if it is not NULL), the 
//! rank of V in *rank.
//! If V does not have full rank, the last n-*rank columns of B and B_1
//! will be zero.
void SUBR(sdMat_qb)(TYPE(sdMat_ptr) B, TYPE(sdMat_ptr) B_1,int* rank, _MT_ rankTol, int* iflag);

//! \brief computes in-place the inverse of a Hermitian and positive semi-definite matrix using Cholesky factorization.
//!
//! If A is singular (actually semi-definite, that is), the pseudo-inverse is computed using rank-revealing Cholesky.
//! The rank of A on input is returned as *rank.
void SUBR(sdMat_inverse)(TYPE(sdMat_ptr) A_hpd, int* rank, _MT_ rankTol, int* iflag);

//! \brief computes in-place the (transposed) Moore-Penrose pseudo-inverse A+ of an arbitrary n x m matrix A.
//!
//! The rank of A on input is returned as *rank.
//! The four defining properties are (where A is the input and B the output matrix):
//!
//! 1. AB'A = A <br>
//! 2. B'AB'= B' <br>
//! 3. (AB')'=AB' <br>
//! 4. (B'A)'=B'A
void SUBR(sdMat_pseudo_inverse)(TYPE(sdMat_ptr) A_gen, int* rank, _MT_ rankTol, int* iflag);

//! \brief singular value decomposition, A = U*Sigma*Vt
//!
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

//! \brief given a m x m matrix A, compute the explicit QR factorization A=Q*R.
//!
//! A should be provided via the R argument.
//! If A is singular, some columns of Q may be 0. The function calls the lapack
//! subroutines XGEQR and XGEMQR to explicitly construct the Q factor. On output,
//! A is overwritten by R.
void SUBR(sdMat_qr)(TYPE(sdMat_ptr) Q, TYPE(sdMat_ptr) R, int* iflag);

//!@}
