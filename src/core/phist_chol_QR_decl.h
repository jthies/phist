//! rank revealing Cholesky-QR with pivoting and null space randomization

//! SVQB algorithm (Stathopoulos & Wu, SISC 23 (6),2165-2182)
//! with rank-revealing pivoted cholesky (SVRR)
//! If the input vector has m columns and rank r, iflag=m-r is
//! returned. The nullspace of V is randomized and orthogonalized
//! against the other columns. If the output matrix is denoted
//! by Q, Q and V are related by Q*R=V. R is 'psychologically'
//! upper triangular, that is: R(:,perm) is upper triangular.
//!
//! This routine is the fallback kernel used by orthog if the kernel
//! library doesn't provide mvec_QR.
//!
//! Before returning, R is synchronized with devices using
//! sdMat_to_device so that it is up-to-date on both host and
//! device (if applicable).
void SUBR(chol_QRp)(TYPE(mvec_ptr) V, TYPE(sdMat_ptr) R, int perm[], int* iflag);

//! this variant of chol_QR discards the perm array. It has the same interface
//! as mvec_QR but returns a column-permuted upper triangular matrix R.
void SUBR(chol_QR)(TYPE(mvec_ptr) V, TYPE(sdMat_ptr) R, int* iflag);
