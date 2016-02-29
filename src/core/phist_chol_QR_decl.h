//! rank revealing Cholesky-QR with null space randomization

//! SVQB algorithm (Stathopoulos & Wu, SISC 23 (6),2165-2182)
//! with rank-revealing pivoted cholesky (SVRR)
//! If the input vector has m columns and rank r, iflag=m-r is
//! returned. The nullspace of V is randomized and orthogonalized
//! against the other columns. If the output matrix is denoted
//! by Q, Q and V are related by Q*R=V.
//!
//! This routine is the fallback kernel used by orthog if the kernel
//! library doesn't provide mvec_QR.
void SUBR(chol_QR)(TYPE(mvec_ptr) V, TYPE(sdMat_ptr) R, int* iflag);

