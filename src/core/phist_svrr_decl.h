//! rank-revealing "Cholesky QR" of an mvec

//! SVQB algorithm (Stathopoulos & Wu, SISC 23 (6),2165-2182)
//! with rank-revealing pivoted cholesky (SVRR)
//! If the input vector has m columns and rank r, iflag=m-r is
//! returned. Columns 0:iflag-1 of V will be zero, the remaining
//! columns will be orthogonal. If the output matrix is denoted
//! by Q, Q and V are related by Q=V*B.
//!
//! Accepts the input argument *iflag=PHIST_ROBUST_REDUCTIONS, in
//! which case double-double arithmetic is used if the kernel lib
//! supports it.
//!
//! \deprecated in contrast to our common practice, this function does
//! not randomize the nullspace (if any). The function chol_QR
//! implements the semantics of the optional kernel function mvec_QR and
//! should therefore be preferred.
void SUBR(svrr)(TYPE(mvec_ptr) V, TYPE(sdMat_ptr) B, int* iflag);
