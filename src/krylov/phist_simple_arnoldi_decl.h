//! a simple Arnoldi process to start up the JaDa iteration.
//! Given a minimum basis size m, compute V(:,1:m+1), 
//! H(1:m+1,1:m) such that A*V(:,1:m) = V(:,1:m+1)*H(1:m+1,1:m)
//! input: v0, V and H allocated with m+1 resp. m columns
//! and nloc resp. m+1 rows.
//!
//! If B_op!=NULL, B is assumed to be symmetric and positive definite 
//! and the B-inner product is used and the Arnoldi-relation is
//! A*V(:,1:m) = B*V(:,1:m+1)*H(1:m+1,1:m).
//!
//! TODO - we may want to include a check if any Ritz values have already
//!        converged in Arnoldi, but this requires some additional programming
//!        effort which we leave for later.
void SUBR(simple_arnoldi)(TYPE(const_linearOp_ptr) A_op, TYPE(const_linearOp_ptr) B_op, TYPE(const_mvec_ptr) v0,
        TYPE(mvec_ptr) V, TYPE(mvec_ptr) AV, TYPE(mvec_ptr) BV, TYPE(sdMat_ptr) H, int m, int* iflag);

//! starts with random vector!
void SUBR(simple_blockArnoldi)(TYPE(const_linearOp_ptr) A_op, TYPE(const_linearOp_ptr) B_op,
                               TYPE(mvec_ptr) V, TYPE(mvec_ptr) AV, TYPE(mvec_ptr) BV,
                               TYPE(sdMat_ptr) H, int m, int blockSize, int* iflag);
