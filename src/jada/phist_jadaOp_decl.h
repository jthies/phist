//! create projected and shifted operator for Jacobi-Davidson,
//! op*X = (I-BV*V')(A*X+B*X*sigma)(I-V*BV')
//! by setting the appropriate pointers. No data is copied.
//! For B_op==NULL, the pre-projection I-V*BV') is
//! omitted to save reductions (and BV is ignored).
//! 
//! We make use of the apply_shifted function in the given A_op object. To implement
//! the operation (A+sigma*B)X for B!=I, the A_op should implement this in apply_shifted.
//! To construct such an operator from two sparseMats, SUBR(linearOp_wrap_sparseMat_pair)
//! can be used.
void SUBR(jadaOp_create)(TYPE(const_linearOp_ptr)    AB_op,
                         TYPE(const_linearOp_ptr)     B_op,
                         TYPE(const_mvec_ptr)  V,       TYPE(const_mvec_ptr)  BV,
                         const _ST_            sigma[], int                   nvec,
                         TYPE(linearOp_ptr)          jdOp,    int*                  iflag);

void SUBR(jadaOp_delete)(TYPE(linearOp_ptr)  jdOp, int *iflag);

//! create a preconditioner for the inner solve in Jacobi-Davidson.
//!
//! Given a linear operator P_op implementing y <- P\x, where P is a preconditioner for A, this function will
//! wrap it up to use apply_shifted when apply() is called. We need this because our implementations
//! of blockedGMRES and MINRES are not aware of the shifts so they can only call apply in the precon-
//! ditioning operator. Obviously not all preconditioners are able to handle varying shifts without
//! recomputing, this is not taken into account by this function: in that case the input P_op must be
//! updated beforehand.
//!
//! If V (and possibly BV) are given, the resulting operator will be
//!
//!  Y <- (I - (P_op\V)*((BV)'P_op\V)^{-1} (BV)') P_op\y
//!
//! Note 1: this operator is always non-symmetric, if a symmetric variant is needed we would have to
//! pre- and postproject.
//!
//! Note 2: Our JaDa implementation does not store P_op\Q in order to save memory. When alling this
//! function, P\V is computed and stored until the operator is deleted. If V=Q is given, this may   
//! be quite a computational overhead because before solving the correction equation, P_op is applied
//! to all locked eigenvectors in Q. We therefore advocate using only the current approximation, V=q,
//! for the projection of the preconditioner.
//!
void SUBR(jadaPrec_create)(TYPE(const_linearOp_ptr) P_op, 
                           TYPE(const_mvec_ptr)  V,       TYPE(const_mvec_ptr)  BV,
                           const _ST_ sigma[], int nvec, 
                           TYPE(linearOp_ptr) jdPrec,
                           int* iflag);
