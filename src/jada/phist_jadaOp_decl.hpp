//! create projected and shifted operator for Jacobi-Davidson,
//! op*X = (I-B*V*V')(A*X+B*X*sigma)
//! by setting the appropriate pointers. No data is copied.
//! TODO: once we do preconditioning we should probably have
//!       the preconditioner handle the projection

void SUBR(jadaOp_create)(TYPE(const_linearOp_ptr)    A_op,    TYPE(const_linearOp_ptr)    B_op,
                         TYPE(const_mvec_ptr)  V,       TYPE(const_mvec_ptr)  BV,
                         const _ST_            sigma[], int                   nvec,
                         TYPE(linearOp_ptr)          jdOp,    int*                  iflag);

void SUBR(jadaOp_delete)(TYPE(linearOp_ptr)  jdOp, int *iflag);

//! create a preconditioner for the inner solve in Jacobi-Davidson.
//!
//! Given a linear operator that is a preconditioner for A, this function will simply
//! wrap it up to use apply_shifted when apply() is called. We need this because our implementations
//! of blockedGMRES and MINRES are not aware of the shifts so they can only call apply in the precon-
//! ditioning operator. Obviously not all preconditioners are able to handle varying shifts without
//! recomputing, this is not taken into account by this function:in that case the input P_op must be
//! updated beforehand.
void SUBR(jadaPrec_create)(TYPE(const_linearOp_ptr) P_op, 
                           const _ST_ sigma[], int nvec, 
                           TYPE(linearOp_ptr) jdPrec,
                           int* iflag);
