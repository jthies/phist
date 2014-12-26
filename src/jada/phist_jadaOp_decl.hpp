//! create projected and shifted operator for Jacobi-Davidson,
//! op*X = (I-B*V*V')(A*X+B*X*sigma)
//! by setting the appropriate pointers. No data is copied.
//! TODO: once we do preconditioning we should probably have
//!       the preconditioner handle the projection

void SUBR(jadaOp_create)(TYPE(const_op_ptr)    A_op,    TYPE(const_op_ptr)    B_op,
                         TYPE(const_mvec_ptr)  V,       TYPE(const_mvec_ptr)  BV,
                         const _ST_            sigma[], int                   nvec,
                         TYPE(op_ptr)          jdOp,    int*                  iflag);

void SUBR(jadaOp_delete)(TYPE(op_ptr)  jdOp, int *iflag);
