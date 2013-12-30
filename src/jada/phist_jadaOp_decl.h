//! create projected and shifted operator for Jacobi-Davidson,
//! op*X = (I-B*V*V')(A*X+B*X*sigma)
//! by setting the appropriate pointers. No data is copied.
//! TODO: once we do preconditioning we should probably have
//!       the preconditioner handle the projection

void SUBR(jadaOp_create)(TYPE(const_op_ptr)    A_op,    TYPE(const_op_ptr)    B_op,
                         TYPE(const_mvec_ptr)  V,       TYPE(const_mvec_ptr)  BV,
                         const _ST_            sigma[], TYPE(sdMat_ptr)       VY,
                         TYPE(mvec_ptr)        AX,      TYPE(mvec_ptr)        BX,
                         TYPE(mvec_ptr)        X_proj,  TYPE(mvec_ptr)        work,
                         TYPE(op_ptr)          jdOp,    int*                  ierr);

void SUBR(jadaOp_delete)(TYPE(op_ptr)  jdOp, int *ierr);

//! access AX from last call to apply (return a view to it)
void SUBR(jadaOp_view_AX)(TYPE(const_op_ptr) jadaOp, TYPE(mvec_ptr)*AX, int* ierr);

//! access BX from last call to apply (return a view to it)
void SUBR(jadaOp_view_BX)(TYPE(const_op_ptr) jadaOp, TYPE(mvec_ptr)*BX, int* ierr);

//! access X_proj from last call to apply (return a view to it)
void SUBR(jadaOp_view_X_proj)(TYPE(const_op_ptr) jadaOp, TYPE(mvec_ptr)*X_proj, int* ierr);

