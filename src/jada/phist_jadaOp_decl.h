// create projected and shifted operator for Jacobi-Davidson,
// op*X = (B-VV')(A-shift*B)
// by setting the appropriate pointers. No data is copied.
// TODO: allow complex shifts and/or 'blocks as shifts' (??)
// TODO: once we do preconditioning we should probably have
//       the preconditioner handle the projection
void SUBR(jadaOp_create)(TYPE(const_op_ptr) A_op, TYPE(const_op_ptr) B_op,
                         TYPE(const_mvec_ptr) V, _ST_ shift, 
                         TYPE(op_ptr)* jdOp, int *ierr);

