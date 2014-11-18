//! The Anasazi package from Trilinos provides block Krylov solvers.
//! We currently provide access to the following methods
//!
//! method=0: block Krylov Schur
//! 
//! TODO: documentation
//!
void SUBR(anasazi)(TYPE(const_op_ptr) Op, 
        TYPE(mvec_ptr) X,
        TYPE(const_mvec_ptr) B, 
        _MT_ tol,int* num_iters, int max_blocks,
        int variant, int *nConv,
        int* ierr);
