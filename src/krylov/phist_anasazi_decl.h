//! The Anasazi package from Trilinos provides block eigensolvers.
//! We currently provide access to the following methods
//!
//! block Krylov Schur  (no choice given to the user)
//! 
//! input parameters: similar to subspacejada (TODO: documentation)
//!
void SUBR(anasazi)(      TYPE(const_op_ptr) A_op, TYPE(const_op_ptr) B_op,
                         TYPE(const_mvec_ptr) v0,  eigSort_t which,
                         _MT_ tol,                 int *nEig,
                         int* nIter,               int blockDim,
                         int numBlocks,
                         bool symmetric,
                         TYPE(mvec_ptr) vX,         _ST_* eigs,
                         _MT_* resNorm,            int* ierr);

