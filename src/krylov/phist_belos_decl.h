//! block GMRES - solve a linear system with multiple RHS
//! using the BlockGMRES method implemented in Belos.
//! AX=B is solved for X. At most max_blocks blocks of vectors are
//! generated. A block consists of as many vectors as there
//! are columns in X and B. If max_blocks is reached, the  
//! method is restarted, until num_iters is reached or     
//! ||r||/||r0||<tol is achieved. *num_iters is overwritten
//! by the actual number of iterations performed.
//!
//! We currently interface two versions of GMRES, selected by
//! the parameter 'variant':
//! 0: block gmres, single subspace expanded by k vectors at a
//!    time (k being the number of rhs)
//! 1: 'pseudo' block gmres, k indep. subspaces built. Equivalent
//!    to running a single GMRES for each rhs.
//! for variant 0, nConv is not touched. For variant 1 it indicates
//! (on entry) the minimum number of systems that should have converged
//! before returning from this function. If nConv==NULL, nConv=k is assumed.
void SUBR(bgmres)(TYPE(const_op_ptr) Op, 
        TYPE(mvec_ptr) X,
        TYPE(const_mvec_ptr) B, 
        _MT_ tol,int* num_iters, int max_blocks,
        int variant, int *nConv,
        int* ierr);
