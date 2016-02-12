//! \addtogroup jada
//@{

//! \defgroup bjdqr Jacobi-Davidson QR method with blocking, locking and restart
//@{


//! Subspace Jacobi-Davidson for extreme eigenvalues

//! Tries to compute a partial schur form $(Q,R)$ of dimension opts.nEigs
//! of the stencil $A*x-\lambda*B*x$ with a general linear operator $A$ and a
//! hermitian positive definite (hpd.) linear operator $B$ using a
//! block-Jacobi-Davidson QR method. <br>
//! The generalized eigenvalues $\lambda_i$ are the diagonal entries of the
//! partial schur form $A*Q = B*Q*R$ returned. <br>
//!
//! solver options are passed in via a jadaOpts_t struct, default parameters can
//! be set using jadaOpts_setDefaults().
//!
//! Input arguments:
//!
//! A_op:             pointer to the operator A
//! B_op:             pointer to the hpd. operator B (if B==NULL, B=I is assumed)
//!
//! for documentation of the other input parameters, see file phist_jadaOpts.h
//! 
//! Output arguments:
//!
//! nConv:     number of converged eigenpairs (e.g. dimension of (Q,R))
//! Q:        orthogonal vectors of the partial schur form (Q,R)
//! R:        small upper triangular matrix of the partial schur form (Q,R)
//! nIter:    number of iterations performed
//! resNorm:  norm of the residua of the schur form $A*q_i-Q*r_i, i=1,nEig$
//! iflag:     return code of the solver (0 on success, negative on error, positive on warning)
//!
void SUBR(subspacejada)( TYPE(const_linearOp_ptr) A_op,  TYPE(const_linearOp_ptr) B_op,
                         phist_jadaOpts_t opts,
                         TYPE(mvec_ptr) Q,         TYPE(sdMat_ptr) R,
                         _CT_* ev,                 _MT_* resNorm,
                         int *nConv,                int *nIter,
                        int* iflag);

//! Subspace Jacobi-Davidson for interior eigenvalues

//! This algorithm is very similar to subspacejada, the main differences are    
//! that eigenvalues near a specified target inside the spectrum are sought.    
//! Instead of starting with some Arnoldi steps, we keep the shift fixed to the 
//! target until the minimum basis size is reached. The scheme computes V and W 
//! such that V'V=W'W=I, and W is an orthogonal basis of A*V (A*V=W*H_A). Eigen-
//! values are approximated using harmonic Ritz values
void SUBR(harmonicjada)( TYPE(const_linearOp_ptr) A_op,  TYPE(const_linearOp_ptr) B_op,
                         phist_jadaOpts_t opts,
                         TYPE(mvec_ptr) Q,         TYPE(sdMat_ptr) R,
                         _CT_* ev,                 _MT_* resNorm,
                         int *nConv,                int *nIter,
                        int* iflag);



//@}

//@}
