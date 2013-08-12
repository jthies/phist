//! tries to compute a given number of eigenpairs (num_eig) of 
//! an hermitian positive definite operator A using the Lanczos
//! method. 
//!
//! input arguments:
//!
//! A: pointer to the operator(-wrapper) A
//! X:  start vectors for the num_eigs desired
//!     eigenvectors (random vectors if you don't know better)
//! tol: convergence tolerance
//! num_eigs: number of desired eigenpairs
//! num_eigs: maximum number of iterations allowed
//! evals and resid should be pre-allocated arrays of size at least num_eigs
//! 
//! output arguments:
//!
//! num_eigs: number of converged eigenpairs
//! num_iters: number of iterations performed
//! X: latest eigenvector approximations
//! evals: latest eigenvalue approximations
//! resid: Ritz residuals
//! ierr: return code of the solver (0 on success, negative on error, positive on warning)
//!
//! Note: This solver is not fully implemented, it just computes a few eigenvalues without
//! keeping the basis, so eigenvectors cannot be obtained, no reorthogonalization or restart
//! is performed and the quality of the resulting eigenvalues may be poor.
//!
void SUBR(lanczos)(TYPE(const_op_ptr) Op, 
        TYPE(mvec_ptr) X, 
        _MT_* evals, _MT_* resid,
        _MT_ tol,int* num_iters, int* num_eigs,
        int* ierr);
         
         
