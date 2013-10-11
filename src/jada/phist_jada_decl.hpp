//! tries to compute a given number of eigenpairs (num_eig) of 
//! an hermitian positive definite operator A using the Lanczos
//! method. 
//!
//! input arguments:
//!
//! A: pointer to the operator(-wrapper) A
//! X:  start vectors for the num_eigs desired
//!     eigenvectors (random vectors if you don't know better)
//!     (currently ignored on input)
//! which - decide at which end of the spectrum to look for eigenvalues
//! tol: convergence tolerance
//! num_eigs: number of desired eigenpairs
//! num_iters: maximum number of iterations allowed
//! is_cmplx,evals and resid should be pre-allocated arrays of size at least num_eigs
//! (in the real case num_eigs*2 so it can hold real- and imaginary part of 
//! complex eigenvalues).
//! 
//! output arguments:
//!
//! num_eigs: number of converged eigenpairs
//! num_iters: number of iterations performed
//! X: latest eigenvector approximations
//! evals: latest eigenvalue approximations
//! resid: Ritz residuals
//! is_cmplx: if is_cmplx[i]==is_cmplx[i+1]=1, then
//! a evals[i] +/- evals[i+1]*I forms a complex conjugate pair,
//! and the corresponding eigenvectors are X(:,i)+/-X(:,i+1)*I.
//! In the complex case is_cmplx is not referenced.
//! ierr: return code of the solver (0 on success, negative on error, positive on warning)
//!
void _SUBR_(jada)(_TYPE_(const_op_ptr) Op, 
        _TYPE_(mvec_ptr) X, eigSort_t which,
        _ST_* evals, _MT_* resid, int* is_cmplx,
        _MT_ tol,int* num_iters, int* num_eigs,
        int minBas, int maxBas,
        int* ierr);
         
         
