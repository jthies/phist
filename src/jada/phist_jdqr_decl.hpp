//! tries to compute a given number of eigenpairs (num_eig) of 
//! an general operator A using the Jacobi-Davidson QR (JDQR)
//! method. We allow for generalized EVPs here, A*x=lambda*B*x,
//! where B should be positive definite.
//!
//! input arguments:
//!
//! A: pointer to the operator(-wrapper) A
//! B: pointer to the "mass-matrix" operator (if B==NULL, B=I is assumed)
//! X:  start vectors for the num_eigs desired                           
//!     eigenvectors (random vectors if you don't know better)           
//!     (currently ignored on input). Should be allocated with           
//!     at least <num_eigs> columns (num_eigs+1 in the real case         
//!     to avoid the situation where the last eigenvalue to con-         
//!     verge is a complex pair).                                        
//!                                                                      
//! which - decide at which end of the spectrum to look for eigenvalues
//! tol: convergence tolerance
//! num_eigs: number of desired eigenpairs
//! num_iters: maximum number of iterations allowed
//! is_cmplx,evals and resid should be pre-allocated arrays of size at least num_eigs
//! (in the real case num_eigs+1, cf. comment for 'X' above).
//! 
//! minBas: start up from a basis consisting of minBas vectors (using Arnoldi)
//! maxBas: when the basis reaches <maxBas> vectors, restart from <minBas> vectors.
//! 
//! output arguments:
//!
//! num_eigs: number of converged eigenpairs
//! num_iters: number of iterations performed
//! X: latest eigenvector approximations
//! evals: latest eigenvalue approximations
//! resid: Ritz residuals
//! is_cmplx: in the real case, if is_cmplx[i]==is_cmplx[i+1]=1, then
//! a evals[i] +/- evals[i+1]*I forms a complex conjugate pair,
//! and the corresponding eigenvectors are X(:,i)+/-X(:,i+1)*I.
//! In the complex case is_cmplx is not referenced.
//! ierr: return code of the solver (0 on success, negative on error, positive on warning)
//!
void SUBR(jdqr)(TYPE(const_op_ptr) A_op, TYPE(const_op_ptr) B_op,
        TYPE(mvec_ptr) X, eigSort_t which,
        _ST_* evals, _MT_* resid, int* is_cmplx,
        _MT_ tol, int* num_eigs, int* num_iters,
        int minBas, int maxBas,
        int* ierr);
         
         
