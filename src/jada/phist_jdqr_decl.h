//! \defgroup jada Jacobi-Davidson style eigenvalue solvers
//@{

//! \defgroup jdqr single-vector JDQR method with deflation and restart
//@{

//! tries to compute a given number of eigenpairs (num_eig) of 
//! an general operator A using the Jacobi-Davidson QR (JDQR)
//! method. We allow for generalized EVPs here, A*x=lambda*B*x,
//! where B should be positive definite.
//!
//! input arguments:
//!
//! A: pointer to the operator(-wrapper) A
//! B: pointer to the "mass-matrix" operator (if B==NULL, B=I is assumed)
//!
//! X:  Should be allocated with           
//!     at least <num_eigs> columns (num_eigs+1 in the real case         
//!     to avoid the situation where the last eigenvalue to con-         
//!     verge is a complex pair).                                        
//!                                                                      

//! is_cmplx,evals and resid should be pre-allocated arrays of size at least opts.numEigs
//! (in the real case opts.numEigs+1)
//! 
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
//! iflag: return code of the solver (0 on success, negative on error, positive on warning)
//!
void SUBR(jdqr)(TYPE(const_linearOp_ptr) A_op, TYPE(const_linearOp_ptr) B_op,
  TYPE(mvec_ptr) X, TYPE(mvec_ptr) Q, TYPE(sdMat_ptr) R,
  _ST_* evals, _MT_* resid, int* is_cmplx,
        phist_jadaOpts_t options, int* num_eigs, int* num_iters,
        int* iflag);
         
//@}

//@}
