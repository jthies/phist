// get type-sepcific macro wrappers for lapack,
// for instance XGEMM <- SGEMM, DGEMM, CGEMM or ZGEMM 
#include "phist_lapack_decl.h"

static int _SUBR_(converged)(_ST_ evmin, _MT_ tol)
{
	static _ST_ oldevmin = -1e9;

	int converged = _ABS_(evmin-oldevmin) < tol;
	oldevmin = evmin;
	return converged;
}
// on entry, vnew = r, vold=v_(j-1).
// on exit, vnew=v_j, vold=r
static void _SUBR_(lanczosStep)(_TYPE_(const_op_ptr) op, 
                _TYPE_(mvec_ptr) vnew, 
                _TYPE_(mvec_ptr) vold,
		_ST_ *alpha, _ST_ *beta, int* ierr)
{
        *ierr=0;
        // vnew = r/beta
        _ST_ ibeta = _ONE_/(*beta);
        _PHIST_ERROR_HANDLER_(_SUBR_(mvec_scale)(vnew,&ibeta,ierr),*ierr);
        // r = A*v - beta*vold 
        _PHIST_ERROR_HANDLER_(op->apply(_ONE_,op->A_,vnew, 
                -(*beta),vold,ierr),*ierr);
        //alpha_j = v_j'r
        _PHIST_ERROR_HANDLER_(_SUBR_(mvec_dot_mvec)(vnew,vold,alpha,ierr),*ierr);
        // r=r-alpha_j*v_j
	_PHIST_ERROR_HANDLER_(_SUBR_(mvec_add_mvec)(-(*alpha),vnew,_ONE_,vold,ierr),*ierr);
	// beta_j = ||r||_2
	_PHIST_ERROR_HANDLER_(_SUBR_(mvec_dot_mvec)(vold,vold,beta,ierr),*ierr);
	*beta=_SQRT_(*beta);
        return;
}

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
void _SUBR_(lanczos)(_TYPE_(const_op_ptr) op, 
        _TYPE_(mvec_ptr) X, 
        _ST_* evals, _MT_* resid,
        _MT_ tol,int* num_iters, int* num_eigs,
        int* ierr)
  {
  _TYPE_(mvec_ptr) vold, vnew, vtmp;
  
  int nIter=*num_iters;
  int it,n;
  _ST_ alpha, beta;
  _MT_ nrm;
 
   // allocate memory for the tridiagonal matrix H = [-beta_(j-1) alpha_j -beta_j]
  _ST_ *alphas  = (_ST_ *)malloc(sizeof(_ST_)*nIter);
  _ST_ *betas   = (_ST_ *)malloc(sizeof(_ST_)*nIter);
  _ST_ *falphas = (_ST_ *)malloc(sizeof(_ST_)*nIter);
  _ST_ *fbetas  = (_ST_ *)malloc(sizeof(_ST_)*nIter);
 
  int ldz=1;
  _MT_ *work=NULL;
  _ST_ *x_ptr; // TODO - we don't compute eigenvectors yet
  //int i,lda,nloc;
  *num_iters=0;
  *num_eigs=0;
  
  if (op->range_map_ != op->domain_map_)
    {
    *ierr=-1; // maps of operator should point to same object (A should be a square matrix)
    }

  _PHIST_ERROR_HANDLER_(_SUBR_(mvec_create)(&vold,op->domain_map_,1,ierr),*ierr);
  _PHIST_ERROR_HANDLER_(_SUBR_(mvec_create)(&vnew,op->domain_map_,1,ierr),*ierr);
  
  // vold = random start vector, r=vold
  _PHIST_ERROR_HANDLER_(_SUBR_(mvec_random)(vold,ierr),*ierr);
  _PHIST_ERROR_HANDLER_(_SUBR_(mvec_normalize)(vold,&nrm,ierr),*ierr);
  _PHIST_ERROR_HANDLER_(_SUBR_(mvec_add_mvec)(_ONE_,vold,_ZERO_,vnew,ierr),*ierr);

  beta = _ONE_;
  betas[0] = beta;

  for(it = 0, n=1; 
      it < nIter; 
      it++, n++)
    {
         _PHIST_ERROR_HANDLER_(_SUBR_(lanczosStep)(op,vnew,vold,&alpha,&beta,ierr),*ierr);
         vtmp=vnew;
         vnew=vold;
         vold=vtmp;
         alphas[it]=alpha;
         betas[it+1]=beta;
/*
         fprintf(stdout,"iter %d, alpha=%8.4f %8.4fi, \t beta=%8.4f+%8.4fi, \n", 
                it,
                _REAL_(alpha),_IMAG_(alpha),
                _REAL_(beta),_IMAG_(beta));
*/
        if (it%5==0)
          {
          memcpy(falphas,alphas,n*sizeof(_ST_)); // alphas and betas will be destroyed 
          memcpy(fbetas,betas,n*sizeof(_ST_));

          _PHIST_ERROR_HANDLER_(XSTEQR("N",&n,falphas,fbetas,x_ptr,&ldz,work,ierr),*ierr);
         
           printf("it %d lambda_min: %f\n",it,_REAL_(falphas[0]));
           if (_SUBR_(converged(falphas[0],tol))) break; 
           }
         }

  *num_iters=it;
  *num_eigs= it<nIter? 1:0;

  evals[0] = falphas[0];
  resid[0] = -_ONE_;
  
  _PHIST_ERROR_HANDLER_(_SUBR_(mvec_delete)(vold,ierr),*ierr);
  _PHIST_ERROR_HANDLER_(_SUBR_(mvec_delete)(vnew,ierr),*ierr);

  return;
  }

