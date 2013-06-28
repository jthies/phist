// get type-sepcific macro wrappers for lapack,
// for instance XGEMM <- SGEMM, DGEMM, CGEMM or ZGEMM 
#include "phist_lapack_decl.h"

static int _SUBR_(converged)(_ST_ evmin)
{
	static _ST_ oldevmin = -1e9;

	int converged = abs(evmin-oldevmin) < 1e-9;
	oldevmin = evmin;

	return converged;
}


static void _SUBR_(lanczosStep)(_TYPE_(const_op_ptr) op, _TYPE_(mvec_ptr) vnew, _TYPE_(mvec_ptr) vold,
		_ST_ *alpha, _ST_ *beta, int* ierr)
{
        _ST_ ibeta;
        //TODO check these operations
        _PHIST_ERROR_HANDLER_(op->apply(_ONE_,op->A_,vold, -(*beta),vnew,ierr),*ierr);
        _PHIST_ERROR_HANDLER_(_SUBR_(mvec_dot_mvec)(vnew,vold,alpha,ierr),*ierr);
	_PHIST_ERROR_HANDLER_(_SUBR_(mvec_add_mvec)(-(*alpha),vold,_ONE_,vnew,ierr),*ierr);
	_PHIST_ERROR_HANDLER_(_SUBR_(mvec_dot_mvec)(vnew,vnew,beta,ierr),*ierr);
	*beta=_SQRT_(*beta);
	ibeta=_ONE_/(*beta);
	_PHIST_ERROR_HANDLER_(_SUBR_(mvec_scale)(vnew,&ibeta,ierr),*ierr);
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
  _TYPE_(mvec_ptr) vold, vnew, vtmp, r0;
  
  int nIter=*num_iters;
  int it,n;
  
  _ST_ alpha=_ZERO_, beta=_ZERO_;
  
  int ldz=1;
  _MT_ *work=NULL;
  _ST_ *x_ptr; // TODO - we don't compute eigenvectors yet
  
  if (op->range_map_ != op->domain_map_)
    {
    *ierr=-1; // maps of operator should point to same object (A should be a square matrix)
    }

  _PHIST_ERROR_HANDLER_(_SUBR_(mvec_create)(&vnew,op->domain_map_,1,ierr),*ierr);
  _PHIST_ERROR_HANDLER_(_SUBR_(mvec_create)(&vold,op->domain_map_,1,ierr),*ierr);
  _PHIST_ERROR_HANDLER_(_SUBR_(mvec_create)(&r0,op->domain_map_,1,ierr),*ierr);
  
  _PHIST_ERROR_HANDLER_(_SUBR_(mvec_random)(r0,ierr),*ierr);

  // normalize the global vector r0
  _MT_ rnorm;
  _PHIST_ERROR_HANDLER_(_SUBR_(mvec_normalize)(r0,&rnorm,ierr),*ierr); 

  _ST_ *alphas  = (_ST_ *)malloc(sizeof(_ST_)*nIter);
  _ST_ *betas   = (_ST_ *)malloc(sizeof(_ST_)*nIter);
  _ST_ *falphas = (_ST_ *)malloc(sizeof(_ST_)*nIter);
  _ST_ *fbetas  = (_ST_ *)malloc(sizeof(_ST_)*nIter);

  betas[0] = beta;

  for (it = 0, n=1; 
       it < nIter && !_SUBR_(converged)(falphas[0]); 
       it++, n++) 
	 {
//         printf("\r");
         printf("Iter %d ",it);
         _PHIST_ERROR_HANDLER_(_SUBR_(lanczosStep)(op,vnew,vold,&alpha,&beta,ierr),*ierr);
         vtmp=vnew;
         vnew=vold;
         vold=vtmp;

         alphas[it] = alpha;
         betas[it+1] = beta;
         memcpy(falphas,alphas,n*sizeof(_ST_)); // alphas and betas will be destroyed 
         memcpy(fbetas,betas,n*sizeof(_ST_));

         _PHIST_ERROR_HANDLER_(XSTEQR("N",&n,falphas,fbetas,x_ptr,&ldz,work,ierr),*ierr);
         printf("\n");
         }

  _PHIST_ERROR_HANDLER_(_SUBR_(mvec_delete)(r0,ierr),*ierr);
  _PHIST_ERROR_HANDLER_(_SUBR_(mvec_delete)(vold,ierr),*ierr);
  _PHIST_ERROR_HANDLER_(_SUBR_(mvec_delete)(vnew,ierr),*ierr);

  return;
  }

