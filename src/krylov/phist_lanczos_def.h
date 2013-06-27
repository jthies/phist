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
  _TYPE_(mvec_ptr) vold, vnew, r0;

  int nIter=*num_iters;  
  int iteration;
  
  if (op->range_map_ != op->domain_map_)
    {
    *ierr=-1; // maps of operator should point to same object (A should be a square matrix)
    }

  _PHIST_ERROR_HANDLER_(_SUBR_(mvec_create)(&vnew,op->domain_map_,1,ierr),*ierr);
  _PHIST_ERROR_HANDLER_(_SUBR_(mvec_create)(&r0,op->domain_map_,1,ierr),*ierr);

/*
	ghost_normalizeVector(r0); // normalize the global vector r0

	vold = ghost_cloneVector(r0); 

	ghost_mdat_t *alphas  = (ghost_mdat_t *)malloc(sizeof(ghost_mdat_t)*nIter);
	ghost_mdat_t *betas   = (ghost_mdat_t *)malloc(sizeof(ghost_mdat_t)*nIter);
	ghost_mdat_t *falphas = (ghost_mdat_t *)malloc(sizeof(ghost_mdat_t)*nIter);
	ghost_mdat_t *fbetas  = (ghost_mdat_t *)malloc(sizeof(ghost_mdat_t)*nIter);

	betas[0] = beta;

	for(iteration = 0, n=1; 
			iteration < nIter && !converged(falphas[0]); 
			iteration++, n++) 
	{
		printf("\r");

		lanczosStep(context,vnew,vold,&alpha,&beta);
		ghost_swapVectors(vnew,vold);

		alphas[iteration] = alpha;
		betas[iteration+1] = beta;
		memcpy(falphas,alphas,n*sizeof(ghost_mdat_t)); // alphas and betas will be destroyed in imtql
		memcpy(fbetas,betas,n*sizeof(ghost_mdat_t));

		// TODO evaluate headers in fortran files
		if (GHOST_MY_MDATATYPE & GHOST_BINCRS_DT_FLOAT) {
			imtql1f_(&n,falphas,fbetas,&ferr);
		} else
			imtql1_(&n,falphas,fbetas,&ferr);

		if(ferr != 0) printf("Error: the %d. ev could not be determined\n",ferr);
		printf("minimal eigenvalue: %f", falphas[0]);
		fflush(stdout);
	}
	printf("\n");

	ghost_freeVector(r0);
	ghost_freeVector(vold);
	ghost_freeVector(vnew);
	ghost_freeContext (context);

	ghost_finish();
*/
  return;
  }
/*
static int converged(ghost_mdat_t evmin)
{
	static ghost_mdat_t oldevmin = -1e9;

	int converged = MABS(evmin-oldevmin) < 1e-9;
	oldevmin = evmin;

	return converged;
}


static void lanczosStep(ghost_context_t *context, ghost_vec_t *vnew, ghost_vec_t *vold,
		ghost_mdat_t *alpha, ghost_mdat_t *beta)
{
	vecscal(vnew,-*beta,context->lnrows(context));
	ghost_spmvm(vnew, context, vold, GHOST_SPMVM_MODE_NOMPI|GHOST_SPMVM_AXPY);
	dotprod(vnew,vold,alpha,context->lnrows(context));
	axpy(vnew,vold,-(*alpha),context->lnrows(context));
	dotprod(vnew,vnew,beta,context->lnrows(context));
	*beta=MSQRT(*beta);
	vecscal(vnew,(ghost_mdat_t)1./(*beta),context->lnrows(context));
}
*/
