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
        _MT_* evals, _MT_* resid,
        _MT_ tol,int* num_iters, int* num_eigs,
        int* ierr)
  {
  _TYPE_(mvec_ptr) vold, vnew, vtmp;
  _TYPE_(mvec_ptr) V;
  _TYPE_(sdMat_ptr) S;
  
  int nIter=*num_iters;
  int nconv=0;
  int i,it,n, inext;
  _ST_ alpha, beta;
  int converged[nIter];
  _MT_ r_est[nIter];
  _MT_ nrm;
 
   // allocate memory for the tridiagonal matrix H = [-beta_(j-1) alpha_j -beta_j]
  _MT_ *alphas  = (_MT_ *)malloc(sizeof(_MT_)*nIter);
  _MT_ *betas   = (_MT_ *)malloc(sizeof(_MT_)*nIter);
  _MT_ *falphas = (_MT_ *)malloc(sizeof(_MT_)*nIter);
  _MT_ *fbetas  = (_MT_ *)malloc(sizeof(_MT_)*nIter);
  _MT_ *work  = (_MT_ *)malloc(sizeof(_MT_)*(2*nIter-2));
  for (i=0;i<nIter;i++) converged[i]=0; 
  int lds;
  _MT_ *S_ptr; // will point to the Ritz values
  //int i,lda,nloc;
  *num_iters=0;
  
  if (op->range_map_ != op->domain_map_)
    {
    *ierr=-1; // maps of operator should point to same object (A should be a square matrix)
    }

  _PHIST_ERROR_HANDLER_(_SUBR_(mvec_create)(&vold,op->domain_map_,1,ierr),*ierr);
  _PHIST_ERROR_HANDLER_(_SUBR_(mvec_create)(&vnew,op->domain_map_,1,ierr),*ierr);
  // S will store the Ritz vectors
#ifdef _IS_DOUBLE_  
  _PHIST_ERROR_HANDLER_(phist_DsdMat_create(&S,nIter,nIter,ierr),*ierr);
  // pointer to the data in S
  _PHIST_ERROR_HANDLER_(phist_DsdMat_extract_view(S,&S_ptr, &lds,ierr),*ierr);
#else
  _PHIST_ERROR_HANDLER_(phist_SsdMat_create(&S,nIter,nIter,ierr),*ierr);
  // pointer to the data in S
  _PHIST_ERROR_HANDLER_(phist_SsdMat_extract_view(S,&S_ptr, &lds,ierr),*ierr);
#endif
  
  // vold = random start vector, r=vold
  _PHIST_ERROR_HANDLER_(_SUBR_(mvec_put_value)(vold,_ZERO_,ierr),*ierr);
  _PHIST_ERROR_HANDLER_(_SUBR_(mvec_random)(vnew,ierr),*ierr);
  // TROET
  _PHIST_ERROR_HANDLER_(_SUBR_(mvec_put_value)(vnew,_ONE_,ierr),*ierr);

  _PHIST_ERROR_HANDLER_(_SUBR_(mvec_normalize)(vnew,&nrm,ierr),*ierr);

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

        if ((it+1)%5==0)
          {
          memcpy(falphas,alphas,n*sizeof(_MT_)); // alphas and betas will be destroyed 
          memcpy(fbetas,betas,n*sizeof(_MT_));
#ifdef _IS_DOUBLE_
          _PHIST_ERROR_HANDLER_(DSTEQR("I",&n,falphas,fbetas+1,S_ptr,&lds,work,ierr),*ierr);
#else
          _PHIST_ERROR_HANDLER_(SSTEQR("I",&n,falphas,fbetas+1,S_ptr,&lds,work,ierr),*ierr);
#endif
           // bound for ||r|| of the first (smallest) Ritz value.
           // Note that we check for the smallest one because the
           // actual eigenvalues of A have a reversed sign.
           nrm =1e99;
           inext=0;
           nconv=0;
           for (i=0;i<n;i++)
             {
//             fprintf(stdout,"S[%d] %4.2e %4.2e %4.2e\n",i,S_ptr[i*lds],S_ptr[i*lds+1],S_ptr[i*lds+2]);
             r_est[i]=fabs(betas[it+1]*S_ptr[i*lds+it]);
//             fprintf(stdout,"lmb[%d] = %e, ||r|| approx. %e\n",i,falphas[i],r_est[i]);
             if (r_est[i]<tol)
               {
               converged[i]=1;
               nconv++;
               }
             else if (r_est[i]<nrm)
               {
               nrm=r_est[i];
               inext=i;
               }
             }
//           printf("it %d, %d converged, next eig to converge: %e\t||r|| est: %e\n",n,nconv,falphas[inext],r_est[inext]);
           if (nconv>=*num_eigs) break;
           }
         }
         
  *num_iters=it;
  inext=0;
  // copy at most num_eigs converged eigenvalues into the user
  // provided array
  for (i=n-1;i>=0, inext<*num_eigs;i--)
    {
    if (converged[i])
      {
      evals[inext] = falphas[i];
      resid[inext] = r_est[i];
      inext++;
      }
    }
  if (nconv<*num_eigs) *num_eigs=nconv;
  
  // TODO - compute eigenvectors (need to store the whole basis V for that)
  
  _PHIST_ERROR_HANDLER_(_SUBR_(mvec_delete)(vold,ierr),*ierr);
  _PHIST_ERROR_HANDLER_(_SUBR_(mvec_delete)(vnew,ierr),*ierr);
  _PHIST_ERROR_HANDLER_(_SUBR_(sdMat_delete)(S,ierr),*ierr);
  
  free(alphas);
  free(betas);
  free(falphas);
  free(fbetas);
  free(work);

  return;
  }

