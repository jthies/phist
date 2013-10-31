// on entry, vnew = r, vold=v_(j-1).
// on exit, vnew=v_j, vold=r
static void SUBR(lanczosStep)(TYPE(const_op_ptr) op, 
                TYPE(mvec_ptr) vnew, 
                TYPE(mvec_ptr) vold,
		_ST_ *alpha, _ST_ *beta, int* ierr)
{
#include "phist_std_typedefs.hpp"

        *ierr=0;
        // vnew = r/beta
        ST ibeta = st::one()/(*beta);
        PHIST_CHK_IERR(SUBR(mvec_scale)(vnew,ibeta,ierr),*ierr);
        // r = A*v - beta*vold 
        PHIST_CHK_IERR(op->apply(st::one(),op->A,vnew, 
                -(*beta),vold,ierr),*ierr);
        //alpha_j = v_j'r
        PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(vnew,vold,alpha,ierr),*ierr);
        // r=r-alpha_j*v_j
	PHIST_CHK_IERR(SUBR(mvec_add_mvec)(-(*alpha),vnew,st::one(),vold,ierr),*ierr);
	// beta_j = ||r||_2
	PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(vold,vold,beta,ierr),*ierr);
	*beta=st::sqrt(*beta);
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
void SUBR(lanczos)(TYPE(const_op_ptr) op, 
        TYPE(mvec_ptr) X,
        _MT_* evals, _MT_* resid,
        _MT_ tol,int* num_iters, int* num_eigs,
        int* ierr)
  {
#include "phist_std_typedefs.hpp"

  st::mvec_t *vold, *vnew, *vtmp;
  st::mvec_t *V;
  mt::sdMat_t *S;
  
  int nIter=*num_iters;
  int nconv=0;
  int i,it,n, inext;
  ST alpha, beta;
  int converged[nIter];
  MT r_est[nIter];
  MT nrm;
 
   // allocate memory for the tridiagonal matrix H = [-beta_(j-1) alpha_j -beta_j]
  MT *alphas  = (MT *)malloc(sizeof(MT)*nIter);
  MT *betas   = (MT *)malloc(sizeof(MT)*nIter);
  MT *falphas = (MT *)malloc(sizeof(MT)*nIter);
  MT *fbetas  = (MT *)malloc(sizeof(MT)*nIter);
  MT *work  = (MT *)malloc(sizeof(MT)*(2*nIter-2));
  for (i=0;i<nIter;i++) converged[i]=0; 
  int lds;
  MT *S_ptr_t; // will point to the Ritz values
  //int i,lda,nloc;
  *num_iters=0;
  
  if (op->range_map != op->domain_map || op->range_map==NULL)
    {
    *ierr=-1; // maps of operator should point to same object (A should be a square matrix)
              // TODO - we can't implement this easily in ghost, maybe the test should be
              // skipped or a function added to the kernel lib to compare maps.
    }
    
  const_comm_ptr_t comm;
  PHIST_CHK_IERR(phist_map_get_comm(op->range_map,&comm,ierr),*ierr);

  PHIST_OUT(1,"create vectors");
  PHIST_CHK_IERR(SUBR(mvec_create)(&vold,op->domain_map,1,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_create)(&vnew,op->domain_map,1,ierr),*ierr);

  PHIST_OUT(3,"vold @ %p\n",vold);
  PHIST_OUT(3,"vnew @ %p\n",vnew);

  // S will store the Ritz vectors
  // it would be nice if we could do these things with a traits class eventually,
  // rather than with macros.
  PHIST_OUT(1,"create matrix for ritz vectors");
#ifdef _IS_DOUBLE_
  PHIST_CHK_IERR(phist_DsdMat_create(&S,nIter,nIter,comm,ierr),*ierr);
  // pointer to the data in S
  PHIST_CHK_IERR(phist_DsdMat_extract_view(S,&S_ptr_t, &lds,ierr),*ierr);
#else
  PHIST_CHK_IERR(phist_SsdMat_create(&S,nIter,nIter,comm,ierr),*ierr);
  // pointer to the data in S
  PHIST_CHK_IERR(phist_SsdMat_extract_view(S,&S_ptr_t, &lds,ierr),*ierr);
#endif
  
  // vold = random start vector, r=vold
  PHIST_OUT(1,"initialize vectors");
  PHIST_CHK_IERR(SUBR(mvec_put_value)(vold,st::zero(),ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_random)(vnew,ierr),*ierr);

  PHIST_OUT(1,"normalize start vector");
  PHIST_CHK_IERR(SUBR(mvec_normalize)(vnew,&nrm,ierr),*ierr);

  beta = st::one();
  betas[0] = st::real(beta);

  for(it = 0, n=1; 
      it < nIter;
      it++, n++)
    {
         PHIST_CHK_IERR(SUBR(lanczosStep)(op,vnew,vold,&alpha,&beta,ierr),*ierr);
         vtmp=vnew;
         vnew=vold;
         vold=vtmp;
         alphas[it]=st::real(alpha);
         betas[it+1]=st::real(beta);

        if ((it+1)%5==0)
          {
          memcpy(falphas,alphas,n*sizeof(MT)); // alphas and betas will be destroyed 
          memcpy(fbetas,betas,n*sizeof(MT));
#ifdef _IS_DOUBLE_
          PHIST_CHK_IERR(DSTEQR("I",&n,falphas,fbetas+1,S_ptr_t,&lds,work,ierr),*ierr);
#else
          PHIST_CHK_IERR(SSTEQR("I",&n,falphas,fbetas+1,S_ptr_t,&lds,work,ierr),*ierr);
#endif
           // bound for ||r|| of the first (smallest) Ritz value.
           // Note that we check for the smallest one because the
           // actual eigenvalues of A have a reversed sign.
           nrm =10000.;
           inext=0;
           nconv=0;
           for (i=0;i<n;i++)
             {
             //fprintf(stdout,"S[%d] %4.2e %4.2e %4.2e\n",i,S_ptr_t[i*lds],S_ptr_t[i*lds+1],S_ptr_t[i*lds+2]);
             r_est[i]=fabs(betas[it+1]*S_ptr_t[i*lds+it]);
             //fprintf(stdout,"lmb[%d] = %e, ||r|| approx. %e\n",i,falphas[i],r_est[i]);
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
           PHIST_OUT(1,"it %d, %d converged, next eig to converge: %e\t||r|| est: %e\n",
                n, nconv, falphas[inext], r_est[inext]);
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
  
  PHIST_CHK_IERR(SUBR(mvec_delete)(vold,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_delete)(vnew,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(S,ierr),*ierr);
  
  free(alphas);
  free(betas);
  free(falphas);
  free(fbetas);
  free(work);

  return;
  }

