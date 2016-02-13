// this implementation is adapted from the simple lanczos solver in essex/ghost-sandbox/ghost-apps/lanczos/
static int SUBR(converged)(_MT_ evmin)
{
#include "phist_std_typedefs.hpp"
    static _MT_ oldevmin = -1e9;

    int converged = mt::abs(evmin-oldevmin) < mt::sqrt(mt::eps())*0.1;
    oldevmin = evmin;

    return converged;
}

static void SUBR(lanczosStep)(TYPE(const_linearOp_ptr) A_op, TYPE(mvec_ptr) vnew, TYPE(const_mvec_ptr) vold,
        _ST_ *alpha, _MT_ *beta, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
// vnew = A*vold-beta*vnew
_ST_ minusbeta=-(_ST_)*beta;
PHIST_CHK_IERR(A_op->apply(st::one(),A_op->A, vold,minusbeta,vnew,iflag),*iflag);
//alpha=vnew'*vold
  PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(vnew,vold,alpha,iflag),*iflag);
  // vnew = vnew - alpha*vold
  PHIST_CHK_IERR(SUBR(mvec_add_mvec)(-*alpha,vold,st::one(),vnew,iflag),*iflag);
  // beta = ||vnew||_2
  PHIST_CHK_IERR(SUBR(mvec_norm2)(vnew,beta,iflag),*iflag);
  // vnew=vnew/beta
  _ST_ recbeta = st::one()/(*beta);
  PHIST_CHK_IERR(SUBR(mvec_scale)(vnew,recbeta,iflag),*iflag);
}


//! a simple Lanczos process to compute the largest eigenvalue and corresponding
//! eigenvector of a linear operator. This code is mainly intended as a testbed 
//! for fault tolerance (FT) capabilities in PHIST right now. Only block size 1 
//! is implemented.
//!
//! input:
//!
//! A: linearOp whose largest eigenpair is sought
//! *numIter, maximum number of iterations allowed
//!
//! output:
//!
//! *lambda_min, lambda_max, (approximate) smallest and largest eigenvalue
//! *numIter, number of iterations performed
//!
void SUBR(simple_lanczos)(TYPE(const_linearOp_ptr) A_op,
        _MT_ *evmin, _MT_ *evmax, int *numIter, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag = 0;
  
  int maxIter=*numIter;
  
  // create the vectors and arrays
  _MT_ alphas[maxIter],falphas[maxIter];
  _MT_ betas[maxIter],fbetas[maxIter];
  
  _ST_ alpha=st::zero();
  _MT_ beta=mt::zero();
  
  TYPE(mvec_ptr) vold, vnew;
  PHIST_CHK_IERR(SUBR(mvec_create)(&vold,A_op->domain_map,1,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_create)(&vnew,A_op->domain_map,1,iflag),*iflag);

  PHIST_CHK_IERR(SUBR(mvec_put_value)(vnew,st::zero(),iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_random)(vold,iflag),*iflag);
  _MT_ nrm;
  PHIST_CHK_IERR(SUBR(mvec_norm2)(vold,&nrm,iflag),*iflag);
  _ST_ scal=st::one()/nrm;
  PHIST_CHK_IERR(SUBR(mvec_scale)(vold,scal,iflag),*iflag);

  *evmin=mt::zero();

// put all iterations in one big compute task; this speeds up the tests with ghost (significantly)
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
  // Lanczos loop
  for(int i = 0; i < maxIter; i++)
  {
        int n=i+1;
        PHIST_CHK_IERR(SUBR(lanczosStep)(A_op,vnew,vold,&alpha,&beta,iflag),*iflag);
        std::swap(vold,vnew);

        alphas[i] = st::real(alpha);
        betas[i] = beta;
        for (int j=0;j<n;j++)
        {
          falphas[j]=alphas[j];
          fbetas[j]=betas[j];
        }

// compute eigenvalues of the tridiagonal matrix using LAPACK

// prohibit parallel execution to assure identical results on different procs
#pragma omp parallel
    {
#pragma omp master
      {
        st::blas_scalar_t z;
        int ldz=1;
        PHIST_TG_PREFIX(STEQR)((blas_char_t*)"N", &n, falphas, fbetas, &z, &ldz, NULL, iflag );
      }
    }
        if (*iflag)
        {
          PHIST_OUT(PHIST_ERROR,"Error %d returned from lapack call XSTEQR!\n",*iflag);
          return;
        }
        *evmin=falphas[0];
        *evmax=falphas[n-1];
        PHIST_SOUT(PHIST_INFO,"LANCZOS step %d\tmin/max eigenvalue: %f/%f\n", i,*evmin,*evmax);
#if 0
         PHIST_SOUT(PHIST_INFO,"D=[ ");
        for (int j=0;j<n;j++) PHIST_SOUT(PHIST_INFO,"%25.16e ",betas[j]);
        PHIST_SOUT(PHIST_INFO,"];\nE=[ ");
        for (int j=0;j<n;j++) PHIST_SOUT(PHIST_INFO,"%25.16e ",alphas[j]);
        PHIST_SOUT(PHIST_INFO,"];\n");
        PHIST_SOUT(PHIST_INFO,"A=spdiags([E',D',[0,E(1:end-1)]'],-1:1,%d,%d)\n",n,n);
#endif       
          
        *numIter=i;
        if (SUBR(converged)(*evmin)) break;
  }
  

PHIST_TASK_END(iflag)


  // delete vectors
  PHIST_CHK_IERR(SUBR(mvec_delete)(vold,  iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_delete)(vnew, iflag), *iflag);
}
