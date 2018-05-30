/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/

// implementation of BiCGStab on several systems simultaneously
extern "C" void SUBR(blockedBiCGStab_iterate)(TYPE(const_linearOp_ptr) Aop, TYPE(const_linearOp_ptr) Pop,
        TYPE(const_mvec_ptr) rhs, TYPE(mvec_ptr) sol_in, TYPE(const_mvec_ptr) V,
        int numSys, int* nIter, _MT_ const tol[], int* iflag)
{
#include "phist_std_typedefs.hpp"
  *iflag = 0;
  if (numSys==0) return; // do not appear in timing stats
  PHIST_ENTER_FCN(__FUNCTION__);

  int maxIter = (*nIter)>0 ? *nIter: 9999999;

  *iflag=0;

  if( Aop == NULL || rhs == NULL || sol_in == NULL )
  {
    *iflag = PHIST_INVALID_INPUT;
    return;
  }
  
  PHIST_CHK_IERR(*iflag=(Pop==NULL)?0:PHIST_NOT_IMPLEMENTED,*iflag);
  {
  }

  // check dimensions
  {
    int ncrhs,ncsol;
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(rhs,&ncrhs,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(sol_in,&ncsol,iflag),*iflag);
    // note: we actually don't care about the dimensions of sol as long as it is large 
    // enough since we just work with views in this function. However, we issue a warning 
    // (positive iflag) if the array dimensions are larger than expected so that the user 
    // doesn't unexpectedly end up with empty rows/cols
    if (ncsol<numSys || ncrhs<numSys) {*iflag = -1; return;}
    if (ncsol>numSys || ncrhs>numSys)
    {
      PHIST_SOUT(PHIST_VERBOSE,"REMARK: input vectors/matrix to BiCGStab are larger than necessary.\n");
      PHIST_SOUT(PHIST_VERBOSE,"        sol has %d columns (expecting %d)\n",
                                        ncsol, ncrhs);
    }
  }

  // allocate temporary storage
  TYPE(mvec_ptr) p=NULL, q=NULL, r=NULL, r0=NULL, s=NULL, t=NULL;
  PHIST_CHK_IERR(SUBR(mvec_create) (&p,  Aop->domain_map, numSys, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_create) (&q,  Aop->domain_map, numSys, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_create) (&r,  Aop->domain_map, numSys, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_create) (&r0, Aop->domain_map, numSys, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_create) (&s, Aop->domain_map, numSys, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_create) (&t, Aop->domain_map, numSys, iflag), *iflag);

  TYPE(mvec_ptr) sol=NULL;
  PHIST_CHK_IERR(SUBR(mvec_view_block)(sol_in,&sol,0,numSys-1,iflag),*iflag);
  MvecOwner<_ST_> _sol(sol), _p(p),_q(q),_r(r),_r0(r0),_s(s),_t(t);

  // set x=0, p=r=r0=b
  PHIST_CHK_IERR(SUBR(mvec_put_value)(sol,st::zero(),iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_get_block)(rhs,r,0,numSys-1,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),r,st::zero(),p,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),r,st::zero(),r0,iflag),*iflag);

  // rho0=(r0*r0), rho=(r*r0), rr=(r*r)
  std::vector<ST> rho(numSys,mt::one()), rho0(numSys,mt::one()), rho_prev(numSys),rr(numSys,mt::one());

  PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(r0,r0,&rho0[0],iflag),*iflag);
  rho=rho0;

  for (*nIter = 0; *nIter <= maxIter; (*nIter)++)
  {
    bool firstConverged = false;
    PHIST_SOUT(PHIST_VERBOSE,"BICGSTAB ITER %d:",*nIter);
    for(int j = 0; j < numSys; j++)
    {
      // we can look at diffrent criterias for convergence
      // this two here get similar results in most cases
      // in cases when BiCGStab needs very many iteration steps to decrease the residuum,
      // the first criteria stops early but with a bigger residuum than the given toleranz
      
      MT rnrm;
      if( maxIter < 40 )
      {
        // this criteria yields to early termination with a bigger residuum
        rnrm=std::sqrt(std::abs(rho[j]/rho0[j]));
      }
      else
      {
        // this criteria yields to termination if the residuum is smal enough relative to rhs
        // (can take many iteration steps)
        PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(r,r,&rr[0],iflag),*iflag);
        rnrm=std::sqrt(std::abs(rr[j]/rho0[j]));
      }
	  
      PHIST_SOUT(PHIST_VERBOSE,"\t%e",rnrm);
      firstConverged = firstConverged || (rnrm < tol[j]);
    }
    PHIST_SOUT(PHIST_VERBOSE,"\n");
    if( firstConverged || *nIter == maxIter )
      break;
    

    // q_i = Op*p_i
    PHIST_CHK_IERR(Aop->apply(st::one(), Aop->A, p, st::zero(), q, iflag), *iflag);

    // alpha_i = rho_i / (q_i^T r0)
    std::vector<ST> alpha(numSys);
    PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(q,r0,&alpha[0],iflag),*iflag);;
    for(int j = 0; j < numSys; j++)
    {
      alpha[j] = rho[j] / alpha[j];
    }

    // s_i = r_i - alpha_i*q_i
    PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(&alpha[0],q,st::zero(),s,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),r,-st::one(),s,iflag),*iflag); 

    // t_i = A*s_i
    PHIST_CHK_IERR(Aop->apply(st::one(), Aop->A, s, st::zero(), t, iflag), *iflag);

    // w_i = (t_i^T s_i) / (t_i^T t_i)
    std::vector<ST> w(numSys);
    std::vector<ST> ts(numSys);
    PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(t,s,&ts[0],iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(t,t,&w[0],iflag),*iflag);
    for(int j = 0; j < numSys; j++)
    {
      w[j] = ts[j] / w[j];   
    }

    // TODO: rprovide a kernel for these 3-term recurrences!

    // x_(i+1) = x_i + alpha_i*p_i + w_i*s_i
    PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(&alpha[0],p,st::one(),sol,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(&w[0],s,st::one(),sol,iflag),*iflag);

    // r_(i+1) = s_i - w_i*t_i
    PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(&w[0],t,st::zero(),r,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),s,-st::one(),r,iflag),*iflag);

    // rho_(i+1) = r_(i+1)^T*r0
    rho_prev = rho;
    PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(r,r0,&rho[0],iflag),*iflag);

    // beta_i = (rho0_(i+1) / rho0_i) x (a_i / w_i)
    std::vector<ST> beta(numSys);
    for(int j = 0; j < numSys; j++)
    {
      beta[j] = (rho[j] / rho_prev[j])*(alpha[j] / w[j]);    
    }
    // p_(i+1) = r_(i+1) + beta_i*(p_i - w_i*q_i)
    // TODO: can we express this in some useful kernel, see comment above
    for(int j = 0; j < numSys; j++) w[j] = -w[j];
    PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(&w[0],q,st::one(),p,iflag),*iflag);    
    PHIST_CHK_IERR(SUBR(mvec_vscale)(p,&beta[0],iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),r,st::one(),p,iflag),*iflag);
  }
}

extern "C" void SUBR( BiCGStab ) (TYPE(const_linearOp_ptr) Aop, TYPE(const_linearOp_ptr) Pop,
		TYPE(const_mvec_ptr) rhs, TYPE(mvec_ptr) sol_in, int* nIter, _MT_ const tol, int* iflag)
{
#include "phist_std_typedefs.hpp"
  *iflag = 0;
  PHIST_ENTER_FCN(__FUNCTION__);
  
  PHIST_CHK_IERR(SUBR(blockedBiCGStab_iterate)(Aop, Pop, rhs, sol_in, NULL, 1, nIter, &tol, iflag),*iflag);

}		


