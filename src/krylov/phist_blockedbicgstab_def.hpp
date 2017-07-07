/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/

// implementation of gmres on several systems simultaneously
void SUBR(blockedBiCGStab_iterate)(TYPE(const_linearOp_ptr) Aop, TYPE(const_linearOp_ptr) Pop,
        TYPE(const_mvec_ptr) rhs, TYPE(mvec_ptr) sol_in, int numSys, int* nIter, _MT_ const tol[], int* iflag)
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

  PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),rhs,st::zero(),r,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_put_value)(sol,st::zero(),iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),r,st::zero(),p,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),r,st::zero(),r0,iflag),*iflag);

  std::vector<ST> rho(numSys,mt::one()), rho0(numSys,mt::one()), rho0_prev(numSys);

  //rho0_0 = r0^T r0
  PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(r0,r0,&rho0[0],iflag),*iflag);

  for (*nIter = 0; *nIter <= maxIter; (*nIter)++)
  {
    // rho_i = norm2(r_i)
    PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(r,r,&rho[0],iflag),*iflag);

    bool firstConverged = false;
    PHIST_SOUT(PHIST_VERBOSE,"BICGSTAB ITER %d: ",*nIter);
    for(int j = 0; j < numSys; j++)
    {
      rho[j]=std::sqrt(rho[j]);
      PHIST_SOUT(PHIST_VERBOSE,"\t%e",std::abs(rho[j]));
      firstConverged = firstConverged || std::abs(rho[j]) < tol[j];
    }
    PHIST_SOUT(PHIST_VERBOSE,"\n");
    if( firstConverged || *nIter == maxIter )
      break;
    

    // q_i = Op*p_i
    PHIST_CHK_IERR(Aop->apply(st::one(), Aop->A, p, st::zero(), q, iflag), *iflag);

    // alpha_i = rho0_i / (q_i^T r0)
    std::vector<ST> alpha(numSys);
    PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(q,r0,&alpha[0],iflag),*iflag);;
    for(int j = 0; j < numSys; j++)
      alpha[j] = rho0[j] / alpha[j];

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
      w[j] = (ST) ts[j] / w[j];   

    // x_(i+1) = x_i + alpha_i*p_i + w_i*s_i
    PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(&alpha[0],p,st::one(),sol,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(&w[0],s,st::one(),sol,iflag),*iflag);

    // r_(i+1) = s_i - w_i*t_i
    PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(&w[0],t,st::zero(),r,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),s,-st::one(),r,iflag),*iflag);

    // rho0_(i+1) = r_(i+1)^T r0*
    rho0_prev = rho0;
    PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(r,r0,&rho0[0],iflag),*iflag);

    // beta_i = (rho0_(i+1) / rho0_i) x (a_i / w_i)
    std::vector<ST> beta(numSys);
    for(int j = 0; j < numSys; j++)
      beta[j] = (ST) (rho0[j] / rho0_prev[j])*(alpha[j] / w[j]);    

    // p_(i+1) = r_(i+1) + beta_i*(p_i - w_i*q_i)
    for(int j = 0; j < numSys; j++)
      w[j] = -w[j];
    PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(&w[0],q,st::one(),p,iflag),*iflag);    
    PHIST_CHK_IERR(SUBR(mvec_vscale)(p,&beta[0],iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),r,st::one(),p,iflag),*iflag);
  }

}


