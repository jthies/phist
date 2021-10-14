/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/

// implementation of CG on several systems simultaneously
extern "C" void SUBR(blockedPCG)(TYPE(const_linearOp_ptr) Aop, TYPE(const_linearOp_ptr) Pop,
        TYPE(const_mvec_ptr) rhs, TYPE(mvec_ptr) sol_in, int numSys, int* nIter,
        _MT_ const tol[], int* iflag)
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
  TYPE(mvec_ptr) p=NULL, q=NULL, r=NULL, r0=NULL, z=NULL;
  PHIST_CHK_IERR(SUBR(mvec_create)(&p,Aop->domain_map,numSys,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_create)(&q,Aop->domain_map,numSys,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_create)(&r,Aop->domain_map,numSys,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_create)(&r0,Aop->domain_map,numSys,iflag),*iflag);
  if (Pop!=NULL)
  {
    PHIST_CHK_IERR(SUBR(mvec_create)(&z,Aop->domain_map,numSys,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_put_value)(z,_ST_(0),iflag),*iflag);
  }
  else
  {
    z=r;
  }
  TYPE(mvec_ptr) sol=NULL;
  PHIST_CHK_IERR(SUBR(mvec_view_block)(sol_in,&sol,0,numSys-1,iflag),*iflag);
  MvecOwner<_ST_> _sol(sol), _p(p),_q(q),_r(r),_r0(r0),_z(z!=r?z:NULL);

  // set x=p=0, r=r0=b, z=P*r
  PHIST_CHK_IERR(SUBR(mvec_put_value)(sol,st::zero(),iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),rhs,st::zero(),r,iflag),*iflag);
  if (Pop!=NULL)
  {
    PHIST_CHK_IERR(Pop->apply(st::one(), Pop->A, r, st::zero(), z, iflag), *iflag);
  }
  PHIST_CHK_IERR(SUBR(mvec_put_value)(p,st::zero(),iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),r,st::zero(),r0,iflag),*iflag);
  
  // rho0=(r0*r0), rho=(r*z) ,rr=(r*r)
  std::vector<ST> rho(numSys,mt::one()),rho0(numSys,mt::one()), rho_prev(numSys),rr(numSys,mt::one());
  PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(r0,r0,&rho0[0],iflag),*iflag);

  for(*nIter = 0; *nIter <= maxIter; (*nIter)++)
  {
    // rho_i = r_i^T z_i
    rho_prev = rho;
    PHIST_SOUT(PHIST_INFO,"CG iteration %d, est. residual norm:",*nIter);
    bool firstConverged = false;
    for(int j = 0; j < numSys; j++)
    {
      PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(r,r,&rr[0],iflag),*iflag);
      MT rnrm=std::sqrt(std::abs(rr[j]/rho0[j]));
      PHIST_SOUT(PHIST_INFO,"\t%e",rnrm);
      firstConverged = firstConverged || (rnrm < tol[j]);
    }
    PHIST_SOUT(PHIST_INFO,"\n");
    if( firstConverged || *nIter == maxIter ) break;
    // rho_i = r_i*z_i
    PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(r,z,&rho[0],iflag),*iflag);

    // beta_i = rho_i / rho_(i-1)
    std::vector<ST> beta(numSys);
    for(int j = 0; j < numSys; j++)
      beta[j] = rho[j] / rho_prev[j];
  
    // p_(i+1) = z_i + beta_i*p_i
    PHIST_CHK_IERR(SUBR(mvec_vscale)(p,&beta[0],iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),z,st::one(),p,iflag),*iflag);

    // q_(i+1) = A*p_(i+1)
    PHIST_CHK_IERR(Aop->apply(st::one(), Aop->A, p, st::zero(), q, iflag), *iflag);

    // alpha_(i+1) = rho_i / p_(i+1)^T q_(i+1)
    std::vector<ST> alpha(numSys);
    PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(p,q,&alpha[0],iflag),*iflag);
    for(int j = 0; j < numSys; j++)
      alpha[j] = rho[j] / alpha[j];

    // x_(i+1) = x_i + alpha*p_(i+1)
    PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(&alpha[0],p,st::one(),sol,iflag),*iflag);
    // r_(i+1) = r_i - alpha*q_(i+1)
    for(int j = 0; j < numSys; j++)
    {
      alpha[j] = -alpha[j];
    }
    PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(&alpha[0],q,st::one(),r,iflag),*iflag);

    // z_(i+1) = P_op*r_(i+1)
    if (Pop!=NULL)
    {
      PHIST_CHK_IERR(Pop->apply(st::one(), Pop->A, r, st::zero(), z, iflag), *iflag);
    }
  }
}

extern "C" void SUBR( PCG ) (TYPE(const_linearOp_ptr) Aop, TYPE(const_linearOp_ptr) Pop,
        TYPE(const_mvec_ptr) rhs, TYPE(mvec_ptr) sol_in, int* nIter, _MT_ const tol, int* iflag)
{
#include "phist_std_typedefs.hpp"
  *iflag = 0;
  PHIST_ENTER_FCN(__FUNCTION__);
  
  int num_sol,num_rhs;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(rhs,&num_rhs,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(sol_in,&num_sol,iflag),*iflag);
  PHIST_CHK_IERR(num_sol==num_rhs?0: PHIST_INVALID_INPUT, *iflag);
  _MT_ vtol[num_rhs];
  for (int i=0; i<num_rhs; i++) vtol[i]=tol;
  PHIST_CHK_IERR(SUBR(blockedPCG)(Aop, Pop, rhs, sol_in, 1, nIter, &tol, iflag),*iflag);

}


