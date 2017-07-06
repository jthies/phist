/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/

// implementation of gmres on several systems simultaneously
void SUBR(blockedQMR_iterate)(TYPE(const_linearOp_ptr) Aop, TYPE(const_linearOp_ptr) Pop,
        TYPE(const_mvec_ptr) rhs, TYPE(mvec_ptr) sol, int numSys, int* nIter, _MT_ const tol[], int* iflag)
{
#include "phist_std_typedefs.hpp"
  *iflag = 0;
  if (numSys==0) return; // do not appear in timing stats
  PHIST_ENTER_FCN(__FUNCTION__);

  int maxIter = (*nIter)>0 ? *nIter: 9999999;

  *iflag=0;

  if( Aop == NULL && rhs == NULL && sol == NULL )
  {
    *iflag = PHIST_NOT_IMPLEMENTED;
    return;
  }

  // check dimensions
  {
    int ncrhs,ncsol;
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(rhs,&ncrhs,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(sol,&ncsol,iflag),*iflag);
    // note: we actually don't care about the dimensions of sol as long as it is large 
    // enough since we just work with views in this function. However, we issue a warning 
    // (positive iflag) if the array dimensions are larger than expected so that the user 
    // doesn't unexpectedly end up with empty rows/cols
    if (ncsol<numSys || ncrhs<numSys) {*iflag = -1; return;}
    if (ncsol>numSys || ncrhs>numSys)
    {
      PHIST_SOUT(PHIST_VERBOSE,"REMARK: input vectors/matrix to qmr are larger than necessary.\n");
      PHIST_SOUT(PHIST_VERBOSE,"        sol has %d columns (expecting %d)\n",
                                        ncsol, ncrhs);
    }
  }

  // allocate temporary storage
  TYPE(mvec_ptr) r = NULL, rp = NULL, v = NULL, t = NULL, q = NULL, p = NULL, u = NULL, d = NULL;
  PHIST_CHK_IERR(SUBR(mvec_create) (&r,  Aop->domain_map, numSys, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_create) (&rp, Aop->domain_map, numSys, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_create) (&v,  Aop->domain_map, numSys, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_create) (&t,  Aop->domain_map, numSys, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_create) (&q,  Aop->domain_map, numSys, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_create) (&p,  Aop->domain_map, numSys, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_create) (&u,  Aop->domain_map, numSys, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_create) (&d,  Aop->domain_map, numSys, iflag), *iflag);

  // compute initial residual
  // assuming x0 = 0, r = rhs
  PHIST_CHK_IERR(SUBR(mvec_set_block)(r, rhs, 0, numSys-1, iflag), *iflag);

  // dp = norm(R)
  _MT_ dp[numSys];
  PHIST_CHK_IERR(SUBR(mvec_norm2)(r, dp, iflag), *iflag);

  // Rp = R  
  PHIST_CHK_IERR(SUBR(mvec_set_block)(rp, r, 0, numSys-1, iflag), *iflag);

  // Set initial conditions
  _ST_ etaold[numSys], psiold[numSys], rhoold[numSys];
  _MT_ tau[numSys], dpold[numSys];
  for (int i=0; i<numSys; i++) {
    etaold[i] = 0.0;
    psiold[i] = 0.0;
    tau   [i] = dp[i];
    dpold [i] = dp[i];
  }

  PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(r, rp, rhoold, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_set_block)(u, r, 0, numSys-1, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_set_block)(p, r, 0, numSys-1, iflag), *iflag);

  // v = AKp
  if (Pop) {
    PHIST_CHK_IERR(Aop->apply(st::one(), Aop->A, p, st::zero(), t, iflag), *iflag);
    PHIST_CHK_IERR(Pop->apply(st::one(), Pop->A, t, st::zero(), v, iflag), *iflag);
  } else {
    PHIST_CHK_IERR(Aop->apply(st::one(), Aop->A, p, st::zero(), v, iflag), *iflag);
  }

  PHIST_CHK_IERR(SUBR(mvec_put_value)(d, st::zero(), iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_put_value)(sol, st::zero(), iflag), *iflag);

  bool conv = false; // true if some rhs has converged
  for (*nIter = 0; *nIter < maxIter; ++*nIter) {
    // s = v'*rp
    _ST_ s[numSys];
    PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(v, rp, s, iflag), *iflag);

    // a = rhoold / s
    _ST_ a[numSys];
    for (int i=0; i<numSys; i++) a[i] = rhoold[i] / s[i];

    // q = -a*v + u
    PHIST_CHK_IERR(SUBR(mvec_set_block)(q, u, 0, numSys-1, iflag), *iflag);
    for (int i=0; i<numSys; i++) a[i] *= -1.0;
    PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(a, v, st::one(), q, iflag), *iflag);
    for (int i=0; i<numSys; i++) a[i] *= -1.0;

    // r = r - a A K (u + q)
    PHIST_CHK_IERR(SUBR(mvec_set_block)(t, u, 0, numSys-1, iflag), *iflag);
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(), q, st::one(), t, iflag), *iflag);
    for (int i=0; i<numSys; i++) a[i] *= -1.0;
    if (Pop) {
      PHIST_CHK_IERR(Aop->apply(st::one(), Aop->A, t, st::zero(), v, iflag), *iflag);
      PHIST_CHK_IERR(Pop->apply(st::one(), Pop->A, v, st::zero(), t, iflag), *iflag);
      PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(a, t, st::one(), r, iflag), *iflag);
    } else {
      PHIST_CHK_IERR(Aop->apply(st::one(), Aop->A, t, st::zero(), v, iflag), *iflag);
      PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(a, v, st::one(), r, iflag), *iflag);
    }
    for (int i=0; i<numSys; i++) a[i] *= -1.0;

    // dp = norm(R)
    _MT_ dp[numSys];
    PHIST_CHK_IERR(SUBR(mvec_norm2)(r, dp, iflag), *iflag);

    for (int m=0; m<2; m++) {
      _ST_ eta[numSys], cf[numSys];
      _MT_ w, psi[numSys], cm;
      for (int i=0; i<numSys; i++) {
        if (!m) w = sqrt(dp[i]*dpold[i]);
        else w = dp[i];
        psi[i] = w / tau[i];
        cm  = 1.0 / sqrt(1.0 + psi[i] * psi[i]);
        tau[i] = tau[i] * psi[i] * cm;
        eta[i] = cm * cm * a[i];
        cf[i]  = psiold[i] * psiold[i] * etaold[i] / a[i];
      }

      // D = cf D + u if m == 0, else cf D + q
      PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(cf, d, st::one(), d, iflag), *iflag);
      PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(), m==0?u:q, st::one(), d, iflag), *iflag);

      // sol = sol + eta d
      PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(eta, d, st::one(), sol, iflag), *iflag);

      PHIST_SOUT(PHIST_VERBOSE,"QMR ITER %d: ",*nIter);
      
      for (int i=0; i<numSys; i++) {
        _MT_ dpest = sqrt(m + 1.0) * tau[i];
             PHIST_SOUT(PHIST_VERBOSE,"\t%8.4g",dpest);

        // Check residual
        if (dpest < tol[i]) {
          conv = true;
          break;
        }

        etaold[i] = eta[i];
        psiold[i] = psi[i];
      }
      PHIST_SOUT(PHIST_VERBOSE,"\n");
      if (conv) break;
    }
    if (conv) break;

    // rho = r'*rp
    _ST_ rho[numSys];
    PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(r, rp, rho, iflag), *iflag);

    // b = rho / rhoold
    _ST_ b[numSys];
    for (int i=0; i<numSys; i++) b[i] = rho[i] / rhoold[i];

    // u = r + b q
    PHIST_CHK_IERR(SUBR(mvec_set_block)(u, r, 0, numSys-1, iflag), *iflag);
    PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(b, q, st::one(), u, iflag), *iflag);

    // p = u + b(q + b p)
    PHIST_CHK_IERR(SUBR(mvec_set_block)(t, q, 0, numSys-1, iflag), *iflag);
    PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(b, p, st::one(), t, iflag), *iflag);
    PHIST_CHK_IERR(SUBR(mvec_set_block)(p, u, 0, numSys-1, iflag), *iflag);
    PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(b, t, st::one(), p, iflag), *iflag);

    // v = A K p
    if (Pop) {
      PHIST_CHK_IERR(Aop->apply(st::one(), Aop->A, p, st::zero(), t, iflag), *iflag);
      PHIST_CHK_IERR(Pop->apply(st::one(), Pop->A, t, st::zero(), v, iflag), *iflag);
    } else {
      PHIST_CHK_IERR(Aop->apply(st::one(), Aop->A, p, st::zero(), v, iflag), *iflag);
    }

    for (int i=0; i<numSys; i++) rhoold[i] = rho[i];
    for (int i=0; i<numSys; i++) dpold[i] = dp[i];
  }

  PHIST_CHK_IERR(SUBR(mvec_delete) (r,  iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_delete) (rp, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_delete) (v,  iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_delete) (t,  iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_delete) (q,  iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_delete) (p,  iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_delete) (u,  iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_delete) (d,  iflag), *iflag);
}
