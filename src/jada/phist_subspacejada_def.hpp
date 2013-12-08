//! Tries to compute a partial schur form $(Q,R)$ of dimension nEig
//! of the stencil $A*x-\lambda*B*x$ with a general linear operator $A$ and a
//! hermitian positive definite (hpd.) linear operator $B$ using a
//! block-Jacobi-Davidson QR method. <br>
//! The generalized eigenvalues $\lambda_i$ are the diagonal entries of the
//! partial schur form $A*Q = B*Q*R$ returned. <br>
//!
//! Input arguments:
//!
//! A_op:     pointer to the operator A
//! B_op:     pointer to the hpd. operator B (if B==NULL, B=I is assumed)
//! v0:       start vector to construct a start basis using <minBase> Arnoldi-iteraions
//! which:    decide at which end of the spectrum to look for eigenvalues
//! tol:      convergence tolerance
//! nEig:     number of desired eigenpairs (e.g. dimension of the partial schur form)
//! nIter:    maximum number of iterations allowed
//! blockDim: block size, calculates <blockDim> corrections in each iteration
//! minBase:  start up from a basis consisting of minBas vectors (using Arnoldi)
//! maxBase:  when the basis reaches <maxBase> vectors, restart from <minBase> vectors.
//! 
//! Output arguments:
//!
//! nEig:     number of converged eigenpairs (e.g. dimension of (Q,R))
//! Q:        orthogonal vectors of the partial schur form (Q,R)
//! R:        small upper triangular matrix of the partial schur form (Q,R)
//! nIter:    number of iterations performed
//! resNorm:  norm of the residua of the schur form $A*q_i-Q*r_i, i=1,nEig$
//! ierr:     return code of the solver (0 on success, negative on error, positive on warning)
//!
#ifndef SUBR  /* ignore the following, just to allow better syntax highlighting/autocompletion etc */
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <iomanip>

#include "phist_macros.h"
#include "phist_subspacejada.h"
#include "phist_kernels.h"
#include "phist_lapack.h"
#include "phist_orthog.h"

#include "phist_ScalarTraits.hpp"
#include "jada_helpers.hpp"
#include "phist_schur_decomp.h"
#include "phist_jadaOp.h"
#include "phist_simple_arnoldi.h"

#include "phist_bgmres.h"

#include "phist_gen_s.h"
#endif
void SUBR(subspacejada)( TYPE(const_op_ptr) A_op,  TYPE(const_op_ptr) B_op,
                         TYPE(const_mvec_ptr) v0,  eigSort_t which,
                         _MT_ tol,                 int nEig,
                         int* nIter,               int blockDim,
                         int minBase,              int maxBase,
                         TYPE(mvec_ptr) Q,         TYPE(sdMat_ptr) R,
                         _MT_* resNorm,            int* ierr)
{
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *ierr = 0;


  //------------------------------- check arguments --------------------------------
  if( nEig < blockDim )
  {
    PHIST_SOUT(PHIST_ERROR, "parameter nEig > minBase!");
    PHIST_CHK_IERR(*ierr = -99, *ierr);
  }
  if( minBase < nEig )
  {
    PHIST_SOUT(PHIST_ERROR, "parameter minBase < nEig!");
    PHIST_CHK_IERR(*ierr = -99, *ierr);
  }
  if( minBase > maxBase )
  {
    PHIST_SOUT(PHIST_ERROR, "parameter minBase > maxBase!");
    PHIST_CHK_IERR(*ierr = -99, *ierr);
  }
  if( maxBase < nEig+blockDim )
  {
    PHIST_SOUT(PHIST_ERROR, "paramater maxBase < nEig+blockDim!");
    PHIST_CHK_IERR(*ierr = -99, *ierr);
  }
  if( B_op != NULL )
  {
    PHIST_SOUT(PHIST_ERROR,"case B_op != NULL (e.g. B != I) not implemented yet!");
    PHIST_CHK_IERR(*ierr = -99, *ierr);
  }

  // set output format for floating point numbers
  std::cout << std::scientific << std::setprecision(4) << std::setw(15);


  //------------------------------- create vectors and matrices --------------------
  // get communicator for sdMats
  const_comm_ptr_t domain_comm;
  PHIST_CHK_IERR(phist_map_get_comm(A_op->domain_map, &domain_comm, ierr), *ierr);
  const_comm_ptr_t range_comm;
  PHIST_CHK_IERR(phist_map_get_comm(A_op->range_map,  &range_comm,  ierr), *ierr);

  // create mvecs and sdMats
  mvec_ptr_t  V_      = NULL;    //< space for V
  mvec_ptr_t  Vtmp_   = NULL;    //< temporary space for V
  mvec_ptr_t  AV_     = NULL;    //< space for AV
  mvec_ptr_t  BV_     = NULL;    //< space for BV
  mvec_ptr_t  BQ      = NULL;    //< B*Q
  mvec_ptr_t  t_      = NULL;    //< space for t
  mvec_ptr_t  res     = NULL;    //< residuum A*Q-Q*R

  sdMat_ptr_t H_      = NULL;    //< space for H
  sdMat_ptr_t Htmp_   = NULL;    //< temporary space for H
  sdMat_ptr_t Q_H_    = NULL;    //< space for Q_H
  sdMat_ptr_t R_H_    = NULL;    //< space for R_H
  _ST_ sigma[nEig];              //< JaDa correction shifts

  _ST_ *Q_H_raw       = NULL;
  _ST_ *R_H_raw       = NULL;
  _ST_ *Htmp_raw 			= NULL;
  lidx_t ldaQ_H, ldaR_H, ldaHtmp;

  PHIST_CHK_IERR(SUBR( mvec_create  ) (&V_,     A_op->domain_map, maxBase,        ierr), *ierr);
  // TODO: remove Vtmp
  PHIST_CHK_IERR(SUBR( mvec_create  ) (&Vtmp_,  A_op->domain_map, maxBase,                ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_create  ) (&AV_,    A_op->range_map,  maxBase,        	      ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_create  ) (&t_,     A_op->domain_map, nEig,                  ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_create  ) (&res,    A_op->range_map,  nEig,                  ierr), *ierr);

  PHIST_CHK_IERR(SUBR( sdMat_create ) (&H_,     maxBase,          maxBase,  range_comm,   ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_create ) (&Htmp_,  maxBase,          maxBase,  range_comm,   ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_create ) (&Q_H_,   maxBase,          maxBase,  range_comm,   ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_create ) (&R_H_,   maxBase,          maxBase,  range_comm,   ierr), *ierr);

  PHIST_CHK_IERR(SUBR( sdMat_extract_view ) (Q_H_,    &Q_H_raw,   &ldaQ_H,   ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_extract_view ) (R_H_,  	&R_H_raw,   &ldaR_H,   ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_extract_view ) (Htmp_,   &Htmp_raw,  &ldaHtmp,  ierr), *ierr);
  if( B_op != NULL )
  {
    PHIST_CHK_IERR(SUBR( mvec_create )(&BV_,    B_op->range_map,  maxBase,                ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_create )(&BQ,    B_op->range_map,  nEig,                  ierr), *ierr);
  }
  else
  {
    BV_ = V_;
    BQ  = Q;
  }
  // array for the (possibly complex) eigenvalues for SchurDecomp
  CT* ev_H = new CT[maxBase];

  // create views on mvecs and sdMats with current dimensions
  int nV  = minBase;          //< current subspace dimension

  mvec_ptr_t  V   = NULL;     //< B-orthogonal basis of the search space
  mvec_ptr_t  Vtmp= NULL;     //< temporary V
  mvec_ptr_t  Vv  = NULL;     //< next columns in V_
  mvec_ptr_t  AV  = NULL;     //< A*V
  mvec_ptr_t  AVv = NULL;     //< next columns in AV_
  mvec_ptr_t  BV  = NULL;     //< B*V
  mvec_ptr_t  BVv = NULL;     //< next columns in BV_
  mvec_ptr_t  t   = NULL;     //< Block-Jacobi-Davidson correction
  mvec_ptr_t  t_res = NULL;   //< part of the residual AQ-QR corresponding to current block t
  mvec_ptr_t  Qtil= NULL;     //< view of part of Q required for the JaDa correction equation
  mvec_ptr_t BQtil= NULL;     //< B*Qtil

  sdMat_ptr_t H   = NULL;     //< projection of A onto H, V'*AV
  sdMat_ptr_t Htmp= NULL;     //< temporary space for H
  sdMat_ptr_t HVv = NULL;     //< next rows in H_
  sdMat_ptr_t HvV = NULL;     //< next columns in H_
  sdMat_ptr_t Hvv = NULL;     //< next block on the diagonal of H_
  sdMat_ptr_t r   = NULL;     //< currently iterated block of R_
  //sdMat_ptr_t Rr  = NULL;     //< [R a; 0 r]
  sdMat_ptr_t Q_H = NULL;     //< schur vectors of H
  sdMat_ptr_t Qq_H = NULL;
  sdMat_ptr_t R_H = NULL;     //< schur matrix of H
  sdMat_ptr_t Rr_H = NULL;

  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (V_,      &V,                       0,     nV-1,      ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (AV_,     &AV,                      0,     nV-1,      ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BV_,     &BV,                      0,     nV-1,      ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,      &H,     0,      nV-1,     0,     nV-1,      ierr), *ierr);



  //------------------------------- initialize subspace etc ------------------------
  // run arnoldi
  nV = minBase;
  PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,      &H,     0,      minBase,  0,     nV-1,      ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (V_,      &V,                       0,     nV,        ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BV_,     &BV,                      0,     nV,        ierr), *ierr);

  //TODO B_op in arnoldi
  // calculates A*V(:,1:m) = V(:,1:m+1)*H(1:m+1,1:m)
  PHIST_CHK_IERR(SUBR( simple_arnoldi ) (A_op, B_op, v0, V, BV, H, nV, ierr), *ierr);

  // calculate AV from V,H
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (AV_,     &AV,                      0,     nV-1,      ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_times_sdMat ) (st::one(), V, H,  st::zero(), AV, ierr), *ierr);

  // calculate H and setup V, BV
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (V_,      &V,                       0,     nV-1,      ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BV_, 		&BV,                      0,     nV-1,      ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,  		&H,     0,      nV-1,     0,     nV-1,      ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvecT_times_mvec ) (st::one(), V, AV, st::zero(), H,  ierr), *ierr);




  //----------------------------------- MAIN LOOP ----------------------------------
  int maxIter = *nIter;
  for(*nIter = 0; *nIter < maxIter; (*nIter)++)
  {
    // update views
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Q_H_,&Q_H, 0,     nV-1,      0,     nV-1,      ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Q_H_,&Qq_H,0,     nV-1,      0,     nEig-1,    ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (R_H_,&R_H, 0,     nV-1,      0,     nV-1,      ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (R_H_,&Rr_H,0,     nEig-1,    0,     nEig-1,    ierr), *ierr);

#ifdef TESTING
{
  // check orthogonality of V, BV, Q
  PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Htmp_,&Htmp,0,    nV-1,      0,     nV-1,      ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvecT_times_mvec ) (st::one(), V, BV, st::zero(), Htmp, ierr), *ierr);
  _MT_ orthEps = std::abs(Htmp_raw[0] - st::one());
  for(int i = 0; i < nV; i++)
    for(int j = 0; j < nV; j++)
      orthEps = std::max(orthEps, std::abs(Htmp_raw[i*ldaHtmp+j] - ((i==j) ? st::one() : st::zero())));
  PHIST_OUT(PHIST_INFO, "B-orthogonality of V: %e", orthEps);

  // check AV = A*V
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (Vtmp_, &Vtmp,                  0,       nV-1,          ierr), *ierr);
  PHIST_CHK_IERR( A_op->apply(st::one(), A_op->A, V, st::zero(), Vtmp, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_add_mvec ) (-st::one(), AV, st::one(), Vtmp, ierr), *ierr);
  _MT_ normVtmp[nV];
  PHIST_CHK_IERR(SUBR( mvec_norm2 ) (Vtmp, normVtmp, ierr), *ierr);
  _MT_ equalEps = normVtmp[0];
  for(int i = 0; i < nV; i++)
    equalEps = std::max(equalEps, normVtmp[i]);
  PHIST_OUT(PHIST_INFO, "AV - A*V: %e", equalEps);
  // check H = V'*A*V
  PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Htmp_,&Htmp,0,    nV-1,      0,     nV-1,      ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvecT_times_mvec ) (st::one(), V, AV, st::zero(), Htmp, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_add_sdMat ) (-st::one(), H, st::one(), Htmp, ierr), *ierr);
  equalEps = std::abs(Htmp_raw[0]);
  for(int i = 0; i < nV; i++)
    for(int j = 0; j < nV; j++)
      equalEps = std::max(equalEps, std::abs(Htmp_raw[i*ldaHtmp+j]));
  PHIST_OUT(PHIST_INFO, "H - V'*AV: %e", equalEps);
}
#endif


    // calculate sorted Schur form of H in (Q_H,R_H)
    PHIST_CHK_IERR(SUBR( sdMat_add_sdMat ) (st::one(), H, st::zero(), R_H, ierr), *ierr);
    int nSort = nEig;
    int nSelect = nSort;
    PHIST_CHK_IERR(SUBR( SchurDecomp ) (R_H_raw, ldaR_H, Q_H_raw, ldaQ_H, nV, nSelect, nSort, which, ev_H, ierr), *ierr);

    // calculate approximate Schur form of A (with deflation of converged vectors Q)
    PHIST_CHK_IERR(SUBR( mvec_times_sdMat ) (st::one(), V,    Qq_H, st::zero(), Q,   ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_add_sdMat  ) (st::one(), Rr_H,       st::zero(), R,   ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_times_sdMat ) (st::one(), AV,   Qq_H, st::zero(), res, ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_times_sdMat ) (st::one(), BV,   Qq_H, st::zero(), BQ,  ierr), *ierr);
    // overwrite res with the residuum: -res = -(Aq - Bqr) = + BQq*r - res
    PHIST_CHK_IERR(SUBR( mvec_times_sdMat ) (st::one(), BQ,   R,    -st::one(), res, ierr), *ierr);
    // calculate norm of the residuum
    PHIST_CHK_IERR(SUBR( mvec_norm2 ) (res, resNorm, ierr), *ierr);

    // check for converged eigenvalues
    int nConvergedEig = 0;
    for(int i = 0; i < nEig; i++)
    {
      PHIST_SOUT(PHIST_INFO,"In iteration %d: Current approximation for eigenvalue %d is %16.8g%+16.8gi with residuum %e", *nIter, i+1, ct::real(ev_H[i]),ct::imag(ev_H[i]), resNorm[i]);
      if( resNorm[i] <= tol )
        nConvergedEig++;
    }
    PHIST_SOUT(PHIST_INFO,"In iteration %d: converged eigenvalues %d", *nIter, nConvergedEig);
    if( nConvergedEig == nEig )
      break;


    // setup matrix of shifts and residuals for the correction equation
    int k = 0;  // current correction block dimension
    int k_ = 0; // 0:k_-1 vectors of Q used for the orthogonal projection in the correction equation
    for(int i = 0; i < nEig && k < blockDim; i++)
    {
#ifndef IS_COMPLEX
      if( std::abs(ct::imag(ev_H[i])) > tol )
      {
        PHIST_SOUT(PHIST_WARNING,"real case with complex conjugate eigenvalues not fully implemented yet!");
      }
#endif

      if( resNorm[i] > tol )
      {
        k_ = i;
#ifndef IS_COMPLEX
        sigma[k] = -ct::real(ev_H[i]);
#else
        sigma[k] = -ev_H[i];
#endif

        // copy residual
        if( i != k )
        {
          PHIST_CHK_IERR(SUBR( mvec_view_block ) (res, &t_res, k, k, ierr), *ierr);
          PHIST_CHK_IERR(SUBR( mvec_get_block  ) (res, t_res,  i, i, ierr), *ierr);
        }
        k++;
      }
    }


    // shrink search space if necessary
    if( nV + k > maxBase )
    {
      PHIST_SOUT(PHIST_INFO,"Shrinking search space from %d to %d", nV, minBase);

      //TODO: this should be done in-place!!!!!!!!!!!!!!
      PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Q_H_,  &Q_H,  0, nV-1,          0, minBase-1,    ierr), *ierr);
      //PHIST_CHK_IERR(SUBR( mvec_view_block  ) (Vtmp_, &Vtmp,                   0, nV-1,         ierr), *ierr);

      PHIST_CHK_IERR(SUBR( mvec_times_sdMat_inplace ) (V, Q_H, 64, ierr), *ierr);
      //PHIST_CHK_IERR(SUBR( mvec_add_mvec    ) (st::one(), V,           st::zero(), Vtmp, ierr), *ierr);
      PHIST_CHK_IERR(SUBR( mvec_view_block  ) (V_,    &V,                      0, minBase-1,    ierr), *ierr);
      //PHIST_CHK_IERR(SUBR( mvec_times_sdMat ) (st::one(), Vtmp, Q_H,   st::zero(), V,    ierr), *ierr);

      PHIST_CHK_IERR(SUBR( mvec_times_sdMat_inplace ) (AV, Q_H, 64, ierr), *ierr);
      //PHIST_CHK_IERR(SUBR( mvec_add_mvec    ) (st::one(), AV,          st::zero(), Vtmp, ierr), *ierr);
      PHIST_CHK_IERR(SUBR( mvec_view_block  ) (AV_,   &AV,                     0, minBase-1,    ierr), *ierr);
      //PHIST_CHK_IERR(SUBR( mvec_times_sdMat ) (st::one(), Vtmp, Q_H,   st::zero(), AV,   ierr), *ierr);

      if( B_op != NULL )
      {
      PHIST_CHK_IERR(SUBR( mvec_times_sdMat_inplace ) (BV, Q_H, 64, ierr), *ierr);
        //PHIST_CHK_IERR(SUBR( mvec_add_mvec    ) (st::one(), BV,        st::zero(), Vtmp, ierr), *ierr);
        PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BV_, &BV,                     0, minBase-1,    ierr), *ierr);
        //PHIST_CHK_IERR(SUBR( mvec_times_sdMat ) (st::one(), Vtmp, Q_H, st::zero(), BV,   ierr), *ierr);
      }
      else
      {
        PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BV_, &BV,                     0, minBase-1, ierr), *ierr);
      }

      // update H <- H_Q' * H * H_Q
      PHIST_CHK_IERR(SUBR( sdMat_view_block  )(Htmp_, &Htmp, 0, nV-1,         0, minBase-1,     ierr), *ierr);
      PHIST_CHK_IERR(SUBR( sdMat_times_sdMat )(st::one(), H,    Q_H,   st::zero(), Htmp, ierr), *ierr);
      PHIST_CHK_IERR(SUBR( sdMat_view_block  )(H_,    &H,    0, minBase-1,    0, minBase-1,     ierr), *ierr);
      PHIST_CHK_IERR(SUBR( sdMatT_times_sdMat)(st::one(), Q_H, Htmp,   st::zero(), H,    ierr), *ierr);

      nV = minBase;
    }


    // calculate corrections
    // setup jadaOp
    // set correction views and temporary jadaOp-storage
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (t_,  &t,     0, k-1,  ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (res, &t_res, 0, k-1,  ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (V_,  &Vv,                    nV,    nV+k-1,    ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (AV_, &AVv,                   nV,    nV+k-1,    ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BV_, &BVv,                   nV,    nV+k-1,    ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (R_H_,&R_H, 0,     k_-1,      0,     k-1,       ierr), *ierr);
    // we only need to view first part of Q
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (Q,   &Qtil,  0, k_-1, ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BQ,  &BQtil, 0, k_-1, ierr), *ierr);

    TYPE(op) jdOp;
    PHIST_CHK_IERR(SUBR( jadaOp_create ) (A_op, B_op, Qtil, BQtil, sigma, R_H, AVv, BVv, Vv, NULL, &jdOp, ierr), *ierr);
    // TODO specify useful bgmresIter and tol per eigenvalue!
    int bgmresIter = 10;
    PHIST_CHK_IERR(SUBR( mvec_put_value )(t, st::zero(), ierr), *ierr);
    PHIST_CHK_NEG_IERR(SUBR( bgmres )    (&jdOp, t, t_res, mt::zero(), &bgmresIter, 10, 1, NULL, ierr), *ierr);
    PHIST_CHK_IERR(SUBR( jadaOp_delete ) (&jdOp, ierr), *ierr);

    // enlarge search space
    // first update views
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,  &HVv, 0,     nV-1,      nV,    nV+k-1,    ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,  &HvV, nV,    nV+k-1,    0,     nV-1,      ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,  &Hvv, nV,    nV+k-1,    nV,    nV+k-1,    ierr), *ierr);
    // orthogonalize t as Vv (reuse R_H)
    PHIST_CHK_IERR(SUBR( mvec_add_mvec ) (st::one(), t, st::zero(), Vv, ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (R_H_,&R_H, 0,     nV-1,      0,     k-1,       ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (R_H_,&Rr_H,nV,    nV+k-1,    nV,    nV+k-1,    ierr), *ierr);
    PHIST_CHK_IERR(SUBR( orthog ) (V, Vv, Rr_H, R_H, 3, ierr), *ierr);
    // TODO: only take non-random vector if *ierr > 0
    // calculate AVv, BVv
    PHIST_CHK_IERR( A_op->apply(st::one(), A_op->A, Vv, st::zero(), AVv, ierr), *ierr);
    if( B_op != NULL )
    {
      PHIST_CHK_IERR( B_op->apply(st::one(), B_op->A, Vv, st::zero(), BVv, ierr), *ierr);
    }
    // update H
    PHIST_CHK_IERR(SUBR( mvecT_times_mvec ) (st::one(), V,  AVv, st::zero(), HVv, ierr), *ierr);
    //TODO: for the symmetric case use AVv*V here, so we don't need AV at all
    PHIST_CHK_IERR(SUBR( mvecT_times_mvec ) (st::one(), Vv, AV,  st::zero(), HvV, ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvecT_times_mvec ) (st::one(), Vv, AVv, st::zero(), Hvv, ierr), *ierr);
    // increase nV
    nV = nV + k;
    // update views
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (V_,  &V,                     0,     nV-1,      ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (AV_, &AV,                    0,     nV-1,      ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BV_, &BV,                    0,     nV-1,      ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,  &H,   0,     nV-1,      0,     nV-1,      ierr), *ierr);
  }



  //------------------------------- delete vectors and matrices --------------------
  // delete views
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (R_H, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Rr_H,ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Q_H, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Qq_H,ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Hvv, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (HvV, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (HVv, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (H,   ierr), *ierr);
  //PHIST_CHK_IERR(SUBR( sdMat_delete ) (Htmp,ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (r,   ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Htmp,ierr),*ierr);

  PHIST_CHK_IERR(SUBR( mvec_delete  ) (t,   ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (res, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (t_res,ierr),*ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (BVv, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (BV,  ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (AVv, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (AV,  ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (Vv,  ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (Vtmp,ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (V,   ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (Qtil,ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (BQtil,ierr), *ierr);

  if( B_op != NULL )
  {
    PHIST_CHK_IERR(SUBR( mvec_delete  ) (BQ,  ierr), *ierr);
  }

  delete[] ev_H;

  // delete mvecs and sdMats
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Q_H_,ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (R_H_,ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Htmp_,ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (H_,  ierr), *ierr);
  if( B_op != NULL )
  {
    PHIST_CHK_IERR(SUBR( mvec_delete )(BV_, ierr), *ierr);
  }
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (t_,  ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (AV_, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (Vtmp_,ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (V_,  ierr), *ierr);
}

