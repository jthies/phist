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

#include "phist_jadaInnerGmres.h"

#include "phist_gen_s.h"
#endif
void SUBR(subspacejada)( TYPE(const_op_ptr) A_op,  TYPE(const_op_ptr) B_op,
                         TYPE(const_mvec_ptr) v0,  eigSort_t which,
                         _MT_ tol,                 int nEig,
                         int* nIter,               int blockDim,
                         int minBase,              int maxBase,
                         int initialShiftIter,     _ST_ initialShift,
                         TYPE(mvec_ptr) Q,         TYPE(sdMat_ptr) R,
                         _MT_* resNorm,            int* ierr)
{
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *ierr = 0;

  // extended number of searched eigenvalues in order to respect the block dimension,
  // this way we always have a fixed blockDim AND it should make the calculation
  // of the last eigenvalues more stable in some cases
  int nEig_ = nEig + blockDim - 1;

  //------------------------------- check arguments --------------------------------
  if( minBase < nEig_ )
  {
    PHIST_SOUT(PHIST_ERROR, "parameter minBase < nEig+blockDim-1!\n");
    PHIST_CHK_IERR(*ierr = -99, *ierr);
  }
  if( minBase+blockDim > maxBase )
  {
    PHIST_SOUT(PHIST_ERROR, "parameter minBase+blockDim > maxBase!\n");
    PHIST_CHK_IERR(*ierr = -99, *ierr);
  }
  if( maxBase < nEig+blockDim )
  {
    PHIST_SOUT(PHIST_ERROR, "paramater maxBase < nEig+blockDim!\n");
    PHIST_CHK_IERR(*ierr = -99, *ierr);
  }
  if( B_op != NULL )
  {
    PHIST_SOUT(PHIST_ERROR,"case B_op != NULL (e.g. B != I) not implemented yet!\n");
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
  _ST_ sigma[nEig_];             //< JaDa correction shifts

  _ST_ *Q_H_raw       = NULL;
  _ST_ *R_H_raw       = NULL;
  _ST_ *Htmp_raw 			= NULL;
  lidx_t ldaQ_H, ldaR_H, ldaHtmp;

  PHIST_CHK_IERR(SUBR( mvec_create  ) (&V_,     A_op->domain_map, maxBase,        ierr), *ierr);
  // TODO: remove Vtmp
  PHIST_CHK_IERR(SUBR( mvec_create  ) (&Vtmp_,  A_op->domain_map, maxBase,                ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_create  ) (&AV_,    A_op->range_map,  maxBase,        	      ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_create  ) (&t_,     A_op->domain_map, nEig_,                 ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_create  ) (&res,    A_op->range_map,  nEig_,                 ierr), *ierr);

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
    PHIST_CHK_IERR(SUBR( mvec_create )(&BQ,    B_op->range_map,  nEig_,                 ierr), *ierr);
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
  sdMat_ptr_t Hh  = NULL;     //< inside view for H
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


  //------------------------------- initialize GMRES solver ------------------------
  TYPE(jadaInnerGmresState_ptr) *gmresState = new TYPE(jadaInnerGmresState_ptr)[blockDim];
  _MT_ *gmresResNorm = new _MT_[blockDim];
  int nTotalGmresIter = 0;
  PHIST_CHK_IERR(SUBR( jadaInnerGmresStates_create ) (gmresState, blockDim, A_op->domain_map, 11, ierr), *ierr);

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
  int nConvergedEig = 0;
  for(*nIter = 0; *nIter < maxIter; (*nIter)++)
  {
#ifdef TESTING
{
  // check orthogonality of V, BV, Q
  PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Htmp_,&Htmp,0,    nV-1,      0,     nV-1,      ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvecT_times_mvec ) (st::one(), V, BV, st::zero(), Htmp, ierr), *ierr);
  _MT_ orthEps = std::abs(Htmp_raw[0] - st::one());
  for(int i = 0; i < nV; i++)
    for(int j = 0; j < nV; j++)
      orthEps = std::max(orthEps, std::abs(Htmp_raw[i*ldaHtmp+j] - ((i==j) ? st::one() : st::zero())));
  PHIST_OUT(PHIST_INFO, "B-orthogonality of V: %e\n", orthEps);

  // check AV = A*V
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (Vtmp_, &Vtmp,                  0,       nV-1,          ierr), *ierr);
  PHIST_CHK_IERR( A_op->apply(st::one(), A_op->A, V, st::zero(), Vtmp, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_add_mvec ) (-st::one(), AV, st::one(), Vtmp, ierr), *ierr);
  _MT_ normVtmp[nV];
  PHIST_CHK_IERR(SUBR( mvec_norm2 ) (Vtmp, normVtmp, ierr), *ierr);
  _MT_ equalEps = normVtmp[0];
  for(int i = 0; i < nV; i++)
    equalEps = std::max(equalEps, normVtmp[i]);
  PHIST_OUT(PHIST_INFO, "AV - A*V: %e\n", equalEps);
  // check H = V'*A*V
  PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Htmp_,&Htmp,0,    nV-1,      0,     nV-1,      ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvecT_times_mvec ) (st::one(), V, AV, st::zero(), Htmp, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_add_sdMat ) (-st::one(), H, st::one(), Htmp, ierr), *ierr);
  equalEps = std::abs(Htmp_raw[0]);
  for(int i = 0; i < nV; i++)
    for(int j = 0; j < nV; j++)
      equalEps = std::max(equalEps, std::abs(Htmp_raw[i*ldaHtmp+j]));
  PHIST_OUT(PHIST_INFO, "H - V'*AV: %e\n", equalEps);
}
#endif


    // calculate sorted Schur form of H in (Q_H,R_H)
    // we only update part of Q_H,R_H, so first set Q_H, R_H to zero
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Q_H_,&Q_H, 0,     nV-1,      0,     nV-1,      ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (R_H_,&R_H, 0,     nV-1,      0,     nV-1,      ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_put_value  ) (Q_H, st::zero(), ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_put_value  ) (R_H, st::zero(), ierr), *ierr);
    // then copy the new block of H
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (R_H_,&R_H, nConvergedEig, nV-1, nConvergedEig, nV-1, ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_get_block  ) (H,   R_H,  nConvergedEig, nV-1, nConvergedEig, nV-1, ierr), *ierr);
    int nSort = nEig_-nConvergedEig;
    int nSelect = nSort;
    lidx_t offR_H = ldaR_H*nConvergedEig+nConvergedEig;
    lidx_t offQ_H = ldaQ_H*nConvergedEig+nConvergedEig;
    PHIST_CHK_IERR(SUBR( SchurDecomp ) (R_H_raw+offR_H, ldaR_H, Q_H_raw+offQ_H, ldaQ_H, nV-nConvergedEig, nSelect, nSort, which, ev_H+nConvergedEig, ierr), *ierr);
    // we still need to add the missing parts of R_H, Q_H
    if( nConvergedEig > 0 )
    {
      // upper left part of Q_H
      for(int i = 0; i < nConvergedEig; i++)
        Q_H_raw[ldaQ_H*i+i] = st::one();

      // upper left part of R_H
      PHIST_CHK_IERR(SUBR( sdMat_view_block ) (R_H_, &R_H, 0, nConvergedEig-1, 0, nConvergedEig-1, ierr), *ierr);
      PHIST_CHK_IERR(SUBR( sdMat_get_block  ) (R,    R_H,  0, nConvergedEig-1, 0, nConvergedEig-1, ierr), *ierr);

      // upper right part of R_H
      PHIST_CHK_IERR(SUBR( sdMat_view_block  ) (Q_H_, &Qq_H, nConvergedEig, nEig_-1,         nConvergedEig, nEig_-1, ierr), *ierr);
      PHIST_CHK_IERR(SUBR( sdMat_view_block  ) (R_H_, &Rr_H, 0,             nConvergedEig-1, nConvergedEig, nEig_-1, ierr), *ierr);
      PHIST_CHK_IERR(SUBR( sdMat_view_block  ) (H_  , &Hh,   0,             nConvergedEig-1, nConvergedEig, nEig_-1, ierr), *ierr);

      PHIST_CHK_IERR(SUBR( sdMat_times_sdMat ) (st::one(), Hh, Qq_H, st::zero(), Rr_H, ierr), *ierr);
    }

    // update views
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Q_H_,&Q_H, 0,     nV-1,      0,     nV-1,      ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Q_H_,&Qq_H,0,     nV-1,      0,     nEig_-1,   ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (R_H_,&R_H, 0,     nV-1,      0,     nV-1,      ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (R_H_,&Rr_H,0,     nEig_-1,   0,     nEig_-1,   ierr), *ierr);



    // calculate approximate Schur form of A
    PHIST_CHK_IERR(SUBR( mvec_times_sdMat ) (st::one(), V,    Qq_H, st::zero(), Q,   ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_add_sdMat  ) (st::one(), Rr_H,       st::zero(), R,   ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_times_sdMat ) (st::one(), AV,   Qq_H, st::zero(), res, ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_times_sdMat ) (st::one(), BV,   Qq_H, st::zero(), BQ,  ierr), *ierr);
    // overwrite res with the residuum: -res = -(Aq - Bqr) = + BQq*r - res
    PHIST_CHK_IERR(SUBR( mvec_times_sdMat ) (st::one(), BQ,   R,    -st::one(), res, ierr), *ierr);
    // calculate norm of the residuum
    PHIST_CHK_IERR(SUBR( mvec_norm2 ) (res, resNorm, ierr), *ierr);



    // check for converged eigenvalues
    int nNewlyConvergedEig = 0;
    for(int i = nConvergedEig; i < nEig; i++)
    {
      PHIST_SOUT(PHIST_INFO,"In iteration %d: Current approximation for eigenvalue %d is %16.8g%+16.8gi with residuum %e\n", *nIter, i+1, ct::real(ev_H[i]),ct::imag(ev_H[i]), resNorm[i]);
      if( resNorm[i] <= tol && i == nConvergedEig+nNewlyConvergedEig ) // only consider eigenvalues from the beginning
        nNewlyConvergedEig++;
    }

    if( nNewlyConvergedEig > 0 )
    {
      PHIST_SOUT(PHIST_INFO,"In iteration %d: locking %d newly converged eigenvalues\n", *nIter, nNewlyConvergedEig);

      // reorder V and H
      PHIST_CHK_IERR(SUBR( transform_searchSpace ) (V, AV, BV, H, Q_H, B_op != NULL, ierr), *ierr);

      nConvergedEig = nConvergedEig+nNewlyConvergedEig;
    }

    if( nConvergedEig >= nEig )
    {
      PHIST_SOUT(PHIST_INFO,"In iteration %d: all eigenvalues converged!\n", *nIter);
      break;
    }

    if( *nIter >= maxIter )
    {
      PHIST_SOUT(PHIST_INFO,"Reached maximum number of iterations!\n");
      break;
    }



    // setup matrix of shifts and residuals for the correction equation
    int k_ = 0; // 0:k_ vectors of Q used for the orthogonal projection in the correction equation
    int k = 0;  // is always == blockDim!
    for(int i = 0; i < nEig_ && k < blockDim; i++)
    {
#ifndef IS_COMPLEX
      if( std::abs(ct::imag(ev_H[i])) > tol )
      {
        PHIST_SOUT(PHIST_WARNING,"real case with complex conjugate eigenvalues not fully implemented yet!\n");
      }
#endif

      if( resNorm[i] > tol )
      {
        k_ = i;

        if( *nIter < initialShiftIter )
        {
          sigma[k] = -initialShift;
        }
        else
        {
#ifndef IS_COMPLEX
          sigma[k] = -ct::real(ev_H[i]);
#else
          sigma[k] = -ev_H[i];
#endif
        }

        // copy residual
        if( i != k )
        {
          PHIST_CHK_IERR(SUBR( mvec_view_block ) (res, &t_res, k, k, ierr), *ierr);
          PHIST_CHK_IERR(SUBR( mvec_get_block  ) (res, t_res,  i, i, ierr), *ierr);
        }
        k++;
      }
    }
    // we should also project out other schur vectors if they have already (nearly) converged
    for(int i = k_+1; i < nEig_; i++)
    {
      if( resNorm[i] <= sqrt(tol) )
      {
        if( i != k_ )
        {
          PHIST_CHK_IERR(SUBR( mvec_view_block ) (Q, &t_res, k_, k_, ierr), *ierr);
          PHIST_CHK_IERR(SUBR( mvec_get_block  ) (Q, t_res,  i,  i,  ierr), *ierr);
        }
        k_++;
      }
    }


    // shrink search space if necessary
    if( nV + k > maxBase )
    {
      PHIST_SOUT(PHIST_INFO,"Shrinking search space from %d to %d\n", nV, minBase);

      // nothing to do if no converged eigenvalues this iteration
      if( nNewlyConvergedEig == 0 )
      {
        PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Q_H_,  &Q_H,  0, nV-1,          0, minBase-1,    ierr), *ierr);

        PHIST_CHK_IERR(SUBR( transform_searchSpace ) (V, AV, BV, H, Q_H, B_op != NULL, ierr), *ierr);
      }

      // update views
      PHIST_CHK_IERR(SUBR( mvec_view_block  ) (V_,    &V,                      0, minBase-1,    ierr), *ierr);
      PHIST_CHK_IERR(SUBR( mvec_view_block  ) (AV_,   &AV,                     0, minBase-1,    ierr), *ierr);
      PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BV_,   &BV,                     0, minBase-1,    ierr), *ierr);
      PHIST_CHK_IERR(SUBR( sdMat_view_block  )(H_,    &H,    0, minBase-1,     0, minBase-1,    ierr), *ierr);

      nV = minBase;
    }
#ifdef PHIST_KERNEL_LIB_FORTRAN
int originalBlockDim = k;
if( k < blockDim )
{
  PHIST_SOUT(PHIST_WARNING,"k < blockDim not supported for fortran kernel lib!\n");
  k = blockDim;
}
#endif

    // calculate corrections
    // setup jadaOp
    // set correction views and temporary jadaOp-storage
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (t_,  &t,     0, k-1,  ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (res, &t_res, 0, k-1,  ierr), *ierr);
    // we only need to view first part of Q
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (Q,   &Qtil,  0, k_, ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BQ,  &BQtil, 0, k_, ierr), *ierr);

    TYPE(op) jdOp;
    PHIST_CHK_IERR(SUBR( jadaOp_create ) (A_op, B_op, Qtil, BQtil, sigma, k, &jdOp, ierr), *ierr);
    // TODO specify useful tol per eigenvalue!
    PHIST_CHK_IERR(SUBR( mvec_put_value )(t, st::zero(), ierr), *ierr);
    for(int i = 0; i < k; i++)
    {
      PHIST_CHK_IERR(SUBR( mvec_view_block  ) (t_,&t, i,i, ierr), *ierr);
      PHIST_CHK_IERR(SUBR( mvec_view_block  ) (res,&t_res, i,i, ierr), *ierr);
      PHIST_CHK_IERR(SUBR( jadaInnerGmresState_reset ) (gmresState[i], t_res, t, ierr), *ierr);
      gmresState[i]->tol = mt::zero();
    }
    PHIST_CHK_NEG_IERR(SUBR( jadaInnerGmresStates_iterate ) (&jdOp, gmresState, k, &nTotalGmresIter, ierr), *ierr);
#ifdef PHIST_KERNEL_LIB_FORTRAN
  k = originalBlockDim;
#endif
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (t_,&t, 0,k-1, ierr), *ierr);
    PHIST_CHK_IERR(SUBR( jadaInnerGmresStates_updateSol ) (gmresState, k, t, NULL, gmresResNorm, ierr), *ierr);
    PHIST_CHK_IERR(SUBR( jadaOp_delete ) (&jdOp, ierr), *ierr);

    // enlarge search space
    // first update views
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,  &HVv, 0,     nV-1,      nV,    nV+k-1,    ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,  &HvV, nV,    nV+k-1,    0,     nV-1,      ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,  &Hvv, nV,    nV+k-1,    nV,    nV+k-1,    ierr), *ierr);
    // orthogonalize t as Vv (reuse R_H)
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (V_,  &Vv,                    nV,    nV+k-1,    ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_add_mvec ) (st::one(), t, st::zero(), Vv, ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (R_H_,&R_H, 0,     nV-1,      0,     k-1,       ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (R_H_,&Rr_H,nV,    nV+k-1,    nV,    nV+k-1,    ierr), *ierr);
    PHIST_CHK_NEG_IERR(SUBR( orthog ) (V, Vv, Rr_H, R_H, 3, ierr), *ierr);
    // TODO: only take non-random vector if *ierr > 0
    // calculate AVv, BVv
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (AV_, &AVv,                   nV,    nV+k-1,    ierr), *ierr);
    PHIST_CHK_IERR( A_op->apply(st::one(), A_op->A, Vv, st::zero(), AVv, ierr), *ierr);
    if( B_op != NULL )
    {
      PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BV_, &BVv,                   nV,    nV+k-1,    ierr), *ierr);
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
  PHIST_CHK_IERR(SUBR( jadaInnerGmresStates_delete ) (gmresState, blockDim, ierr), *ierr);
  delete[] gmresState;
  delete[] gmresResNorm;

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

