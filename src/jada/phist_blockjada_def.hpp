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
//! Q_:       orthogonal vectors of the partial schur form (Q,R)
//! R_:       small upper triangular matrix of the partial schur form (Q,R)
//! nIter:    number of iterations performed
//! resNorm:  norm of the residua of the schur form $A*q_i-Q*r_i, i=1,nEig$
//! ierr:     return code of the solver (0 on success, negative on error, positive on warning)
//!
void SUBR(blockjada)( TYPE(const_op_ptr) A_op,  TYPE(const_op_ptr) B_op,
                      TYPE(const_mvec_ptr) v0,  eigSort_t which,
                      _MT_ tol,                 int* nEig,
                      int* nIter,               int blockDim,
                      int minBase,              int maxBase,
                      TYPE(mvec_ptr) Q_,        TYPE(sdMat_ptr) R_,
                      _MT_* resNorm,            int* ierr)
{
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *ierr = 0;


  //------------------------------- check arguments --------------------------------
  if( B_op != NULL )
  {
    PHIST_SOUT(PHIST_ERROR,"case B_op != NULL (e.g. B != I) not implemented yet!");
    PHIST_CHK_IERR(*ierr = 99, *ierr);
  }



  //------------------------------- create vectors and matrices --------------------
  // create mvecs and sdMats
  mvec_ptr_t  V_   = NULL;    //< space for V
  mvec_ptr_t  AV_  = NULL;    //< space for AV
  mvec_ptr_t  BV_  = NULL;    //< space for BV
  mvec_ptr_t  BQ_  = NULL;    //< space for BQ
  mvec_ptr_t  t_   = NULL;    //< space for t_
  sdMat_ptr_t H_   = NULL;    //< space for H
  sdMat_ptr_t Q_H_ = NULL;    //< space for Q_H
  sdMat_ptr_t R_H_ = NULL;    //< space for R_H
  _ST_ *Q_H_raw = NULL;
  _ST_ *R_H_raw = NULL;
  lidx_t ldaQ_H, ldaR_H;
  PHIST_CHK_IERR(SUBR( mvec_create  ) (&V_,   A_op->domain_map, maxBase,  ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_create  ) (&AV_,  A_op->range_map,  maxBase,  ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_create  ) (&t_,   A_op->domain_map, blockDim, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_create ) (&H_,   maxBase, maxBase, NULL,     ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_create ) (&Q_H_, maxBase, maxBase, NULL,     ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_create ) (&R_H_, maxBase, maxBase, NULL,     ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_extract_view ) (Q_H_, &Q_H_raw, &ldaQ_H, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_extract_view ) (R_H_, &R_H_raw, &ldaR_H, ierr), *ierr);
  if( B_op != NULL )
  {
    PHIST_CHK_IERR(SUBR( mvec_create )(&BV_, B_op->range_map,  maxBase,  ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_create )(&BQ_, B_op->range_map,  *nEig,    ierr), *ierr);
  }
  else
  {
    BV_ = V_;
    BQ_ = Q_;
  }
  // array for the (possibly complex) eigenvalues for SchurDecomp
  CT* ev_H = new CT[maxBase];


  // create views on mvecs and sdMats with current dimensions
  int nV  = minBase;          //< current subspace dimension
  int k   = blockDim;         //< current block size (can be if only few eigenvalues remain)
  int maxEig = *nEig;         //< number of eigenvalues to compute
  *nEig = 0;                  //< already converged number of eigenvalues

  mvec_ptr_t  V   = NULL;     //< B-orthogonal basis of the search space
  mvec_ptr_t  Vv  = NULL;     //< next columns in V_
  mvec_ptr_t  AV  = NULL;     //< A*V
  mvec_ptr_t  AVv = NULL;     //< next columns in AV_
  mvec_ptr_t  BV  = NULL;     //< B*V
  mvec_ptr_t  BVv = NULL;     //< next columns in BV_
  mvec_ptr_t  t   = NULL;     //< Block-Jacobi-Davidson correction
  mvec_ptr_t  Q   = NULL;     //< already converged schur vectors
  mvec_ptr_t  Qq  = NULL;     //< currently iterated block of Q_
  mvec_ptr_t  QQ  = NULL;     //< [Q Qq]
  mvec_ptr_t  BQ  = NULL;     //< B*Q
  mvec_ptr_t  BQq = NULL;     //< B*Qq
  mvec_ptr_t  BQQ = NULL;     //< B*QQ

  sdMat_ptr_t H   = NULL;     //< projection of A onto H, V'*AV
  sdMat_ptr_t HVv = NULL;     //< next rows in H_
  sdMat_ptr_t HvV = NULL;     //< next columns in H_
  sdMat_ptr_t Hvv = NULL;     //< next block on the diagonal of H_
  sdMat_ptr_t r   = NULL;     //< currently iterated block of R_
  sdMat_ptr_t R   = NULL;     //< already converged schur matrix
  sdMat_ptr_t a   = NULL;     //< part of R which is blended out by deflation (e.g. BQ'*q, before q<-q-B*a)
  sdMat_ptr_t Rr  = NULL;     //< [R a; 0 r]
  sdMat_ptr_t Q_H = NULL;     //< schur vectors of H
  sdMat_ptr_t R_H = NULL;     //< schur matrix of H

  // set R_ to zero because we don't explicitly consider the lower left part
  PHIST_CHK_IERR(SUBR( sdMat_put_value  ) (R_, st::zero(), ierr), *ierr);




  //------------------------------- initialize subspace etc ------------------------
  // run arnoldi
  nV = minBase;
  PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,  &H,  0, minBase,        0,     nV-1,      ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (V_,  &V,                     0,     nV,        ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BV_, &BV,                    0,     nV,        ierr), *ierr);

  //TODO B_op in arnoldi
  // calculates A*V(:,1:m) = V(:,1:m+1)*H(1:m+1,1:m)
  PHIST_CHK_IERR(SUBR( simple_arnoldi ) (A_op, B_op, v0, V, BV, H, nV, ierr), *ierr);

  // calculate AV from V,H
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (AV_, &AV,                    0,     nV-1,      ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_times_sdMat ) (st::one(), V, H,  st::zero(), AV, ierr), *ierr);

  // calculate H and setup V, BV
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (V_, &V,                      0,     nV-1,      ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (V_, &BV,                     0,     nV-1,      ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,  &H,   0,     nV-1,      0,     nV-1,      ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvecT_times_mvec ) (st::one(), V, AV, st::zero(), H,  ierr), *ierr);




  //----------------------------------- MAIN LOOP ----------------------------------
  int maxIter = *nIter;
  for(*nIter = 0; *nIter < maxIter; *nIter++)
  {
    // update views for the current iteration
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (V_,  &V,                     0,     nV-1,      ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (V_,  &Vv,                    nV,    nV+k-1,    ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (AV_, &AV,                    0,     nV-1,      ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (AV_, &AVv,                   nV,    nV+k-1,    ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BV_, &BV,                    0,     nV-1,      ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BV_, &BVv,                   nV,    nV+k-1,    ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (t_,  &t,                     0,     k-1,       ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (Q_,  &Qq,                    *nEig, *nEig+k-1, ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (Q_,  &QQ,                    0,     *nEig+k-1, ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BQ_, &BQq,                   *nEig, *nEig+k-1, ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BQ_, &BQQ,                   0,     *nEig+k-1, ierr), *ierr);

    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (R_,  &r,   *nEig, *nEig+k-1, *nEig, *nEig+k-1, ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,  &H,   0,     nV-1,      0,     nV-1,      ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,  &HVv, 0,     nV-1,      nV,    nV+k-1,    ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,  &HvV, nV,    nV+k-1,    0,     nV-1,      ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,  &Hvv, nV,    nV+k-1,    nV,    nV+k-1,    ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Q_H_,&Q_H, 0,     nV-1,      0,     nV-1,      ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (R_H_,&R_H, 0,     nV-1,      0,     nV-1,      ierr), *ierr);

    if( *nEig > 0 )     // only valid views after first converged ev
    {
      PHIST_CHK_IERR(SUBR( mvec_view_block  ) (Q_,  &Q,                     0,     *nEig-1,   ierr), *ierr);
      PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BQ_, &BQ,                    0,     *nEig-1,   ierr), *ierr);
      PHIST_CHK_IERR(SUBR( sdMat_view_block ) (R_,  &R,   0,     *nEig-1,   0,     *nEig-1,   ierr), *ierr);
      PHIST_CHK_IERR(SUBR( sdMat_view_block ) (R_,  &a,   0,     *nEig-1,   *nEig, *nEig+k-1, ierr), *ierr);
    }


    // copy H to Q and calculate sorted schur form of H
    PHIST_CHK_IERR(SUBR( sdMat_add_sdMat ) (st::one(), H, st::zero(), Q_H, ierr), *ierr);
    PHIST_CHK_IERR(SUBR( SchurDecomp ) (Q_H_raw, ldaQ_H, R_H_raw, ldaR_H, nV, nV, nV, which, ev_H, ierr), *ierr);
  }



  //------------------------------- delete vectors and matrices --------------------
  // delete views
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Hvv, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (HvV, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (HVv, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (H,   ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (r,   ierr), *ierr);

  PHIST_CHK_IERR(SUBR( mvec_delete  ) (BQQ, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (BQq, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (QQ,  ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (Qq,  ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (t,   ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (BVv, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (BV,  ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (AVv, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (AV,  ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (Vv,  ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (V,   ierr), *ierr);

  if( *nEig > 0 )
  {
    PHIST_CHK_IERR(SUBR( sdMat_delete ) (a,   ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_delete ) (R,   ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_delete  ) (BQ,  ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_delete  ) (Q,   ierr), *ierr);
  }

  delete[] ev_H;

  // delete mvecs and sdMats
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (H_,  ierr), *ierr);
  if( B_op != NULL )
  {
    PHIST_CHK_IERR(SUBR( mvec_delete )(BV_, ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_delete )(BQ_, ierr), *ierr);
  }
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (t_,  ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (AV_, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (V_,  ierr), *ierr);

}

