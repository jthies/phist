//! small helper function that checks if an sdMat is symmetric
void SUBR(sdMat_check_symmetrie)(TYPE(const_sdMat_ptr) mat, _MT_ tol, int*ierr)
{
  PHIST_ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"

  // check dimensions
  int m, n;
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(mat, &m, ierr), *ierr);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(mat, &n, ierr), *ierr);
  if( m != n )
  {
    *ierr = 1;
    return;
  }

  // create tmp storage
  TYPE(sdMat_ptr) tmp = NULL;
  PHIST_CHK_IERR(SUBR(sdMat_create)(&tmp, m, n, NULL, ierr), *ierr);

  // construct identity matrix
  TYPE(sdMat_ptr) id = NULL;
  PHIST_CHK_IERR(SUBR(sdMat_create)(&id, m, n, NULL, ierr), *ierr);
  {
    _ST_ *id_raw = NULL;
    lidx_t lda;
    PHIST_CHK_IERR(SUBR( sdMat_put_value ) (id, st::zero(), ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_extract_view) (id, &id_raw, &lda, ierr), *ierr);
    for(int i = 0; i < m; i++)
      id_raw[i*lda+i] = st::one();
  }


  // set tmp to transposed
  PHIST_CHK_IERR(SUBR(sdMatT_times_sdMat)(st::one(), mat, id, st::zero(), tmp, ierr), *ierr);
  // subtract mat
  PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(-st::one(), mat, st::one(), tmp, ierr), *ierr);
  // calc max. abs. value of mat^T-mat
  _MT_ maxVal = mt::zero();
  {
    _ST_ *tmp_raw = NULL;
    lidx_t lda;
    PHIST_CHK_IERR(SUBR(sdMat_extract_view)(tmp, &tmp_raw, &lda, ierr), *ierr);
    for(int i = 0; i < m; i++)
      for(int j = 0; j < n; j++)
        maxVal = std::max(maxVal, st::abs(tmp_raw[i*lda+j]));
  }
  PHIST_SOUT(PHIST_VERBOSE, "Symmetrie deviation of projection %e\n", maxVal);

  PHIST_CHK_IERR(SUBR(sdMat_delete)(id, ierr), *ierr);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(tmp, ierr), *ierr);

  if( maxVal < tol )
    *ierr = 0;
  else
    *ierr = 1;
}


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
//! innerBlockDim: block dimension used in the inner GMRES itersion
//! innerMaxBase:  restart inner GMRES after this number of iterations
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
void SUBR(subspacejada)( TYPE(const_op_ptr) A_op,  TYPE(const_op_ptr) B_op,
                         TYPE(const_mvec_ptr) v0,  eigSort_t which,
                         _MT_ tol,                 int nEig,
                         int* nIter,               int blockDim,
                         int minBase,              int maxBase,
                         int innerBlockDim,        int innerMaxBase,
                         int initialShiftIter,     _ST_ initialShift,
                         bool innerIMGS,           bool innerGMRESabortAfterFirstConverged,
                         bool symmetric,
                         TYPE(mvec_ptr) Q__,       TYPE(sdMat_ptr) R_,

                         _MT_* resNorm,            int* ierr)
{
  PHIST_ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *ierr = 0;

  // extended number of searched eigenvalues in order to respect the block dimension,
  // this way we always have a fixed blockDim AND it should make the calculation
  // of the last eigenvalues more stable in some cases
  int nEig_ = nEig + blockDim - 1;

  //------------------------------- check arguments --------------------------------
  if( blockDim < 1 )
  {
    PHIST_SOUT(PHIST_ERROR, "parameter blockDim < 1!\n");
    PHIST_CHK_IERR(*ierr = -99, *ierr);
  }
  if( innerBlockDim > blockDim || innerBlockDim < 1)
  {
    PHIST_SOUT(PHIST_ERROR, "parameter innerBlockDim > blockDim || innerBlockDim < 1!\n");
    PHIST_CHK_IERR(*ierr = -99, *ierr);
  }
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
/*
  if( minBase % innerBlockDim != 0 )
  {
    PHIST_SOUT(PHIST_WARNING, "minBase is not a multiple of innerBlockDim, switching to single-vector arnoldi for initiali subspace!\n");
  }
*/

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
  mvec_ptr_t  Q_      = NULL;    //< Q, enlarged dynamically
  mvec_ptr_t  BQ_     = NULL;    //< B*Q, enlarged dynamically
  mvec_ptr_t  t_      = NULL;    //< space for t
  mvec_ptr_t  At_     = NULL;    //< space for A*t
  mvec_ptr_t  res     = NULL;    //< residuum A*Q-Q*R

  sdMat_ptr_t H_      = NULL;    //< space for H
  sdMat_ptr_t Htmp_   = NULL;    //< temporary space for H
  sdMat_ptr_t Q_H_    = NULL;    //< space for Q_H
  sdMat_ptr_t R_H_    = NULL;    //< space for R_H
  sdMat_ptr_t sdMI_   = NULL;    //< identity matrix
  _ST_ sigma[nEig_];             //< JaDa correction shifts

  _ST_ *Q_H_raw       = NULL;
  _ST_ *R_H_raw       = NULL;
  _ST_ *Htmp_raw 			= NULL;
  lidx_t ldaQ_H, ldaR_H, ldaHtmp;

  PHIST_CHK_IERR(SUBR( mvec_create  ) (&V_,     A_op->domain_map, maxBase,        ierr), *ierr);
  // TODO: remove Vtmp
#ifdef TESTING
  PHIST_CHK_IERR(SUBR( mvec_create  ) (&Vtmp_,  A_op->domain_map, maxBase,                ierr), *ierr);
#endif
  PHIST_CHK_IERR(SUBR( mvec_create  ) (&AV_,    A_op->range_map,  maxBase,        	      ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_create  ) (&t_,     A_op->domain_map, blockDim,               ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_create  ) (&At_,    A_op->range_map,  blockDim,               ierr), *ierr);
  int resDim = std::min(2*blockDim, nEig+blockDim-1);
  PHIST_CHK_IERR(SUBR( mvec_create  ) (&res,    A_op->range_map,  resDim,                 ierr), *ierr);

  PHIST_CHK_IERR(SUBR( sdMat_create ) (&H_,     maxBase,          maxBase,  range_comm,   ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_create ) (&Htmp_,  maxBase,          maxBase,  range_comm,   ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_create ) (&Q_H_,   maxBase,          maxBase,  range_comm,   ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_create ) (&R_H_,   maxBase,          maxBase,  range_comm,   ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_create ) (&sdMI_,  maxBase,          maxBase,  range_comm,   ierr), *ierr);
  // construct identity matrix
  {
    _ST_ *sdMI_raw = NULL;
    lidx_t lda;
    PHIST_CHK_IERR(SUBR( sdMat_put_value ) (sdMI_, st::zero(), ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_extract_view) (sdMI_, &sdMI_raw, &lda, ierr), *ierr);
    for(int i = 0; i < maxBase; i++)
      sdMI_raw[i*lda+i] = st::one();
  }

  PHIST_CHK_IERR(SUBR( sdMat_extract_view ) (Q_H_,    &Q_H_raw,   &ldaQ_H,   ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_extract_view ) (R_H_,  	&R_H_raw,   &ldaR_H,   ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_extract_view ) (Htmp_,   &Htmp_raw,  &ldaHtmp,  ierr), *ierr);
  if( B_op != NULL )
  {
    PHIST_CHK_IERR(SUBR( mvec_create )(&BV_,    B_op->range_map,  maxBase,                ierr), *ierr);
  }
  else
  {
    BV_ = V_;
  }
  // array for the (possibly complex) eigenvalues for SchurDecomp
  CT* ev_H = new CT[maxBase];

  // create views on mvecs and sdMats with current dimensions
  int nV  = minBase;          //< current subspace dimension

  mvec_ptr_t  V   = NULL;     //< B-orthogonal basis of the search space
  mvec_ptr_t  Vful= NULL;     //< B-orthogonal basis of the search space + already locked Schur vectors
  mvec_ptr_t  Vtmp= NULL;     //< temporary V
  mvec_ptr_t  Vv  = NULL;     //< next columns in V_
  mvec_ptr_t  AV  = NULL;     //< A*V
  mvec_ptr_t  AVful= NULL;     //< A*Vful
  mvec_ptr_t  AVv = NULL;     //< next columns in AV_
  mvec_ptr_t  BV  = NULL;     //< B*V
  mvec_ptr_t  BVful  = NULL;     //< B*Vful
  mvec_ptr_t  BVv = NULL;     //< next columns in BV_
  mvec_ptr_t  t   = NULL;     //< Block-Jacobi-Davidson correction
  mvec_ptr_t  t_res = NULL;   //< part of the residual AQ-QR corresponding to current block t
  mvec_ptr_t  Qtil= NULL;     //< view of part of Q required for the JaDa correction equation
  mvec_ptr_t BQtil= NULL;     //< B*Qtil
  mvec_ptr_t Q = NULL, BQ = NULL, R = NULL;

  sdMat_ptr_t H   = NULL;     //< projection of A onto H, V'*AV
  sdMat_ptr_t Hful= NULL;     //< projection of A onto H, Vful'*A*Vful
  sdMat_ptr_t Hh  = NULL;     //< inside view for H
  sdMat_ptr_t Htmp= NULL;     //< temporary space for H
  sdMat_ptr_t HVv = NULL;     //< next rows in H_
  sdMat_ptr_t HvV = NULL;     //< next columns in H_
  sdMat_ptr_t Hvv = NULL;     //< next block on the diagonal of H_
  //sdMat_ptr_t Rr  = NULL;     //< [R a; 0 r]
  sdMat_ptr_t Q_H = NULL;     //< schur vectors of H
  sdMat_ptr_t Qq_H = NULL;
  mvec_ptr_t  Qq  = NULL;
  mvec_ptr_t BQq = NULL;
  sdMat_ptr_t R_H = NULL;     //< schur matrix of H
  sdMat_ptr_t Rr_H = NULL;
  sdMat_ptr_t sdMI = NULL;



  //------------------------------- initialize correction equation solver solver ------------------------
  TYPE(jadaCorrectionSolver_ptr) innerSolv = NULL;
  PHIST_CHK_IERR(SUBR(jadaCorrectionSolver_create)(&innerSolv, innerBlockDim, A_op->domain_map, GMRES, innerMaxBase, symmetric, ierr), *ierr);
  std::vector<_MT_> innerTol(nEig_,mt::one());
  std::vector<_MT_> lastOuterRes(nEig_,mt::zero());

  //------------------------------- initialize subspace etc ------------------------
  // run arnoldi
  nV = minBase;
  //TODO B_op in arnoldi
  // calculates A*V(:,1:m) = V(:,1:m+1)*H(1:m+1,1:m)
  // also outputs A*V (DON'T recalculate it from V*H, because this may not be accurate enough!)
/* // seems numerically less useful than unblocked arnoldi!
  if( blockDim > 0 && nV % innerBlockDim == 0 )
  {
    int bs = innerBlockDim;
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (V_,      &V,                       0,     nV+bs-1,   ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (AV_,     &AV,                      0,     nV+bs-1,   ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BV_,     &BV,                      0,     nV+bs-1,   ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,      &H,     0,      nV+bs-1,  0,     nV-1,      ierr), *ierr);

    PHIST_CHK_IERR(SUBR( simple_blockArnoldi ) (A_op, B_op, V, AV, BV, H, nV, innerBlockDim, ierr), *ierr);
  }
  else
*/
  {
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (V_,      &V,                       0,     nV,        ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (AV_,     &AV,                      0,     nV,        ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BV_,     &BV,                      0,     nV,        ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,      &H,     0,      minBase,  0,     nV-1,      ierr), *ierr);

    PHIST_CHK_IERR(SUBR( simple_arnoldi ) (A_op, B_op, v0, V, AV, BV, H, nV, ierr), *ierr);
  }

  // set views
  int nConvEig = 0;
#define UPDATE_SUBSPACE_VIEWS \
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (V_,      &V,                         nConvEig, nV-1,      ierr), *ierr); \
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (AV_,     &AV,                        nConvEig, nV-1,      ierr), *ierr); \
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BV_, 		&BV,                        nConvEig, nV-1,      ierr), *ierr); \
  PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,  		&H,     nConvEig, nV-1,     nConvEig, nV-1,      ierr), *ierr); \
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (V_,      &Vful,                      0,        nV-1,      ierr), *ierr); \
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (AV_,     &AVful,                     0,        nV-1,      ierr), *ierr); \
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BV_,     &BVful,                     0,        nV-1,      ierr), *ierr); \
  PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,      &Hful,  0,        nV-1,     0,        nV-1,      ierr), *ierr);

  UPDATE_SUBSPACE_VIEWS;



#ifdef TESTING
#define TESTING_CHECK_SUBSPACE_INVARIANTS \
{ \
  /* check orthogonality of V, BV, Q */ \
  PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Htmp_,&Htmp,0,    nV-1,      0,     nV-1,      ierr), *ierr); \
  PHIST_CHK_IERR(SUBR( mvecT_times_mvec ) (st::one(), Vful, BVful, st::zero(), Htmp, ierr), *ierr); \
  _MT_ orthEps = std::abs(Htmp_raw[0] - st::one()); \
  for(int i = 0; i < nV; i++) \
    for(int j = 0; j < nV; j++) \
      orthEps = std::max(orthEps, std::abs(Htmp_raw[i*ldaHtmp+j] - ((i==j) ? st::one() : st::zero()))); \
  PHIST_OUT(PHIST_INFO, "Line %d: B-orthogonality of V: %e\n", __LINE__, orthEps); \
 \
  /* check AV = A*V */ \
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (Vtmp_, &Vtmp,                  0,       nV-1,          ierr), *ierr); \
  PHIST_CHK_IERR( A_op->apply(st::one(), A_op->A, Vful, st::zero(), Vtmp, ierr), *ierr); \
  PHIST_CHK_IERR(SUBR( mvec_add_mvec ) (-st::one(), AVful, st::one(), Vtmp, ierr), *ierr); \
  _MT_ normVtmp[nV]; \
  PHIST_CHK_IERR(SUBR( mvec_norm2 ) (Vtmp, normVtmp, ierr), *ierr); \
  _MT_ equalEps = normVtmp[0]; \
  for(int i = 0; i < nV; i++) \
    equalEps = std::max(equalEps, normVtmp[i]); \
  PHIST_OUT(PHIST_INFO, "Line %d: AV - A*V: %e\n", __LINE__, equalEps); \
  /* check H = V'*A*V */ \
  PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Htmp_,&Htmp,0,    nV-1,      0,     nV-1,      ierr), *ierr); \
  PHIST_CHK_IERR(SUBR( mvecT_times_mvec ) (st::one(), Vful, AVful, st::zero(), Htmp, ierr), *ierr); \
  PHIST_CHK_IERR(SUBR( sdMat_add_sdMat ) (-st::one(), Hful, st::one(), Htmp, ierr), *ierr); \
  equalEps = std::abs(Htmp_raw[0]); \
  for(int i = 0; i < nV; i++) \
    for(int j = 0; j < nV; j++) \
      equalEps = std::max(equalEps, std::abs(Htmp_raw[i*ldaHtmp+j])); \
  PHIST_OUT(PHIST_INFO, "Line %d: H - V'*AV: %e\n", __LINE__, equalEps); \
  /* PHIST_SOUT(PHIST_INFO, "H:\n"); */ \
  /* PHIST_CHK_IERR(SUBR( sdMat_print )(H, ierr), *ierr); */ \
  /* PHIST_SOUT(PHIST_INFO, "H - V'*AV:\n"); */ \
  /* PHIST_CHK_IERR(SUBR( sdMat_print )(Htmp, ierr), *ierr); */ \
}
#else
#define TESTING_CHECK_SUBSPACE_INVARIANTS
#endif


TESTING_CHECK_SUBSPACE_INVARIANTS;


  //----------------------------------- MAIN LOOP ----------------------------------
  int maxIter = *nIter;
  int Qsize = 0;
  for(*nIter = 0; *nIter < maxIter; (*nIter)++)
  {
// dynamically adjust current number of "sought" eigenvalues, so we do not consider the 20th eigenvalue if the 1st is not calculated yet!
nEig_ = std::min(nConvEig + 2*blockDim, nEig+blockDim-1);
// dynamically adjust the size of the buffer for Q_, so we don't work with a huge stride at the beginning of the calculation!
if( Qsize < nEig_ )
{
  mvec_ptr_t newQ = NULL;
  int newQsize = std::min(std::max(nEig_,2*Qsize), nEig+blockDim-1);
  if( newQsize % 2 != 0 )
    newQsize++;
  while( newQsize % innerBlockDim != 0 )
    newQsize++;
  PHIST_CHK_IERR(SUBR(mvec_create)(&newQ, A_op->range_map, newQsize, ierr), *ierr);
  if( Qsize > 0 )
  {
    PHIST_CHK_IERR(SUBR(mvec_set_block)(newQ, Q_, 0, Qsize-1, ierr), *ierr);
  }
  PHIST_CHK_IERR(SUBR(mvec_delete)(Q_, ierr), *ierr);
  Q_ = newQ;

  mvec_ptr_t newBQ = NULL;
  if( B_op != NULL )
  {
    PHIST_CHK_IERR(SUBR(mvec_create)(&newBQ, A_op->range_map, newQsize, ierr), *ierr);
    if( Qsize > 0 )
    {
      PHIST_CHK_IERR(SUBR(mvec_set_block)(newBQ, BQ_, 0, Qsize-1, ierr), *ierr);
    }
    PHIST_CHK_IERR(SUBR(mvec_delete)(Q_, ierr), *ierr);
    BQ_ = newBQ;
  }
  else
  {
    BQ_ = Q_;
  }
  Qsize = newQsize;
}
// update views
PHIST_CHK_IERR(SUBR( mvec_view_block ) (Q_,  &Q,  0, nEig_-1, ierr), *ierr);
PHIST_CHK_IERR(SUBR( mvec_view_block ) (BQ_, &BQ, 0, nEig_-1, ierr), *ierr);
PHIST_CHK_IERR(SUBR( sdMat_view_block ) (R_,  &R, 0, nEig_-1, 0, nEig_-1, ierr), *ierr);

    // for convenience check symmetrie of H (shouldn't much hurt the performance)
    if( symmetric )
    {
      PHIST_CHK_IERR(SUBR(sdMat_check_symmetrie)(Hful, tol, ierr), *ierr);
    }

    // calculate sorted Schur form of H in (Q_H,R_H)
    // we only update part of Q_H,R_H, so first set Q_H, R_H to zero
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Q_H_,&Q_H, 0,     nV-1,      0,     nV-1,      ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (R_H_,&R_H, 0,     nV-1,      0,     nV-1,      ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_put_value  ) (Q_H, st::zero(), ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_put_value  ) (R_H, st::zero(), ierr), *ierr);
    // then copy the new block of H
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (R_H_,&R_H, nConvEig, nV-1, nConvEig, nV-1, ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_get_block  ) (H_,   R_H,  nConvEig, nV-1, nConvEig, nV-1, ierr), *ierr);
    int nSort = minBase-nConvEig; //nEig_-nConvEig;
    int nSelect = nSort;
    lidx_t offR_H = ldaR_H*nConvEig+nConvEig;
    lidx_t offQ_H = ldaQ_H*nConvEig+nConvEig;
    PHIST_CHK_IERR(SUBR( SchurDecomp ) (R_H_raw+offR_H, ldaR_H, Q_H_raw+offQ_H, ldaQ_H, nV-nConvEig, nSelect, nSort, which, tol, ev_H+nConvEig, ierr), *ierr);
    // we still need to add the missing parts of R_H, Q_H
    if( nConvEig > 0 )
    {
      // upper left part of Q_H
      for(int i = 0; i < nConvEig; i++)
        Q_H_raw[ldaQ_H*i+i] = st::one();

      // left part of R_H
      PHIST_CHK_IERR(SUBR( sdMat_view_block ) (R_H_, &R_H, 0, nConvEig-1, 0, nConvEig-1, ierr), *ierr);
      PHIST_CHK_IERR(SUBR( sdMat_get_block  ) (R,    R_H,  0, nConvEig-1, 0, nConvEig-1, ierr), *ierr);

      // upper right part of R_H
      PHIST_CHK_IERR(SUBR( sdMat_view_block  ) (Q_H_, &Qq_H, nConvEig, nV-1,            nConvEig, nEig_-1, ierr), *ierr);
      PHIST_CHK_IERR(SUBR( sdMat_view_block  ) (R_H_, &Rr_H, 0,             nConvEig-1, nConvEig, nEig_-1, ierr), *ierr);
      PHIST_CHK_IERR(SUBR( sdMat_view_block  ) (H_  , &Hh,   0,             nConvEig-1, nConvEig, nV-1,    ierr), *ierr);

      PHIST_CHK_IERR(SUBR( sdMat_times_sdMat ) (st::one(), Hh, Qq_H, st::zero(), Rr_H, ierr), *ierr);
    }
#ifdef TESTING
{
  // check that H Q_H = Q_H R_H
  PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Q_H_,&Q_H, 0, nV-1,    0, nEig_-1, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_view_block ) (R_H_,&R_H, 0, nEig_-1, 0, nEig_-1, ierr), *ierr);
  PHIST_CHK_IERR(SUBR(sdMat_view_block)(Htmp_, &Htmp, 0, nV-1,    0, nEig_-1, ierr), *ierr);
  PHIST_CHK_IERR(SUBR(sdMat_times_sdMat)(st::one(), Hful, Q_H, st::zero(), Htmp, ierr), *ierr);
  PHIST_CHK_IERR(SUBR(sdMat_times_sdMat)(-st::one(), Q_H, R_H, st::one(), Htmp, ierr), *ierr);
  // get max norm
  _MT_ absErr = mt::zero();
  for(int i = 0; i < nEig_; i++)
    for(int j = 0; j < nV; j++)
      absErr = std::max(absErr, st::abs(Htmp_raw[i*ldaHtmp+j]));
  PHIST_SOUT(PHIST_INFO, "H*Q_H - Q_H*R_H: %8.4e\n", absErr);
  //PHIST_CHK_IERR(SUBR(sdMat_print)(H, ierr), *ierr);
  //PHIST_CHK_IERR(SUBR(sdMat_print)(Q_H, ierr), *ierr);
  //PHIST_CHK_IERR(SUBR(sdMat_print)(R_H, ierr), *ierr);
  //PHIST_CHK_IERR(SUBR(sdMat_print)(Htmp, ierr), *ierr);
}
#endif

    // update views
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Q_H_,&Q_H, 0,     nV-1,      0,             nV-1,      ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (R_H_,&R_H, 0,     nV-1,      0,             nV-1,      ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Q_H_,&Qq_H,nConvEig,nV-1,      nConvEig, nEig_-1,   ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (R_H_,&Rr_H,0,     nEig_-1,   nConvEig, nEig_-1,   ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (Q_,   &Qq,                   nConvEig, nEig_-1,   ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BQ_,  &BQq,                  nConvEig, nEig_-1,   ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (res, &t_res,                 0, nEig_-nConvEig-1,   ierr), *ierr);

    // update approximate Schur form of A (keeping the already computed part locked)
    PHIST_CHK_IERR(SUBR( mvec_times_sdMat ) (st::one(), V,    Qq_H, st::zero(), Qq,  ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_set_block ) (Q, Qq, nConvEig, nEig_-1, ierr), *ierr);
    PHIST_CHK_IERR(SUBR( sdMat_set_block  ) (R, Rr_H, 0, nEig_-1, nConvEig, nEig_-1, ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_times_sdMat ) (st::one(), AV,   Qq_H, st::zero(), t_res, ierr), *ierr);
    if( B_op != NULL )
    {
      PHIST_CHK_IERR(SUBR( mvec_times_sdMat ) (st::one(), BV,   Qq_H, st::zero(), BQq,  ierr), *ierr);
      PHIST_CHK_IERR(SUBR( mvec_set_block ) (BQ, BQq, nConvEig, nEig_-1, ierr), *ierr);
    }
    // overwrite res with the residual: -res = -(Aq - Bqr) = + BQq*r - res
    PHIST_CHK_IERR(SUBR( mvec_times_sdMat ) (st::one(), BQ,   Rr_H,    -st::one(), t_res, ierr), *ierr);
    // calculate norm of the residual
    PHIST_CHK_IERR(SUBR( mvec_norm2 ) (t_res, resNorm+nConvEig, ierr), *ierr);
#ifdef TESTING
{
  // check that the residual is orthogonal to Q (should be by construction!)
  PHIST_CHK_IERR( SUBR( sdMat_view_block ) (Htmp_, &Htmp, 0, nEig_-1, nConvEig, nEig_-1, ierr), *ierr);
  PHIST_CHK_IERR( SUBR( mvecT_times_mvec ) (st::one(), Q, t_res, st::zero(), Htmp, ierr), *ierr);
  _MT_ equalEps = std::abs(Htmp_raw[0]);
  for(int i = 0; i < nEig_; i++)
    for(int j = nConvEig; j < nEig_; j++)
      equalEps = std::max(equalEps, std::abs(Htmp_raw[i*ldaHtmp+j]));
  PHIST_OUT(PHIST_INFO, "Res orthogonality wrt. Q: %e\n", equalEps);
  // calculate explicit residual
  // AQ - BQR
  PHIST_CHK_IERR( SUBR( mvec_view_block ) (Vtmp_, &Vtmp, 0, nEig_-1, ierr), *ierr);
  PHIST_CHK_IERR( A_op->apply(st::one(), A_op->A, Q, st::zero(), Vtmp, ierr), *ierr);
  PHIST_CHK_IERR( SUBR( mvec_times_sdMat ) (-st::one(), BQ, R, st::one(), Vtmp, ierr), *ierr);
  _MT_ expRes[nEig_];
  PHIST_CHK_IERR( SUBR( mvec_norm2 )(Vtmp, expRes, ierr), *ierr);
  PHIST_SOUT(PHIST_INFO, "est. residual: ");
  for(int i = 0; i < nEig_; i++)
    PHIST_SOUT(PHIST_INFO, "\t%8.4e", resNorm[i]);
  PHIST_SOUT(PHIST_INFO, "\nexp. residual: ");
  for(int i = 0; i < nEig_; i++)
    PHIST_SOUT(PHIST_INFO, "\t%8.4e", expRes[i]);
  PHIST_SOUT(PHIST_INFO, "\n");
}
#endif


    // reorder multiple eigenvalues in schur form by residual norm
    // only possible/implemented for the symmetric case
    std::vector<int> resPermutation(nEig_);
    for(int i = 0; i < nEig_; i++)
      resPermutation[i] = i;
    if( false ) // symmetric ) // doesn't work always?
    {
      PHIST_CHK_IERR( SUBR(ReorderPartialSchurDecomp)(R_H_raw+offR_H, ldaR_H, Q_H_raw+offQ_H, ldaQ_H, nV-nConvEig, nEig_-nConvEig, which, sqrt(tol), resNorm+nConvEig, ev_H+nConvEig, &resPermutation[nConvEig], ierr), *ierr);
      for(int i = nConvEig; i < nEig_; i++)
        resPermutation[i] += nConvEig;
#ifdef TESTING
PHIST_SOUT(PHIST_INFO,"resPermutation: ");
for(int i = 0; i < nEig_; i++)
  PHIST_SOUT(PHIST_INFO,"\t%d", resPermutation[i]);
PHIST_SOUT(PHIST_INFO,"\n");
#endif
    }
    // check if we need to adapt Q to the new ordering (e.g. there were duplicate eigenvalues not sorted by their residual norm)
    // res is handled implicitly later
    for(int i = 0; i < nEig_; i++)
    {
      if( resPermutation[i] != i )
      {
        // setup permutation matrix
        //PHIST_CHK_IERR(SUBR(sdMat_view_block)(Htmp_, &Htmp, i, nEig_-1, i, nEig_-1, ierr), *ierr);
        //PHIST_CHK_IERR(SUBR(sdMat_put_value)(Htmp, st::zero(), ierr), *ierr);
        //for(int j = i; j < nEig_; j++)
        //{
          //int j_ = resPermutation[j];
          //Htmp_raw[j_+j*ldaHtmp] = st::one();
        //}
        //PHIST_CHK_IERR(SUBR(mvec_view_block)(Q_, &Qq, i, nEig_-1, ierr), *ierr);
        //PHIST_CHK_IERR(SUBR(mvec_times_sdMat_inplace)(Qq, Htmp, ierr), *ierr);
        //if( B_op != NULL )
        //{
          //PHIST_CHK_IERR(SUBR(mvec_view_block)(BQ_, &BQq, i, nEig_-1, ierr), *ierr);
          //PHIST_CHK_IERR(SUBR(mvec_times_sdMat_inplace)(BQq, Htmp, ierr), *ierr);
        //}
        PHIST_CHK_IERR(SUBR( mvec_times_sdMat ) (st::one(), V,    Qq_H, st::zero(), Qq,  ierr), *ierr);
        PHIST_CHK_IERR(SUBR( sdMat_set_block  ) (R, Rr_H, 0, nEig_-1, nConvEig, nEig_-1, ierr), *ierr);
        break;
      }
    }
#ifdef TESTING
{
  // check that Q is in correct order
  PHIST_CHK_IERR( SUBR(mvec_view_block)(Vtmp_, &t, 0, nEig_-1, ierr), *ierr);
  PHIST_CHK_IERR(A_op->apply(st::one(), A_op->A, Q, st::zero(), t, ierr), *ierr);
  PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(-st::one(), BQ, R, st::one(), t, ierr), *ierr);
  _MT_ reorderedResNorm[nEig_];
  PHIST_CHK_IERR(SUBR(mvec_norm2)(t, reorderedResNorm, ierr), *ierr);
  PHIST_SOUT(PHIST_INFO,"est. residual norm after reordering:");
  for(int i = 0; i < nEig_; i++)
    PHIST_SOUT(PHIST_INFO,"\t%8.4e",resNorm[i]);
  PHIST_SOUT(PHIST_INFO,"\nexp. residual norm after reordering:");
  bool err = false;
  for(int i = 0; i < nEig_; i++)
  {
    PHIST_SOUT(PHIST_INFO,"\t%8.4e",reorderedResNorm[i]);
    // we don't want it to abort in case of some inaccuracies...
    if( mt::abs(reorderedResNorm[i] - resNorm[i]) > mt::sqrt(tol) )
      err = true;
  }
  PHIST_SOUT(PHIST_INFO,"\n");
  if( err )
  {
    PHIST_SOUT(PHIST_WARNING,"strong deviation of estimated and explicit residuals (see above)!\n");
  }
}
#endif

    // check for converged eigenvalues
    int nNewConvEig = 0;
    for(int i = nConvEig; i < nEig; i++)
    {
      PHIST_SOUT(PHIST_INFO,"In iteration %d: Current approximation for eigenvalue %d is %16.8g%+16.8gi with residuum %e\n", *nIter, i+1, ct::real(ev_H[i]),ct::imag(ev_H[i]), resNorm[i]);
#ifndef IS_COMPLEX
      if( std::abs(ct::imag(ev_H[i])) > tol && blockDim == 1 )
      {
        PHIST_SOUT(PHIST_WARNING, "Detected possible complex-conjugate eigenpair, but blockDim == 1 (you need at least blockDim=2 to detect complex-conjugate eigenpairs correctly)!\n");
      }
#endif
      if( resNorm[i] <= tol && i == nConvEig+nNewConvEig )
      {
#ifndef IS_COMPLEX
        // detect complex conjugate eigenvalue pairs
        if( std::abs(ct::imag(ev_H[i])) > tol && i+1 < nEig_ )
        {
          if( ct::abs(ct::conj(ev_H[i]) - ev_H[i+1]) < sqrt(tol) && resNorm[i+1] <= tol )
            nNewConvEig+=2;
          continue;
        }
#endif
        nNewConvEig++;
      }
      else if( blockDim == 1 )
        break;
    }

    if( nNewConvEig > 0 )
    {
      PHIST_SOUT(PHIST_INFO,"In iteration %d: locking %d newly converged eigenvalues\n", *nIter, nNewConvEig);
    }

    if( nConvEig + nNewConvEig >= nEig )
    {
      nConvEig += nNewConvEig;
      PHIST_SOUT(PHIST_INFO,"In iteration %d: all eigenvalues converged!\n", *nIter);
      break;
    }

    if( *nIter >= maxIter )
    {
      PHIST_SOUT(PHIST_INFO,"Reached maximum number of iterations!\n");
      break;
    }

    if( nNewConvEig > 0 )
    {
      // to avoid unnecessary subspace transformations, shrink the searchspace one iteration earlier...
      if( nV + 2*blockDim > maxBase )
      {
        PHIST_SOUT(PHIST_INFO,"Shrinking search space (one iteration earlier) from %d to %d\n", nV, minBase);
        if( symmetric )
        {
          PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Q_H_,  &Q_H,  nConvEig, nV-1,          nConvEig, minBase-1,    ierr), *ierr);
        }
        else
        {
          PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Q_H_,  &Q_H,  0, nV-1,          0, minBase-1,    ierr), *ierr);
        }
      }
      else
      {
        if( symmetric )
        {
          PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Q_H_,  &Q_H,  nConvEig, nV-1,          nConvEig, nV-1,    ierr), *ierr);
        }
        else
        {
          PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Q_H_,  &Q_H,  0, nV-1,          0, nV-1,    ierr), *ierr);
        }
      }

      // reorder V and H
      if( symmetric )
      {
        PHIST_CHK_IERR(SUBR( transform_searchSpace ) (V, AV, BV, H, Q_H, B_op != NULL, ierr), *ierr);
      }
      else
      {
        PHIST_CHK_IERR(SUBR( transform_searchSpace ) (Vful, AVful, BVful, Hful, Q_H, B_op != NULL, ierr), *ierr);
      }

      nConvEig = nConvEig+nNewConvEig;

      // update views if necessary
      if( nV + 2*blockDim > maxBase )
      {
        nV = minBase;

        UPDATE_SUBSPACE_VIEWS;
      }
TESTING_CHECK_SUBSPACE_INVARIANTS;
    }



    // setup matrix of shifts and residuals for the correction equation
    int k_ = 0; // 0:k_ vectors of Q used for the orthogonal projection in the correction equation
    int k = 0;  // is always <= blockDim!
    std::vector<int> selectedRes(blockDim);
    for(int i = 0; i < nEig_ && k < blockDim; i++)
    {
      // only allow i >= nEig for multiple eigenvalues and complex conjugated eigenpairs
      if( i >= nEig )
      {
        if( ct::abs(ev_H[i] - ev_H[i-1]) > mt::sqrt(tol) && ct::abs(ct::conj(ev_H[i]) - ev_H[i-1]) > mt::sqrt(tol) )
          break;
      }

      if( resNorm[i] > tol )
      {
        k_ = i;

        if( *nIter < initialShiftIter )
        {
          sigma[k] = initialShift;
        }
        else
        {
#ifndef IS_COMPLEX
          sigma[k] = ct::real(ev_H[i]);
#else
          sigma[k] = ev_H[i];
#endif
        }

        // select the correct (unpermuted) residual
        selectedRes[k] = resPermutation[i];
        k++;
      }
    }
    // deflate with more vectors if there are multiple, partly converged eigenvalues
    while( k_+1 < nEig_ && (ct::abs(ev_H[k_+1]-ev_H[k_]) < 10*ct::abs(ev_H[k_+1])*mt::sqrt(tol) || ct::abs(ct::conj(ev_H[k_+1])-ev_H[k_]) < mt::sqrt(tol)) )
      k_++;

PHIST_SOUT(PHIST_INFO,"selectedRes: ");
for(int i = 0; i < k; i++)
  PHIST_SOUT(PHIST_INFO,"\t%d (%e)", selectedRes[i], resNorm[selectedRes[i]]);
PHIST_SOUT(PHIST_INFO,"\n");


    // shrink search space if necessary
    if( nV + k > maxBase )
    {
      PHIST_SOUT(PHIST_INFO,"Shrinking search space from %d to %d\n", nV, minBase);

      // nothing to do if converged eigenvalues this iteration
      if( nNewConvEig == 0 )
      {
        if( symmetric )
        {
          PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Q_H_,  &Q_H,  nConvEig, nV-1,          nConvEig, minBase-1,    ierr), *ierr);

          PHIST_CHK_IERR(SUBR( transform_searchSpace ) (V, AV, BV, H, Q_H, B_op != NULL, ierr), *ierr);
        }
        else
        {
          PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Q_H_,  &Q_H,  0, nV-1,          0, minBase-1,    ierr), *ierr);

          PHIST_CHK_IERR(SUBR( transform_searchSpace ) (Vful, AVful, BVful, Hful, Q_H, B_op != NULL, ierr), *ierr);
        }
      }

      nV = minBase;

      UPDATE_SUBSPACE_VIEWS;
TESTING_CHECK_SUBSPACE_INVARIANTS;
    }

    if( k > 0 ) 
    {
      // calculate corrections
      // setup jadaOp
      // set correction views and temporary jadaOp-storage
      PHIST_CHK_IERR(SUBR( mvec_view_block  ) (t_,  &t,     0, k-1,  ierr), *ierr);
      // we only need to view first part of Q
      PHIST_CHK_IERR(SUBR( mvec_view_block  ) (Q_,   &Qtil,  0, k_, ierr), *ierr);
      PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BQ_,  &BQtil, 0, k_, ierr), *ierr);
      // set tolerances
      for(int i = 0; i < k; i++)
      {
        if( resNorm[nConvEig+i] > 4*lastOuterRes[nConvEig+i] )
          innerTol[nConvEig+i] = 1.;
        innerTol[nConvEig+i] *= 0.5;
        lastOuterRes[nConvEig+i] = resNorm[nConvEig+i];
      }

      for(int i = 0; i < blockDim; i++)
        selectedRes[i] -= nConvEig-nNewConvEig;
      PHIST_CHK_NEG_IERR(SUBR(jadaCorrectionSolver_run)(innerSolv, A_op, B_op, Qtil, BQtil, sigma, res, &selectedRes[0],
                                                    &innerTol[nConvEig], innerMaxBase, t, innerIMGS, innerGMRESabortAfterFirstConverged, ierr), *ierr);

      // get solution and reuse res for At
      PHIST_CHK_IERR(SUBR( mvec_view_block  ) (t_, &Vv,  0, k-1, ierr), *ierr);
      PHIST_CHK_IERR(SUBR( mvec_view_block  ) (At_,&AVv, 0, k-1, ierr), *ierr);

      // enlarge search space
      // first update views
      PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,  &HVv, 0,     nV-1,      nV,    nV+k-1,    ierr), *ierr);
      PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,  &HvV, nV,    nV+k-1,    0,     nV-1,      ierr), *ierr);
      PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,  &Hvv, nV,    nV+k-1,    nV,    nV+k-1,    ierr), *ierr);
      // orthogonalize t as Vv (reuse R_H)
      PHIST_CHK_IERR(SUBR( sdMat_view_block ) (R_H_,&R_H, 0,     nV-1,      0,     k-1,       ierr), *ierr);
      PHIST_CHK_IERR(SUBR( sdMat_view_block ) (R_H_,&Rr_H,nV,    nV+k-1,    nV,    nV+k-1,    ierr), *ierr);
      int rankV;
      PHIST_CHK_NEG_IERR(SUBR( orthog ) (Vful, Vv, Rr_H, R_H, 5, &rankV, ierr), *ierr);
      // TODO: only take non-random vector if *ierr > 0
      // calculate AVv, BVv
      PHIST_CHK_IERR( A_op->apply(st::one(), A_op->A, Vv, st::zero(), AVv, ierr), *ierr);
      if( B_op != NULL )
      {
        PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BV_, &BVv,                   nV,    nV+k-1,    ierr), *ierr);
        PHIST_CHK_IERR( B_op->apply(st::one(), B_op->A, Vv, st::zero(), BVv, ierr), *ierr);
      }
      // update H
      PHIST_CHK_IERR(SUBR( mvecT_times_mvec ) (st::one(), Vful,  AVv, st::zero(), HVv, ierr), *ierr);
      // for the symmetric case use AVv*V here, so we don't need AV at all
      if( !symmetric )
      {
        PHIST_CHK_IERR(SUBR( mvecT_times_mvec ) (st::one(), Vv, AVful,  st::zero(), HvV, ierr), *ierr);
      }
      else
      {
        PHIST_CHK_IERR(SUBR(sdMat_view_block) (sdMI_, &sdMI, 0, nV-1, 0, nV-1, ierr), *ierr);
        PHIST_CHK_IERR(SUBR(sdMatT_times_sdMat)(st::one(), HVv, sdMI, st::zero(), HvV, ierr), *ierr);
      }
      PHIST_CHK_IERR(SUBR( mvecT_times_mvec ) (st::one(), Vv, AVv, st::zero(), Hvv, ierr), *ierr);
      // use set block to put Vv and AVv really into V and AV
      PHIST_CHK_IERR(SUBR( mvec_set_block ) (V_,  Vv,  nV, nV+k-1, ierr), *ierr);
      PHIST_CHK_IERR(SUBR( mvec_set_block ) (AV_, AVv, nV, nV+k-1, ierr), *ierr);
      // increase nV
      nV = nV + k;
    } // k > 0

    UPDATE_SUBSPACE_VIEWS;
TESTING_CHECK_SUBSPACE_INVARIANTS;
  }

  // copy result to Q_
  PHIST_CHK_IERR(SUBR(mvec_set_block)(Q__, Q, 0, nEig_-1, ierr), *ierr);


  //------------------------------- delete vectors and matrices --------------------
  PHIST_CHK_IERR(SUBR(jadaCorrectionSolver_delete)(innerSolv, ierr), *ierr);

  // delete views
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (R_H, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Rr_H,ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Q_H, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Qq_H,ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Hvv, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (HvV, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (HVv, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (H,   ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Hful,ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Htmp,ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (sdMI,ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (R,   ierr), *ierr);

  PHIST_CHK_IERR(SUBR( mvec_delete  ) (Qq,  ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (BQq, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (t,   ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (t_res,ierr),*ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (BVv, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (BV,  ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (AVv, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (AV,  ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (Vv,  ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (Vtmp,ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (V,   ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (Vful,ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (AVful,ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (BVful,ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (Qtil,ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (BQtil,ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (Q,   ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (BQ,  ierr), *ierr);

  delete[] ev_H;

  // delete mvecs and sdMats
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Q_H_,ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (R_H_,ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Htmp_,ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (H_,  ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (sdMI_,ierr), *ierr);
  if( B_op != NULL )
  {
    PHIST_CHK_IERR(SUBR( mvec_delete )(BV_, ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_delete )(BQ_, ierr), *ierr);
  }
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (Q_,  ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (t_,  ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (At_, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (res, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (AV_, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (Vtmp_,ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (V_,  ierr), *ierr);
}

