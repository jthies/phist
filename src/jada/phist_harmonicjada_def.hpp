
//! Tries to compute a partial schur form $(Q,R)$ of dimension opts.numEigs
//! of the stencil $A*x-\lambda*B*x$ with a general linear operator $A$ and a
//! hermitian positive definite (hpd.) linear operator $B$ using a
//! block-Jacobi-Davidson QR method based on harmonic Ritz values. <br>
//! The generalized eigenvalues $\lambda_i$ are the diagonal entries of the
//! partial schur form $A*Q = B*Q*R$ returned. <br>
//!
//! see header file for further documentation of the parameters
//!
void SUBR(harmonicjada)( TYPE(const_op_ptr) A_op,  TYPE(const_op_ptr) B_op,
                         phist_jadaOpts_t opts,
                         TYPE(mvec_ptr) Q__,       TYPE(sdMat_ptr) R_,
                         _CT_* ev,                 _MT_* resNorm,
                         int* nConv,               int* nIter,
                         int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag = 0;

  // copy options
  TYPE(const_mvec_ptr) v0=opts.v0;
  eigSort_t which=opts.which;
  _MT_ tol=opts.convTol;
  int nEig=opts.numEigs;
  int blockDim=opts.blockSize;
  int maxIter =opts.maxIters;
                         
int minBase=opts.minBas;
int maxBase=opts.maxBas;
int innerBlockDim=opts.innerSolvBlockSize;        
int innerMaxBase=opts.innerSolvMaxBas;
int initialShiftIter=opts.initialShiftIters;   
_ST_ initialShift   =(_ST_)opts.initialShift_r;
                    +(_ST_)opts.initialShift_i*st::cmplx_I();
                         
bool innerIMGS=(opts.innerSolvRobust!=0);
bool innerGMRESabortAfterFirstConverged=opts.innerSolvStopAfterFirstConverged;
bool symmetric=opts.symmetry==HERMITIAN;
#ifndef IS_COMPLEX
symmetric=symmetric||(opts.symmetry==COMPLEX_SYMMETRIC);
#endif

  eigExtr_t how=opts.how;
  if (how==STANDARD)
  {
    PHIST_SOUT(PHIST_ERROR,"if you want to use standard Ritz values, please use the subspacejada routine instead\n");
  }
  if (how!=HARMONIC)
  {
    PHIST_SOUT(PHIST_ERROR,"only Harmonic Ritz extraction is implemented (jadaOpts.how=%s), found %s\n",
        eigExtr2str(HARMONIC),eigExtr2str(how));
    *iflag=PHIST_NOT_IMPLEMENTED;
    return;
  }


  // extended number of searched eigenvalues in order to respect the block dimension,
  // this way we always have a fixed blockDim AND it should make the calculation
  // of the last eigenvalues more stable in some cases
  int nEig_ = nEig + blockDim - 1;

  // initialize residual norms to -1 to indicate that they haven't been computed
  for (int i=0;i<nEig; i++) resNorm[i]=-mt::one();

  //------------------------------- check arguments --------------------------------
  if( blockDim < 1 )
  {
    PHIST_SOUT(PHIST_ERROR, "parameter blockDim < 1!\n");
    PHIST_CHK_IERR(*iflag = PHIST_NOT_IMPLEMENTED, *iflag);
  }
  if( innerBlockDim > blockDim || innerBlockDim < 1)
  {
    PHIST_SOUT(PHIST_ERROR, "parameter innerBlockDim > blockDim || innerBlockDim < 1!\n");
    PHIST_CHK_IERR(*iflag = PHIST_NOT_IMPLEMENTED, *iflag);
  }
  if( minBase < nEig_ )
  {
    PHIST_SOUT(PHIST_ERROR, "parameter minBase < nEig+blockDim-1!\n");
    PHIST_CHK_IERR(*iflag = PHIST_NOT_IMPLEMENTED, *iflag);
  }
  if( minBase+blockDim > maxBase )
  {
    PHIST_SOUT(PHIST_ERROR, "parameter minBase+blockDim > maxBase!\n");
    PHIST_CHK_IERR(*iflag = PHIST_NOT_IMPLEMENTED, *iflag);
  }
  if( maxBase < nEig+blockDim )
  {
    PHIST_SOUT(PHIST_ERROR, "paramater maxBase < nEig+blockDim!\n");
    PHIST_CHK_IERR(*iflag = PHIST_NOT_IMPLEMENTED, *iflag);
  }
  if( B_op != NULL )
  {
    PHIST_SOUT(PHIST_ERROR,"case B_op != NULL (e.g. B != I) not implemented yet!\n");
    PHIST_CHK_IERR(*iflag = PHIST_NOT_IMPLEMENTED, *iflag);
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
  PHIST_CHK_IERR(phist_map_get_comm(A_op->domain_map, &domain_comm, iflag), *iflag);
  const_comm_ptr_t range_comm;
  PHIST_CHK_IERR(phist_map_get_comm(A_op->range_map,  &range_comm,  iflag), *iflag);

  // create mvecs and sdMats
  mvec_ptr_t  V_      = NULL;    //< space for V (current search space)
  mvec_ptr_t  Vtmp_   = NULL;    //< temporary space for V, only used for checking invariants
  mvec_ptr_t  W_     = NULL;    //< space for the orthogonal basis of AV
  mvec_ptr_t  BV_     = NULL;    //< space for BV
  mvec_ptr_t  BW_     = NULL;    //< space for BW
  mvec_ptr_t  Q_      = NULL;    //< Q, enlarged dynamically (converged eigenspace)
  mvec_ptr_t  BQ_     = NULL;    //< B*Q, enlarged dynamically
  mvec_ptr_t  t_      = NULL;    //< space for t (new correction vector(s))
  mvec_ptr_t  At_     = NULL;    //< space for A*t
  mvec_ptr_t  res     = NULL;    //< residuum A*Q-Q*R

  // For harmonic Ritz values (approximating interior eigenvalues) we have
  // V'V=I, W'W=I, A*V=W*H_A, V'Q=W'Q=0 and the (generalized) Schur decomposition
  // H=W'V, H_A = S_L T_A S_R', H = S_L T S_R'. Since H_A is upper triangular
  // by construction, this simplifies to S_L=I, H=T*S_R' and T_A=H_A*S_R'

  sdMat_ptr_t H_      = NULL; //< H=W'V
  sdMat_ptr_t Htmp_   = NULL; //< temporary space used only for checking invariants
  sdMat_ptr_t H_A_    = NULL; //< identity matrix
  sdMat_ptr_t sdMI_   = NULL; //< identity matrix
  // QZ decomposition
  sdMat_ptr_t S_L_    = NULL;
  sdMat_ptr_t S_R_    = NULL;
  sdMat_ptr_t T_      = NULL;
  sdMat_ptr_t T_A_    = NULL;

  _ST_ sigma[nEig_];             //< JaDa correction shifts

  _ST_ *H_raw       = NULL;
  _ST_ *Htmp_raw      = NULL;
  _ST_ *H_A_raw       = NULL;
  _ST_ *S_L_raw       = NULL;
  _ST_ *S_R_raw       = NULL;
  _ST_ *T_raw       = NULL;
  _ST_ *T_A_raw       = NULL;
  lidx_t ldaH, ldaH_A, ldaHtmp, ldaS_L, ldaS_R, ldaT, ldaT_A;

  PHIST_CHK_IERR(SUBR( mvec_create  ) (&V_,     A_op->domain_map, maxBase,        iflag), *iflag);
  // TODO: remove Vtmp
#ifdef TESTING
  PHIST_CHK_IERR(SUBR( mvec_create  ) (&Vtmp_,  A_op->domain_map, maxBase,                iflag), *iflag);
#endif
  PHIST_CHK_IERR(SUBR( mvec_create  ) (&W_,    A_op->range_map,  maxBase,                iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_create  ) (&t_,     A_op->domain_map, blockDim,               iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_create  ) (&At_,    A_op->range_map,  blockDim,               iflag), *iflag);
  int resDim = std::min(2*blockDim, nEig+blockDim-1);
  PHIST_CHK_IERR(SUBR( mvec_create  ) (&res,    A_op->range_map,  resDim,                 iflag), *iflag);

  PHIST_CHK_IERR(SUBR( sdMat_create ) (&H_,     maxBase,          maxBase,  range_comm,   iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_create ) (&H_A_,     maxBase,          maxBase,  range_comm,   iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_create ) (&Htmp_,  maxBase,          maxBase,  range_comm,   iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_create ) (&sdMI_,  maxBase,          maxBase,  range_comm,   iflag), *iflag);

  PHIST_CHK_IERR(SUBR( sdMat_create ) (&S_L_,   maxBase,          maxBase,  range_comm,   iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_create ) (&S_R_,   maxBase,          maxBase,  range_comm,   iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_create ) (&T_,   maxBase,          maxBase,  range_comm,   iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_create ) (&T_A_,   maxBase,          maxBase,  range_comm,   iflag), *iflag);
  // construct identity matrix
  PHIST_CHK_IERR(SUBR( sdMat_identity ) (sdMI_, iflag), *iflag);

  PHIST_CHK_IERR(SUBR( sdMat_extract_view ) (Htmp_,   &Htmp_raw,  &ldaHtmp,  iflag), *iflag);

  PHIST_CHK_IERR(SUBR( sdMat_extract_view ) (H_,   &Htmp_raw,  &ldaH,  iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_extract_view ) (H_A_,   &H_A_raw,  &ldaHtmp,  iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_extract_view ) (S_L_,    &T_raw,   &ldaT,   iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_extract_view ) (S_L_,    &T_A_raw,   &ldaT_A,   iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_extract_view ) (S_L_,    &S_L_raw,   &ldaS_L,   iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_extract_view ) (S_R_,    &S_R_raw,   &ldaS_R,   iflag), *iflag);

  if( B_op != NULL )
  {
    PHIST_CHK_IERR(SUBR( mvec_create )(&BV_,    B_op->range_map,  maxBase,                iflag), *iflag);
    PHIST_CHK_IERR(SUBR( mvec_create )(&BW_,    B_op->range_map,  maxBase,                iflag), *iflag);
  }
  else
  {
    BV_ = V_;
    BW_ = W_;
  }
  // array for the (possibly complex) eigenvalues for SchurDecomp
  CT* ev_H = new CT[maxBase];

  // create views on mvecs and sdMats with current dimensions
  int nV  = minBase;          //< current subspace dimension

  mvec_ptr_t  V   = NULL;     //< B-orthogonal basis of the search space
  mvec_ptr_t  Vful= NULL;     //< B-orthogonal basis of the search space + already locked Schur vectors Q
  mvec_ptr_t  Vtmp= NULL;     //< temporary V
  mvec_ptr_t  Vv  = NULL;     //< next columns in V_
  mvec_ptr_t  W  = NULL;     //< B-orthogonal basis of A*V
  mvec_ptr_t  Wful= NULL;     //< W plus orth basis of A*Q
  mvec_ptr_t  Wv = NULL;     //< next columns in W_
  mvec_ptr_t  BV  = NULL;     //< B*V
  mvec_ptr_t  BVful  = NULL;     //< B*Vful
  mvec_ptr_t  BVv = NULL;     //< next columns in BV_
  mvec_ptr_t  t   = NULL;     //< Block-Jacobi-Davidson correction
  mvec_ptr_t  t_res = NULL;   //< part of the residual AQ-QR corresponding to current block t
  mvec_ptr_t  Qtil= NULL;     //< view of part of Q required for the JaDa correction equation
  mvec_ptr_t BQtil= NULL;     //< B*Qtil
  mvec_ptr_t Q = NULL, BQ = NULL, R = NULL;

  sdMat_ptr_t H   = NULL;     //< projection of A onto V, H=W'*V
  sdMat_ptr_t Hful= NULL;     //< Hful=Wful'*Vful
  sdMat_ptr_t H_Aful= NULL;     //< A*Vful=Wful*H_Aful
  sdMat_ptr_t Htmp= NULL;     //< temporary space for checking invariants

  sdMat_ptr_t S_L = NULL;     //< current left schur vectors of (H,H_A)
  sdMat_ptr_t S_R = NULL;     //< current right schur vectors of (H,H_A)
  sdMat_ptr_t T = NULL;       //< current Schur matrix
  sdMat_ptr_t T_A = NULL;     //< upper triangular matrix of current generalized Schur form

//TODO - adjust sdMats for harmonic extraction

  sdMat_ptr_t Hh  = NULL;     //< inside view for H TODO which views do we need
  sdMat_ptr_t HVv = NULL;     //< next rows in H_
  sdMat_ptr_t HvV = NULL;     //< next columns in H_
  sdMat_ptr_t Hvv = NULL;     //< next block on the diagonal of H_

  sdMat_ptr_t Qq_H = NULL;    //TODO
  mvec_ptr_t  Qq  = NULL;     //TODO
  mvec_ptr_t BQq = NULL;      //TODO
  sdMat_ptr_t Rr_H = NULL;
  sdMat_ptr_t sdMI = NULL;


  //------------------------------- initialize correction equation solver solver ------------------------
  TYPE(jadaCorrectionSolver_ptr) innerSolv = NULL;
  linSolv_t method = opts.innerSolvType;
  PHIST_CHK_IERR(SUBR(jadaCorrectionSolver_create)(&innerSolv, opts, A_op->domain_map, iflag), *iflag);
  std::vector<_MT_> innerTol(nEig_,0.1);
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
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (V_,      &V,                       0,     nV+bs-1,   iflag), *iflag);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (AV_,     &AV,                      0,     nV+bs-1,   iflag), *iflag);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BV_,     &BV,                      0,     nV+bs-1,   iflag), *iflag);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,      &H,     0,      nV+bs-1,  0,     nV-1,      iflag), *iflag);

    PHIST_CHK_IERR(SUBR( simple_blockArnoldi ) (A_op, B_op, V, AV, BV, H, nV, innerBlockDim, iflag), *iflag);
  }
  else
*/
  {
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (V_,      &V,                       0,     nV,        iflag), *iflag);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (W_,     &W,                      0,     nV,        iflag), *iflag);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BV_,     &BV,                      0,     nV,        iflag), *iflag);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,      &H,     0,      minBase,  0,     nV-1,      iflag), *iflag);

    PHIST_CHK_IERR(SUBR( simple_arnoldi ) (A_op, B_op, v0, V, AV, BV, H, nV, iflag), *iflag);
  }

  // set views
  int nConvEig = 0;
#define UPDATE_SUBSPACE_VIEWS \
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (V_,      &V,                         nConvEig, nV-1,      iflag), *iflag); \
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (W_,     &W,                        nConvEig, nV-1,      iflag), *iflag); \
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BV_,     &BV,                        nConvEig, nV-1,      iflag), *iflag); \
  PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,      &H,     nConvEig, nV-1,     nConvEig, nV-1,      iflag), *iflag); \
  PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_A_,      &H_A,     nConvEig, nV-1,     nConvEig, nV-1,      iflag), *iflag); \
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (V_,      &Vful,                      0,        nV-1,      iflag), *iflag); \
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (W_,     &Wful,                     0,        nV-1,      iflag), *iflag); \
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BV_,     &BVful,                     0,        nV-1,      iflag), *iflag); \
  PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,      &Hful,  0,        nV-1,     0,        nV-1,      iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_A_,      &H_Aful,  0,        nV-1,     0,        nV-1,      iflag), *iflag);

  UPDATE_SUBSPACE_VIEWS;



#ifdef TESTING
// V'V=I, W'W=I, A*V=W*H_A, V'Q=W'Q=0 
// H=W'V, H_A S_R = S_L T_A, H S_R = S_L T
#define TESTING_CHECK_SUBSPACE_INVARIANTS \
{ \
  /* check orthogonality of V, BV, Q */ \
  PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Htmp_,&Htmp,0,    nV-1,      0,     nV-1,      iflag), *iflag); \
  PHIST_CHK_IERR(SUBR( mvecT_times_mvec ) (st::one(), Vful, BVful, st::zero(), Htmp, iflag), *iflag); \
  _MT_ orthEps = st::abs(Htmp_raw[0] - st::one()); \
  for(int i = 0; i < nV; i++) \
    for(int j = 0; j < nV; j++) \
      orthEps = mt::max(orthEps, st::abs(Htmp_raw[i*ldaHtmp+j] - ((i==j) ? st::one() : st::zero()))); \
  PHIST_OUT(PHIST_INFO, "Line %d: B-orthogonality of V: %e\n", __LINE__, orthEps); \
 \
  /* check orthogonality of W, BW, Q */ \
  PHIST_CHK_IERR(SUBR( mvecT_times_mvec ) (st::one(), Wful, BWful, st::zero(), Htmp, iflag), *iflag); \
  _MT_ orthEps = st::abs(Htmp_raw[0] - st::one()); \
  for(int i = 0; i < nV; i++) \
    for(int j = 0; j < nV; j++) \
      orthEps = mt::max(orthEps, st::abs(Htmp_raw[i*ldaHtmp+j] - ((i==j) ? st::one() : st::zero()))); \
  PHIST_OUT(PHIST_INFO, "Line %d: B-orthogonality of W: %e\n", __LINE__, orthEps); \
 \
  /* check A*V = W*H_A */ \
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (Vtmp_, &Vtmp,                  0,       nV-1,          iflag), *iflag); \
  PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(st::one(),Wful,H_Aful,st::zero(),iflag),*iflag);
  PHIST_CHK_IERR( A_op->apply(st::one(), A_op->A, Vful, -st::one(), Vtmp, iflag), *iflag); \
  _MT_ normVtmp[nV]; \
  PHIST_CHK_IERR(SUBR( mvec_norm2 ) (Vtmp, normVtmp, iflag), *iflag); \
  _MT_ equalEps = normVtmp[0]; \
  for(int i = 0; i < nV; i++) \
    equalEps = mt::max(equalEps, normVtmp[i]); \
  PHIST_OUT(PHIST_INFO, "Line %d: A*V - W*H_A: %e\n", __LINE__, equalEps); \
  /* check H = W'*V */ \
  PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Htmp_,&Htmp,0,    nV-1,      0,     nV-1,      iflag), *iflag); \
  PHIST_CHK_IERR(SUBR( mvecT_times_mvec ) (st::one(), Wful, Vful, st::zero(), Htmp, iflag), *iflag); \
  PHIST_CHK_IERR(SUBR( sdMat_add_sdMat ) (-st::one(), Hful, st::one(), Htmp, iflag), *iflag); \
  equalEps = st::abs(Htmp_raw[0]); \
  for(int i = 0; i < nV; i++) \
    for(int j = 0; j < nV; j++) \
      equalEps = mt::max(equalEps, st::abs(Htmp_raw[i*ldaHtmp+j])); \
  PHIST_OUT(PHIST_INFO, "Line %d: H - W'*V: %e\n", __LINE__, equalEps); \
  /* PHIST_SOUT(PHIST_INFO, "H:\n"); */ \
  /* PHIST_CHK_IERR(SUBR( sdMat_print )(H, iflag), *iflag); */ \
  /* PHIST_SOUT(PHIST_INFO, "H - W'*V:\n"); */ \
  /* PHIST_CHK_IERR(SUBR( sdMat_print )(Htmp, iflag), *iflag); */ \
}
#else
#define TESTING_CHECK_SUBSPACE_INVARIANTS
#endif


TESTING_CHECK_SUBSPACE_INVARIANTS;


  //----------------------------------- MAIN LOOP ----------------------------------
  int Qsize = 0;
  *nConv=0;
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
  PHIST_CHK_IERR(SUBR(mvec_create)(&newQ, A_op->range_map, newQsize, iflag), *iflag);
  if( Qsize > 0 )
  {
    PHIST_CHK_IERR(SUBR(mvec_set_block)(newQ, Q_, 0, Qsize-1, iflag), *iflag);
  }
  PHIST_CHK_IERR(SUBR(mvec_delete)(Q_, iflag), *iflag);
  Q_ = newQ;

  mvec_ptr_t newBQ = NULL;
  if( B_op != NULL )
  {
    PHIST_CHK_IERR(SUBR(mvec_create)(&newBQ, A_op->range_map, newQsize, iflag), *iflag);
    if( Qsize > 0 )
    {
      PHIST_CHK_IERR(SUBR(mvec_set_block)(newBQ, BQ_, 0, Qsize-1, iflag), *iflag);
    }
    PHIST_CHK_IERR(SUBR(mvec_delete)(Q_, iflag), *iflag);
    BQ_ = newBQ;
  }
  else
  {
    BQ_ = Q_;
  }
  Qsize = newQsize;
}
// update views
PHIST_CHK_IERR(SUBR( mvec_view_block ) (Q_,  &Q,  0, nEig_-1, iflag), *iflag);
PHIST_CHK_IERR(SUBR( mvec_view_block ) (BQ_, &BQ, 0, nEig_-1, iflag), *iflag);
PHIST_CHK_IERR(SUBR( sdMat_view_block ) (R_,  &R, 0, nEig_-1, 0, nEig_-1, iflag), *iflag);

    // for convenience check symmetrie of H (shouldn't much hurt the performance)
    if( symmetric )
    {
      PHIST_CHK_NEG_IERR(SUBR(sdMat_check_symmetrie)(Hful, tol, iflag), *iflag);
    }

    // calculate sorted (generalized) Schur form of (H,H_A) in (S_L,S_R,T,T_A)
    // to get (H,H_A) = (S_L T S_R', S_L T_A S_R')
    //
    // Since we have Vful = [Q V], Wful*H_Aful=A*Vful
    //... TODO ...
    //
    
    // we only update part of S/T, so first set them to zero
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (S_L_,&S_L, 0,     nV-1,      0,     nV-1,      iflag), *iflag);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (S_R_,&S_R, 0,     nV-1,      0,     nV-1,      iflag), *iflag);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (T_,&T, 0,     nV-1,      0,     nV-1,      iflag), *iflag);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (T_A_,&T_A, 0,     nV-1,      0,     nV-1,      iflag), *iflag);
    PHIST_CHK_IERR(SUBR( sdMat_put_value  ) (S_L, st::zero(), iflag), *iflag);
    PHIST_CHK_IERR(SUBR( sdMat_put_value  ) (S_R, st::zero(), iflag), *iflag);
    PHIST_CHK_IERR(SUBR( sdMat_put_value  ) (T, st::zero(), iflag), *iflag);
    PHIST_CHK_IERR(SUBR( sdMat_put_value  ) (T_A, st::zero(), iflag), *iflag);
    // then copy the new block of H, H_A
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (T_,&T, nConvEig, nV-1, nConvEig, nV-1, iflag), *iflag);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (T_A_,&T_A, nConvEig, nV-1, nConvEig, nV-1, iflag), *iflag);
    PHIST_CHK_IERR(SUBR( sdMat_get_block  ) (H_,   T,  nConvEig, nV-1, nConvEig, nV-1, iflag), *iflag);
    PHIST_CHK_IERR(SUBR( sdMat_get_block  ) (H_A_,   T_A,  nConvEig, nV-1, nConvEig, nV-1, iflag), *iflag);
    int nSort = minBase-nConvEig; //nEig_-nConvEig;
    int nSelect = nSort;
    lidx_t offR_H = ldaR_H*nConvEig+nConvEig;
    lidx_t offQ_H = ldaQ_H*nConvEig+nConvEig;
    PHIST_CHK_IERR(SUBR( GenSchurDecomp ) (T_raw+offT, ldaT, T_A_raw+offT_A, ldaT_A, 
                                           S_L_raw+offS, ldaS_L, S_R_raw+offS, ldaS_R,
                                           nV-nConvEig, nSelect, nSort, which, tol, ev_H+nConvEig, iflag), *iflag);
    // we still need to add the missing parts of R_H, Q_H
    if( nConvEig > 0 )
    {
      // upper left part of Q_H
      for(int i = 0; i < nConvEig; i++)
        Q_H_raw[ldaQ_H*i+i] = st::one();

      // left part of R_H
      PHIST_CHK_IERR(SUBR( sdMat_view_block ) (R_H_, &R_H, 0, nConvEig-1, 0, nConvEig-1, iflag), *iflag);
      PHIST_CHK_IERR(SUBR( sdMat_get_block  ) (R,    R_H,  0, nConvEig-1, 0, nConvEig-1, iflag), *iflag);

      // upper right part of R_H
      PHIST_CHK_IERR(SUBR( sdMat_view_block  ) (Q_H_, &Qq_H, nConvEig, nV-1,            nConvEig, nEig_-1, iflag), *iflag);
      PHIST_CHK_IERR(SUBR( sdMat_view_block  ) (R_H_, &Rr_H, 0,             nConvEig-1, nConvEig, nEig_-1, iflag), *iflag);
      PHIST_CHK_IERR(SUBR( sdMat_view_block  ) (H_  , &Hh,   0,             nConvEig-1, nConvEig, nV-1,    iflag), *iflag);

      PHIST_CHK_IERR(SUBR( sdMat_times_sdMat ) (st::one(), Hh, Qq_H, st::zero(), Rr_H, iflag), *iflag);
    }
#ifdef TESTING
{
  // check that H Q_H = Q_H R_H
  PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Q_H_,&Q_H, 0, nV-1,    0, nEig_-1, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_view_block ) (R_H_,&R_H, 0, nEig_-1, 0, nEig_-1, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(sdMat_view_block)(Htmp_, &Htmp, 0, nV-1,    0, nEig_-1, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(sdMat_times_sdMat)(st::one(), Hful, Q_H, st::zero(), Htmp, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(sdMat_times_sdMat)(-st::one(), Q_H, R_H, st::one(), Htmp, iflag), *iflag);
  // get max norm
  _MT_ absErr = mt::zero();
  for(int i = 0; i < nEig_; i++)
    for(int j = 0; j < nV; j++)
      absErr = mt::max(absErr, st::abs(Htmp_raw[i*ldaHtmp+j]));
  PHIST_SOUT(PHIST_INFO, "H*Q_H - Q_H*R_H: %8.4e\n", absErr);
  //PHIST_CHK_IERR(SUBR(sdMat_print)(H, iflag), *iflag);
  //PHIST_CHK_IERR(SUBR(sdMat_print)(Q_H, iflag), *iflag);
  //PHIST_CHK_IERR(SUBR(sdMat_print)(R_H, iflag), *iflag);
  //PHIST_CHK_IERR(SUBR(sdMat_print)(Htmp, iflag), *iflag);
}
#endif

    // update views
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Q_H_,&Q_H, 0,     nV-1,      0,             nV-1,      iflag), *iflag);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (R_H_,&R_H, 0,     nV-1,      0,             nV-1,      iflag), *iflag);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Q_H_,&Qq_H,nConvEig,nV-1,      nConvEig, nEig_-1,   iflag), *iflag);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (R_H_,&Rr_H,0,     nEig_-1,   nConvEig, nEig_-1,   iflag), *iflag);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (Q_,   &Qq,                   nConvEig, nEig_-1,   iflag), *iflag);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BQ_,  &BQq,                  nConvEig, nEig_-1,   iflag), *iflag);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (res, &t_res,                 0, nEig_-nConvEig-1,   iflag), *iflag);

    // update approximate Schur form of A (keeping the already computed part locked)
    PHIST_CHK_IERR(SUBR( mvec_times_sdMat ) (st::one(), V,    Qq_H, st::zero(), Qq,  iflag), *iflag);
    PHIST_CHK_IERR(SUBR( mvec_set_block ) (Q, Qq, nConvEig, nEig_-1, iflag), *iflag);
    PHIST_CHK_IERR(SUBR( sdMat_set_block  ) (R, Rr_H, 0, nEig_-1, nConvEig, nEig_-1, iflag), *iflag);
    PHIST_CHK_IERR(SUBR( mvec_times_sdMat ) (st::one(), AV,   Qq_H, st::zero(), t_res, iflag), *iflag);
    if( B_op != NULL )
    {
      PHIST_CHK_IERR(SUBR( mvec_times_sdMat ) (st::one(), BV,   Qq_H, st::zero(), BQq,  iflag), *iflag);
      PHIST_CHK_IERR(SUBR( mvec_set_block ) (BQ, BQq, nConvEig, nEig_-1, iflag), *iflag);
    }
    // overwrite res with the residual: -res = -(Aq - Bqr) = + BQq*r - res
    PHIST_CHK_IERR(SUBR( mvec_times_sdMat ) (st::one(), BQ,   Rr_H,    -st::one(), t_res, iflag), *iflag);
    // calculate norm of the residual
    PHIST_CHK_IERR(SUBR( mvec_norm2 ) (t_res, resNorm+nConvEig, iflag), *iflag);
#ifdef TESTING
{
  // check that the residual is orthogonal to Q (should be by construction!)
  PHIST_CHK_IERR( SUBR( sdMat_view_block ) (Htmp_, &Htmp, 0, nEig_-1, nConvEig, nEig_-1, iflag), *iflag);
  PHIST_CHK_IERR( SUBR( mvecT_times_mvec ) (st::one(), Q, t_res, st::zero(), Htmp, iflag), *iflag);
  _MT_ equalEps = st::abs(Htmp_raw[0]);
  for(int i = 0; i < nEig_; i++)
    for(int j = nConvEig; j < nEig_; j++)
      equalEps = mt::max(equalEps, st::abs(Htmp_raw[i*ldaHtmp+j]));
  PHIST_OUT(PHIST_INFO, "Res orthogonality wrt. Q: %e\n", equalEps);
  // calculate explicit residual
  // AQ - BQR
  PHIST_CHK_IERR( SUBR( mvec_view_block ) (Vtmp_, &Vtmp, 0, nEig_-1, iflag), *iflag);
  PHIST_CHK_IERR( A_op->apply(st::one(), A_op->A, Q, st::zero(), Vtmp, iflag), *iflag);
  PHIST_CHK_IERR( SUBR( mvec_times_sdMat ) (-st::one(), BQ, R, st::one(), Vtmp, iflag), *iflag);
  _MT_ expRes[nEig_];
  PHIST_CHK_IERR( SUBR( mvec_norm2 )(Vtmp, expRes, iflag), *iflag);
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
      PHIST_CHK_IERR( SUBR(ReorderPartialSchurDecomp)(R_H_raw+offR_H, ldaR_H, Q_H_raw+offQ_H, ldaQ_H, nV-nConvEig, nEig_-nConvEig, which, sqrt(tol), resNorm+nConvEig, ev_H+nConvEig, &resPermutation[nConvEig], iflag), *iflag);
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
        //PHIST_CHK_IERR(SUBR(sdMat_view_block)(Htmp_, &Htmp, i, nEig_-1, i, nEig_-1, iflag), *iflag);
        //PHIST_CHK_IERR(SUBR(sdMat_put_value)(Htmp, st::zero(), iflag), *iflag);
        //for(int j = i; j < nEig_; j++)
        //{
          //int j_ = resPermutation[j];
          //Htmp_raw[j_+j*ldaHtmp] = st::one();
        //}
        //PHIST_CHK_IERR(SUBR(mvec_view_block)(Q_, &Qq, i, nEig_-1, iflag), *iflag);
        //PHIST_CHK_IERR(SUBR(mvec_times_sdMat_inplace)(Qq, Htmp, iflag), *iflag);
        //if( B_op != NULL )
        //{
          //PHIST_CHK_IERR(SUBR(mvec_view_block)(BQ_, &BQq, i, nEig_-1, iflag), *iflag);
          //PHIST_CHK_IERR(SUBR(mvec_times_sdMat_inplace)(BQq, Htmp, iflag), *iflag);
        //}
        PHIST_CHK_IERR(SUBR( mvec_times_sdMat ) (st::one(), V,    Qq_H, st::zero(), Qq,  iflag), *iflag);
        PHIST_CHK_IERR(SUBR( sdMat_set_block  ) (R, Rr_H, 0, nEig_-1, nConvEig, nEig_-1, iflag), *iflag);
        break;
      }
    }
#ifdef TESTING
{
  // check that Q is in correct order
  PHIST_CHK_IERR( SUBR(mvec_view_block)(Vtmp_, &t, 0, nEig_-1, iflag), *iflag);
  PHIST_CHK_IERR(A_op->apply(st::one(), A_op->A, Q, st::zero(), t, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(-st::one(), BQ, R, st::one(), t, iflag), *iflag);
  _MT_ reorderedResNorm[nEig_];
  PHIST_CHK_IERR(SUBR(mvec_norm2)(t, reorderedResNorm, iflag), *iflag);
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
    for(int i = nConvEig; i < std::min(nEig,nEig_); i++)
    {
      if (resNorm[i]>0)
      {
        PHIST_SOUT(PHIST_INFO,"In iteration %d: Current approximation for eigenvalue %d is %16.8g%+16.8gi with residual %e\n", 
        *nIter, i+1, ct::real(ev_H[i]),ct::imag(ev_H[i]), resNorm[i]);
      }
      else
      {
        PHIST_SOUT(PHIST_INFO,"In iteration %d: Current approximation for eigenvalue %d is %16.8g%+16.8gi (residual not yet available)\n", 
        *nIter, i+1, ct::real(ev_H[i]),ct::imag(ev_H[i]));
      }
#ifndef IS_COMPLEX
      if( mt::abs(ct::imag(ev_H[i])) > tol && blockDim == 1 )
      {
        PHIST_SOUT(PHIST_WARNING, "Detected possible complex-conjugate eigenpair, but blockDim == 1 (you need at least blockDim=2 to detect complex-conjugate eigenpairs correctly)!\n");
      }
#endif
      if( resNorm[i] <= tol && i == nConvEig+nNewConvEig )
      {
#ifndef IS_COMPLEX
        // detect complex conjugate eigenvalue pairs
        if( mt::abs(ct::imag(ev_H[i])) > tol && i+1 < nEig_ )
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
      *nConv=nConvEig;
      PHIST_SOUT(PHIST_INFO,"In iteration %d: all eigenvalues converged!\n", *nIter);
      break;
    }

    if( *nIter >= maxIter )
    {
      *nConv=nConvEig;
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
          PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Q_H_,  &Q_H,  nConvEig, nV-1,          nConvEig, minBase-1,    iflag), *iflag);
        }
        else
        {
          PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Q_H_,  &Q_H,  0, nV-1,          0, minBase-1,    iflag), *iflag);
        }
      }
      else
      {
        if( symmetric )
        {
          PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Q_H_,  &Q_H,  nConvEig, nV-1,          nConvEig, nV-1,    iflag), *iflag);
        }
        else
        {
          PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Q_H_,  &Q_H,  0, nV-1,          0, nV-1,    iflag), *iflag);
        }
      }

      // reorder V and H
      if( symmetric )
      {
        PHIST_CHK_IERR(SUBR( transform_searchSpace ) (V, AV, BV, H, Q_H, B_op != NULL, iflag), *iflag);
      }
      else
      {
        PHIST_CHK_IERR(SUBR( transform_searchSpace ) (Vful, AVful, BVful, Hful, Q_H, B_op != NULL, iflag), *iflag);
      }

      nConvEig = nConvEig+nNewConvEig;
      *nConv=nConvEig;

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


    // shrink search space if necessary
    if( nV + k > maxBase )
    {
      PHIST_SOUT(PHIST_INFO,"Shrinking search space from %d to %d\n", nV, minBase);

      // nothing to do if converged eigenvalues this iteration
      if( nNewConvEig == 0 )
      {
        if( symmetric )
        {
          PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Q_H_,  &Q_H,  nConvEig, nV-1,          nConvEig, minBase-1,    iflag), *iflag);

          PHIST_CHK_IERR(SUBR( transform_searchSpace ) (V, AV, BV, H, Q_H, B_op != NULL, iflag), *iflag);
        }
        else
        {
          PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Q_H_,  &Q_H,  0, nV-1,          0, minBase-1,    iflag), *iflag);

          PHIST_CHK_IERR(SUBR( transform_searchSpace ) (Vful, AVful, BVful, Hful, Q_H, B_op != NULL, iflag), *iflag);
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
      PHIST_CHK_IERR(SUBR( mvec_view_block  ) (t_,  &t,     0, k-1,  iflag), *iflag);
      // we only need to view first part of Q
      PHIST_CHK_IERR(SUBR( mvec_view_block  ) (Q_,   &Qtil,  0, k_, iflag), *iflag);
      PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BQ_,  &BQtil, 0, k_, iflag), *iflag);
      // set tolerances
      for(int i = 0; i < k; i++)
      {
        if( resNorm[nConvEig+i] > 4*lastOuterRes[nConvEig+i] )
          innerTol[nConvEig+i] = 0.1;
        if( innerTol[nConvEig+i] > mt::eps() )
          innerTol[nConvEig+i] = mt::max(mt::eps(), innerTol[nConvEig+i]*0.1);
        innerTol[nConvEig+i] = mt::max(innerTol[nConvEig+i], 0.1*tol/(mt::eps()+resNorm[nConvEig+i]));
        innerTol[nConvEig+i] = mt::min(innerTol[nConvEig+i], 0.1);
        lastOuterRes[nConvEig+i] = resNorm[nConvEig+i];
      }
PHIST_SOUT(PHIST_INFO,"selectedRes: ");
for(int i = 0; i < k; i++)
  PHIST_SOUT(PHIST_INFO,"\t%d (%e, tol %e)", selectedRes[i], resNorm[selectedRes[i]], innerTol[nConvEig+i]);
PHIST_SOUT(PHIST_INFO,"\n");


      for(int i = 0; i < blockDim; i++)
        selectedRes[i] -= nConvEig-nNewConvEig;
      PHIST_CHK_NEG_IERR(SUBR(jadaCorrectionSolver_run)(innerSolv, A_op, B_op, Qtil, BQtil, sigma, res, &selectedRes[0],
                                                    &innerTol[nConvEig], innerMaxBase, t, innerIMGS, innerGMRESabortAfterFirstConverged, iflag), *iflag);

      // get solution and reuse res for At
      PHIST_CHK_IERR(SUBR( mvec_view_block  ) (t_, &Vv,  0, k-1, iflag), *iflag);
      PHIST_CHK_IERR(SUBR( mvec_view_block  ) (At_,&AVv, 0, k-1, iflag), *iflag);

      // enlarge search space
      // first update views
      PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,  &HVv, 0,     nV-1,      nV,    nV+k-1,    iflag), *iflag);
      PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,  &HvV, nV,    nV+k-1,    0,     nV-1,      iflag), *iflag);
      PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,  &Hvv, nV,    nV+k-1,    nV,    nV+k-1,    iflag), *iflag);
      // orthogonalize t as Vv (reuse R_H)
      PHIST_CHK_IERR(SUBR( sdMat_view_block ) (R_H_,&R_H, 0,     nV-1,      0,     k-1,       iflag), *iflag);
      PHIST_CHK_IERR(SUBR( sdMat_view_block ) (R_H_,&Rr_H,nV,    nV+k-1,    nV,    nV+k-1,    iflag), *iflag);
      int rankV;
      PHIST_CHK_NEG_IERR(SUBR( orthog ) (Vful, Vv, NULL, Rr_H, R_H, 5, &rankV, iflag), *iflag);
      // TODO: only take non-random vector if *iflag > 0
      // calculate AVv, BVv
      PHIST_CHK_IERR( A_op->apply(st::one(), A_op->A, Vv, st::zero(), AVv, iflag), *iflag);
      if( B_op != NULL )
      {
        PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BV_, &BVv,                   nV,    nV+k-1,    iflag), *iflag);
        PHIST_CHK_IERR( B_op->apply(st::one(), B_op->A, Vv, st::zero(), BVv, iflag), *iflag);
      }
      // update H
      PHIST_CHK_IERR(SUBR( mvecT_times_mvec ) (st::one(), Vful,  AVv, st::zero(), HVv, iflag), *iflag);
      // for the symmetric case use AVv*V here, so we don't need AV at all
      if( !symmetric )
      {
        PHIST_CHK_IERR(SUBR( mvecT_times_mvec ) (st::one(), Vv, AVful,  st::zero(), HvV, iflag), *iflag);
      }
      else
      {
        PHIST_CHK_IERR(SUBR(sdMat_view_block) (sdMI_, &sdMI, 0, nV-1, 0, nV-1, iflag), *iflag);
        PHIST_CHK_IERR(SUBR(sdMatT_times_sdMat)(st::one(), HVv, sdMI, st::zero(), HvV, iflag), *iflag);
      }
      PHIST_CHK_IERR(SUBR( mvecT_times_mvec ) (st::one(), Vv, AVv, st::zero(), Hvv, iflag), *iflag);
      // use set block to put Vv and AVv really into V and AV
      PHIST_CHK_IERR(SUBR( mvec_set_block ) (V_,  Vv,  nV, nV+k-1, iflag), *iflag);
      PHIST_CHK_IERR(SUBR( mvec_set_block ) (AV_, AVv, nV, nV+k-1, iflag), *iflag);
      // increase nV
      nV = nV + k;
    } // k > 0

    UPDATE_SUBSPACE_VIEWS;
TESTING_CHECK_SUBSPACE_INVARIANTS;
  }

  // copy result to Q_
  PHIST_CHK_IERR(SUBR(mvec_set_block)(Q__, Q, 0, nEig_-1, iflag), *iflag);
  // copy resulting eigenvalues to ev
  for(int i = 0; i < nEig_; i++)
    ev[i] = ev_H[i];


  //------------------------------- delete vectors and matrices --------------------
  PHIST_CHK_IERR(SUBR(jadaCorrectionSolver_delete)(innerSolv, iflag), *iflag);

  // delete views
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (R_H, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Rr_H,iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Q_H, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Qq_H,iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Hvv, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (HvV, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (HVv, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Hh,  iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (H,   iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Hful,iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Htmp,iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (sdMI,iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (R,   iflag), *iflag);

  PHIST_CHK_IERR(SUBR( mvec_delete  ) (Qq,  iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (BQq, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (t,   iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (t_res,iflag),*iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (BVv, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (BV,  iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (AVv, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (AV,  iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (Vv,  iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (Vtmp,iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (V,   iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (Vful,iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (AVful,iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (BVful,iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (Qtil,iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (BQtil,iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (Q,   iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (BQ,  iflag), *iflag);

  delete[] ev_H;

  // delete mvecs and sdMats
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Q_H_,iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (R_H_,iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Htmp_,iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (H_,  iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (sdMI_,iflag), *iflag);
  if( B_op != NULL )
  {
    PHIST_CHK_IERR(SUBR( mvec_delete )(BV_, iflag), *iflag);
    PHIST_CHK_IERR(SUBR( mvec_delete )(BQ_, iflag), *iflag);
  }
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (Q_,  iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (t_,  iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (At_, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (res, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (AV_, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (Vtmp_,iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (V_,  iflag), *iflag);
}

