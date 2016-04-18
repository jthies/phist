namespace {
//! small helper function that checks if an sdMat is symmetric
void SUBR(sdMat_check_symmetry)(TYPE(const_sdMat_ptr) mat, _MT_ tol, int*iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"

  // check dimensions
  int m, n;
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(mat, &m, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(mat, &n, iflag), *iflag);
  if( m != n )
  {
    *iflag = 1;
    return;
  }

  // create tmp storage
  TYPE(sdMat_ptr) tmp = NULL;
  PHIST_CHK_IERR(SUBR(sdMat_create)(&tmp, m, n, NULL, iflag), *iflag);

  // construct identity matrix
  TYPE(sdMat_ptr) id = NULL;
  PHIST_CHK_IERR(SUBR(sdMat_create)(&id, m, n, NULL, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_identity) (id, iflag), *iflag);

  // set tmp to transposed
  PHIST_CHK_IERR(SUBR(sdMatT_times_sdMat)(st::one(), mat, id, st::zero(), tmp, iflag), *iflag);
  // subtract mat
  PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(-st::one(), mat, st::one(), tmp, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(sdMat_from_device)(tmp,iflag),*iflag);
  // calc max. abs. value of mat^T-mat
  _MT_ maxVal = mt::zero();
  {
    _ST_ *tmp_raw = NULL;
    phist_lidx lda;
    PHIST_CHK_IERR(SUBR(sdMat_extract_view)(tmp, &tmp_raw, &lda, iflag), *iflag);
    for(int i = 0; i < m; i++)
      for(int j = 0; j < n; j++)
        maxVal = mt::max(maxVal, st::abs(tmp_raw[i*lda+j]));
  }

  PHIST_CHK_IERR(SUBR(sdMat_delete)(id, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(tmp, iflag), *iflag);

  if( maxVal < tol )
  {
    PHIST_SOUT(PHIST_VERBOSE, "Symmetry deviation of projection %e\n", maxVal);
    *iflag = 0;
  }
  else
  {
    PHIST_SOUT(PHIST_WARNING, "Symmetry deviation of projection %e exceeds tolerance (%e)\n", maxVal, tol);
    *iflag = 1;
  }
}
}

//! Subspace Jacobi-Davidson for interior eigenvalues
//!
//! Tries to compute a partial schur form $(Q,R)$ of dimension opts.numEigs
//! of the stencil $A*x-\lambda*B*x$ with a general linear operator $A$ and a
//! hermitian positive definite (hpd.) linear operator $B$ using a
//! block-Jacobi-Davidson QR method based on harmonic Ritz values. <br>
//! The generalized eigenvalues $\lambda_i$ are the diagonal entries of the
//! partial schur form $A*Q = B*Q*R$ returned. <br>
//!
//! see header file for further documentation of the parameters
//!
void SUBR(harmonicjada)( TYPE(const_linearOp_ptr) A_op,  TYPE(const_linearOp_ptr) B_op,
                         phist_jadaOpts opts,
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
  phist_EeigSort which=opts.which;
  _MT_ tol=opts.convTol;  
  int nEig=opts.numEigs;
  int blockDim=opts.blockSize;
  int maxIter =opts.maxIters;
                         
int minBase=opts.minBas;
int maxBase=opts.maxBas;
int innerBlockDim=opts.innerSolvBlockSize;        
int innerMaxBase=opts.innerSolvMaxBas;

int arno=opts.arno;
int initialShiftIter=opts.initialShiftIters;   
_ST_ initialShift   =(_ST_)opts.initialShift_r;
                    +(_ST_)opts.initialShift_i*st::cmplx_I();
                         
bool innerIMGS=(opts.innerSolvRobust!=0);
bool innerGMRESabortAfterFirstConverged=opts.innerSolvStopAfterFirstConverged;
bool symmetric=opts.symmetry==phist_HERMITIAN;
#ifndef IS_COMPLEX
symmetric=symmetric||(opts.symmetry==phist_COMPLEX_SYMMETRIC);
#endif

  phist_EeigExtr how=opts.how;
  if (how==phist_STANDARD)
  {
    PHIST_SOUT(PHIST_ERROR,"if you want to use standard Ritz values, please use the subspacejada routine instead\n");
  }
  if (how!=phist_HARMONIC)
  {
    PHIST_SOUT(PHIST_ERROR,"only Harmonic Ritz extraction is implemented (jadaOpts.how=%s), found %s\n",
        eigExtr2str(phist_HARMONIC),eigExtr2str(how));
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
  phist_const_comm_ptr domain_comm;
  PHIST_CHK_IERR(phist_map_get_comm(A_op->domain_map, &domain_comm, iflag), *iflag);
  phist_const_comm_ptr range_comm;
  PHIST_CHK_IERR(phist_map_get_comm(A_op->range_map,  &range_comm,  iflag), *iflag);

  // create mvecs and sdMats
  mvec_ptr  V_      = NULL;    //< space for V (current search space)
  mvec_ptr  Vtmp_   = NULL;    //< temporary space for V, only used for checking invariants
  mvec_ptr  W_     = NULL;    //< space for the orthogonal basis of AV
  mvec_ptr  BV_     = NULL;    //< space for BV TODO - B-orthogonalization not fully implemented
  mvec_ptr  Q_      = NULL;    //< Q, enlarged dynamically (converged eigenspace)
  mvec_ptr  BQ_     = NULL;    //< B*Q, enlarged dynamically
  mvec_ptr  t_      = NULL;    //< space for t (new correction vector(s))
  mvec_ptr  At_     = NULL;    //< space for A*t
  mvec_ptr  res     = NULL;    //< residuum A*Q-Q*R

  // For harmonic Ritz values (approximating interior eigenvalues) we have
  // V'V=I, W'W=I, H=W'V, A*V=W*H_A, V'Q=W'Q=0 and the (generalized) Schur decomposition
  // H_A = S_L T_A S_R', H = S_L T S_R'.
  
  sdMat_ptr H_      = NULL; //< H=W'V
  sdMat_ptr Htmp_   = NULL; //< temporary space used only for checking invariants
  sdMat_ptr H_A_    = NULL; //< identity matrix
  sdMat_ptr H_Atmp_ = NULL; //< temporary space used only for checking invariants
  // QZ decomposition
  sdMat_ptr S_L_    = NULL;
  sdMat_ptr S_R_    = NULL;
  sdMat_ptr T_      = NULL;
  sdMat_ptr T_A_    = NULL;  

  _ST_ sigma[nEig_];             //< JaDa correction shifts

  _ST_ *H_raw       = NULL;
  _ST_ *Htmp_raw      = NULL;
  _ST_ *H_A_raw       = NULL;
  _ST_ *H_Atmp_raw    = NULL;
  _ST_ *S_L_raw       = NULL;
  _ST_ *S_R_raw       = NULL;
  _ST_ *T_raw       = NULL;
  _ST_ *T_A_raw       = NULL;
  phist_lidx ldaH, ldaH_A, ldaHtmp, ldaH_Atmp,ldaS_L, ldaS_R, ldaT, ldaT_A;

  PHIST_CHK_IERR(SUBR( mvec_create  ) (&V_,     A_op->domain_map, maxBase,        iflag), *iflag);
#ifdef TESTING
  PHIST_CHK_IERR(SUBR( mvec_create  ) (&Vtmp_,  A_op->domain_map, maxBase,                iflag), *iflag);
#endif
  PHIST_CHK_IERR(SUBR( mvec_create  ) (&W_,    A_op->range_map,  maxBase,                iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_create  ) (&t_,     A_op->domain_map, blockDim,               iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_create  ) (&At_,    A_op->range_map,  blockDim,               iflag), *iflag);
  int resDim = std::min(2*blockDim, nEig+blockDim-1);
  PHIST_CHK_IERR(SUBR( mvec_create  ) (&res,    A_op->range_map,  resDim,                 iflag), *iflag);

  PHIST_CHK_IERR(SUBR( sdMat_create ) (&H_,      maxBase,          maxBase,  range_comm,   iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_create ) (&H_A_,    maxBase,          maxBase,  range_comm,   iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_create ) (&Htmp_,   maxBase,          maxBase,  range_comm,   iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_create ) (&H_Atmp_, maxBase,          maxBase,  range_comm,   iflag), *iflag);

  PHIST_CHK_IERR(SUBR( sdMat_create ) (&S_L_,   maxBase,          maxBase,  range_comm,   iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_create ) (&S_R_,   maxBase,          maxBase,  range_comm,   iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_create ) (&T_,   maxBase,          maxBase,  range_comm,   iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_create ) (&T_A_,   maxBase,          maxBase,  range_comm,   iflag), *iflag);

  PHIST_CHK_IERR(SUBR( sdMat_extract_view ) (Htmp_,   &Htmp_raw,   &ldaHtmp,    iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_extract_view ) (H_Atmp_, &H_Atmp_raw, &ldaH_Atmp,  iflag), *iflag);

  PHIST_CHK_IERR(SUBR( sdMat_extract_view ) (H_,     &H_raw,  &ldaH,  iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_extract_view ) (H_A_,   &H_A_raw,  &ldaHtmp,  iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_extract_view ) (S_L_,    &T_raw,   &ldaT,   iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_extract_view ) (S_L_,    &T_A_raw,   &ldaT_A,   iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_extract_view ) (S_L_,    &S_L_raw,   &ldaS_L,   iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_extract_view ) (S_R_,    &S_R_raw,   &ldaS_R,   iflag), *iflag);

  if( B_op != NULL )
  {
    PHIST_CHK_IERR(SUBR( mvec_create )(&BV_,    B_op->range_map,  maxBase,                iflag), *iflag);
  }
  else
  {
    BV_ = V_;
  }
  // array for the (possibly complex) harmonic Ritz values (retruned by GenSchurDecomp)
  CT* ev_H = new CT[maxBase];
  
  // As discussed in the eigentemplates book (Alg. 7.19, step 16, p. 226/227)
  // we take as shifts the Rayleigh quotients associated with the first few  
  // harmonic Ritz vectors instead of the harmonic Ritz values. We get them  
  // from the requirement (Au-theta*U) orth u, which gives theta_j = conj(T_jj)*T_A_jj
  CT* rq_H = new CT[maxBase];

  // create views on mvecs and sdMats with current dimensions
  int nV;                     //< current subspace dimension

  mvec_ptr  V   = NULL;     //< B-orthogonal basis of the search space
  mvec_ptr  Vful= NULL;     //< B-orthogonal basis of the search space + already locked Schur vectors Q
  mvec_ptr  Vtmp= NULL;     //< temporary V
  mvec_ptr  Vv  = NULL;     //< next columns in V_
  mvec_ptr  W  = NULL;     //< B-orthogonal basis of A*V
  mvec_ptr  Wful= NULL;     //< W plus orth basis of A*Q
  mvec_ptr  Wv = NULL;     //< next columns in W_
  mvec_ptr  BV  = NULL;     //< B*V
  mvec_ptr  BVful  = NULL;     //< B*Vful
  mvec_ptr  BVv = NULL;     //< next columns in BV_
  mvec_ptr  t   = NULL;     //< Block-Jacobi-Davidson correction
  mvec_ptr  t_res = NULL;   //< part of the residual AQ-QR corresponding to current block t
  mvec_ptr  Qtil= NULL;     //< view of part of Q required for the JaDa correction equation
  mvec_ptr BQtil= NULL;     //< B*Qtil
  mvec_ptr Q = NULL, BQ = NULL, R = NULL;

  mvec_ptr  Qq  = NULL;     // view of the 'locked' part of Q, the next few eigenvectors to converge
  mvec_ptr BQq = NULL;      // view of B*Qq

  sdMat_ptr H   = NULL;     //< projection of A onto V, H=W'*V
  sdMat_ptr Hful= NULL;     //< Hful=Wful'*Vful
  sdMat_ptr H_A = NULL;     //< A*V=W*H_A
  sdMat_ptr H_Aful= NULL;     //< A*Vful=Wful*H_Aful

  sdMat_ptr Htmp= NULL;     //< temporary space for checking invariants
  sdMat_ptr H_Atmp= NULL;     //< temporary space for checking invariants

  sdMat_ptr S_L = NULL;     //< current left schur vectors of (H,H_A)
  sdMat_ptr S_R = NULL;     //< current right schur vectors of (H,H_A)
  sdMat_ptr T = NULL;       //< current Schur matrix
  sdMat_ptr T_A = NULL;     //< upper triangular matrix of current generalized Schur form

//views of parts of sdMats
  sdMat_ptr HVv = NULL;     //< next rows in H_A_
  sdMat_ptr HvV = NULL;     //< next columns in H_A_
  sdMat_ptr Hvv = NULL;     //< next block on the diagonal of H_A

  sdMat_ptr Hq  = NULL;     //< inside view for H(iq,jv), with iq=0:nconv-1, jv=nconv:nV-1
  sdMat_ptr H_Aq  = NULL;   //< inside view for H_A(iq,jv)
  sdMat_ptr Tq  = NULL;     //< inside view for T(iq,jq), jq=nconv:nEig_-1
                              //< nEig_<=nEig is used to consider only the next few sought eigenvalues
                              //< at a time.
  sdMat_ptr T_Aq  = NULL;   //< inside view for T_A(iq,jq)

  //------------------------------- initialize correction equation solver solver ------------------------
  TYPE(jadaCorrectionSolver_ptr) innerSolv = NULL;
  phist_ElinSolv method = opts.innerSolvType;
  PHIST_CHK_IERR(SUBR(jadaCorrectionSolver_create)(&innerSolv, opts, A_op->domain_map, iflag), *iflag);
  std::vector<_MT_> innerTol(nEig_,0.1);
  std::vector<_MT_> lastOuterRes(nEig_,mt::zero());

  if (!arno && false)
  {
    // This driver is intended for computing inner eigenvalues near 0, so
    // an Arnoldi process is not the best idea. Instead, we will lock the shifts
    // to 0 until nV=minBase has been reached.
    PHIST_CHK_IERR(SUBR(mvec_put_value)(V_,st::zero(),iflag),*iflag);
    nV=0;
    if (v0!=NULL)
    {
      PHIST_CHK_IERR(SUBR(mvec_num_vectors)(v0,&nV,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(mvec_set_block)(V_,v0,0,nV-1,iflag),*iflag);
    }

    // we need at least min(2*blockDim,nEig+blockDim) vectors because this is the initial size of Q
    // (some unconverged vectors are deflated along with the converged ones in this algorithm).
    // The missing vectors (all of them if v0==NULL) are generated automatically by the orthog
    // routine, which always creates some orthogonal full rank matrix Q.
    nV=std::max(nV, std::min(2*blockDim,nEig+blockDim));

    // orthogonalize the initial vector block
    PHIST_CHK_IERR(SUBR(mvec_view_block)(V_,&Vv,0,nV-1,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(Htmp_,  &Htmp,0,  nV-1,0,nV-1,iflag),*iflag);
    int rank_v0;
    *iflag=PHIST_ROBUST_REDUCTIONS;
    // note: ignore positive iflag, it just means the input was rank deficient and the
    // corresponding columns were randomized.
    PHIST_CHK_NEG_IERR(SUBR(orthog)(NULL,Vv,B_op,Htmp,NULL,3,&rank_v0,iflag),*iflag);

    // W=AV
    PHIST_CHK_IERR(SUBR(mvec_view_block)(V_,&Vv,0,nV-1,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_view_block)(W_,&Wv,0,nV-1,iflag),*iflag);
    PHIST_CHK_IERR(A_op->apply(st::one(),A_op->A,Vv,st::zero(),Wv,iflag),*iflag);
  }
  else
  {
    nV=minBase;
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (V_,      &V,                       0,     nV,        iflag), *iflag);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (W_,     &W,                      0,     nV,        iflag), *iflag);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BV_,     &BV,                      0,     nV,        iflag), *iflag);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,      &H,     0,      minBase,  0,     nV-1,      iflag), *iflag);

    PHIST_CHK_IERR(SUBR( simple_arnoldi ) (A_op, B_op, v0, V, W, BV, H, nV, iflag), *iflag);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (W_, &Wv, 0,     nV-1, iflag), *iflag);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (V_, &Vv, 0,     nV-1, iflag), *iflag);
  }
    
    // we now have V n x nV: V'BV=I, and W=AV

  PHIST_CHK_IERR(SUBR(sdMat_view_block)(H_A_,  &H_A,0,  nV-1,0,nV-1,iflag),*iflag);

  int rank_w0;
  *iflag=PHIST_ROBUST_REDUCTIONS;
  PHIST_CHK_NEG_IERR(SUBR(orthog)(NULL,Wv,B_op,H_A,NULL,3,&rank_w0,iflag),*iflag);
  // initial H=W'V
  PHIST_CHK_IERR(SUBR(sdMat_view_block)(H_,  &H,0,  nV-1,0,nV-1,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),Wv,Vv,st::zero(),H,iflag),*iflag);
  
  initialShiftIter=std::max(minBase-nV,initialShiftIter);
  if (initialShiftIter>0) PHIST_SOUT(PHIST_VERBOSE,"perform %d iterations with fixed shift 0\n",initialShiftIter);

  // set views
  int nConvEig = 0;

#ifdef UPDATE_SUBSPACE_VIEWS
#undef UPDATE_SUBSPACE_VIEWS
#endif
#define UPDATE_SUBSPACE_VIEWS \
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (V_,   &V,      nConvEig, nV-1,                 iflag), *iflag); \
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (W_,   &W,      nConvEig, nV-1,                 iflag), *iflag); \
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BV_,  &BV,     nConvEig, nV-1,                 iflag), *iflag); \
  PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,   &H,      nConvEig, nV-1, nConvEig, nV-1, iflag), *iflag); \
  PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_A_, &H_A,    nConvEig, nV-1, nConvEig, nV-1, iflag), *iflag); \
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (V_,   &Vful,   0,        nV-1,                 iflag), *iflag); \
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (W_,   &Wful,   0,        nV-1,                 iflag), *iflag); \
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BV_,  &BVful,  0,        nV-1,                 iflag), *iflag); \
  PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,   &Hful,   0,        nV-1, 0,        nV-1, iflag), *iflag); \
  PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_A_, &H_Aful, 0,        nV-1, 0,        nV-1, iflag), *iflag);

  UPDATE_SUBSPACE_VIEWS;


#ifdef TESTING_CHECK_SUBSPACE_INVARIANTS
#undef TESTING_CHECK_SUBSPACE_INVARIANTS
#endif
#ifdef TESTING
// V'V=I, W'W=I, A*V=W*H_A, V'Q=W'Q=0 
// H=W'V, H_A S_R = S_L T_A, H S_R = S_L T
#define TESTING_CHECK_SUBSPACE_INVARIANTS \
{ \
  PHIST_CHK_IERR(SUBR( mvec_view_block  ) (Vtmp_, &Vtmp, 0, nV-1,          iflag), *iflag); \
  PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Htmp_, &Htmp, 0, nV-1, 0, nV-1, iflag), *iflag); \
  /* check (B-)orthogonality of V, Q */ \
  PHIST_CHK_IERR(SUBR( mvecT_times_mvec ) (st::one(), Vful, BVful, st::zero(), Htmp, iflag), *iflag); \
  _MT_ orthEps = st::abs(Htmp_raw[0] - st::one()); \
  for(int i = 0; i < nV; i++) \
    for(int j = 0; j < nV; j++) \
      orthEps = mt::max(orthEps, st::abs(Htmp_raw[i*ldaHtmp+j] - ((i==j) ? st::one() : st::zero()))); \
  PHIST_OUT(PHIST_INFO, "Line %d: B-orthogonality of V: %e\n", __LINE__, orthEps); \
 \
  /* check (B-)orthogonality of W, Q */ \
  if (B_op==NULL) {\
    PHIST_CHK_IERR(SUBR( mvecT_times_mvec ) (st::one(), Wful, Wful, st::zero(), Htmp, iflag), *iflag); \
  } else { \
    PHIST_CHK_IERR(B_op->apply(st::one(),B_op->A,Wful,st::zero(),Vtmp,iflag),*iflag); \
    PHIST_CHK_IERR(SUBR(mvecT_times_mvec) (st::one(), Wful, Vtmp, st::zero(), Htmp, iflag), *iflag); \
  } \
  orthEps = st::abs(Htmp_raw[0] - st::one()); \
  for(int i = 0; i < nV; i++) \
    for(int j = 0; j < nV; j++) \
      orthEps = mt::max(orthEps, st::abs(Htmp_raw[i*ldaHtmp+j] - ((i==j) ? st::one() : st::zero()))); \
  PHIST_OUT(PHIST_INFO, "Line %d: B-orthogonality of W: %e\n", __LINE__, orthEps); \
 \
  /* check A*V = W*H_A */ \
  PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(st::one(),Wful,H_Aful,st::zero(),Vtmp,iflag),*iflag); \
  PHIST_CHK_IERR( A_op->apply(st::one(), A_op->A, Vful, -st::one(), Vtmp, iflag), *iflag); \
  _MT_ normVtmp[nV]; \
  PHIST_CHK_IERR(SUBR( mvec_norm2 ) (Vtmp, normVtmp, iflag), *iflag); \
  _MT_ equalEps = normVtmp[0]; \
  for(int i = 0; i < nV; i++) \
    equalEps = mt::max(equalEps, normVtmp[i]); \
  PHIST_OUT(PHIST_INFO, "Line %d: A*V - W*H_A: %e\n", __LINE__, equalEps); \
  /* check H = W'*V */ \
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
  mvec_ptr newQ = NULL;
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

  mvec_ptr newBQ = NULL;
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

    if( symmetric)
    {
      // we have V'AV=V'WH_A=H'H_A symmetric, check this:
      PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Htmp_,&Htmp,0,    nV-1,      0,     nV-1,      iflag), *iflag); \
      PHIST_CHK_IERR(SUBR( sdMatT_times_sdMat) (st::one(), Hful, H_Aful, st::zero(), Htmp, iflag),*iflag);
      PHIST_CHK_NEG_IERR(SUBR(sdMat_check_symmetry)(Htmp, tol, iflag), *iflag);
    }

    // calculate sorted (generalized) Schur form of (H,H_A) in (S_L,S_R,T,T_A)
    // to get (H,H_A) = (S_L T S_R', S_L T_A S_R')
    //
    // We have  Vful    = [Q V]
    //          Wful    = [Q W]                  (because AQ=QR)
    //          H_Aful  = [R H_Aq; 0 H_A]        (where the coefficients H_Aq and H_A come from
    //                                            orthogonalizing AV against Q and itself, resp.)
    //          Hful    = Wful'Vful=[I Hq; 0 H]  (where Hq, H come from orthogonalizing V against
    //                                            Q and itself, respectively).
    //
    // So we get the complete decomposition
    //
    // |I Hq|   |I  0 ||I Tq||I 0   |
    // |0 H | = |0 S_L||0 T ||0 S_R'|
    //
    // |R H_Aq|   |I  0 ||R T_Aq||I 0   |
    // |0 H_A | = |0 S_L||0 T_A ||0 S_R'|
    //
    // and hence T_Aq H_Aq*S_R, Tq=Hq*S_R
    
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
    phist_lidx offT = ldaT*nConvEig+nConvEig;
    phist_lidx offT_A = ldaT_A*nConvEig+nConvEig;
    phist_lidx offS_L = ldaS_L*nConvEig+nConvEig;
    phist_lidx offS_R = ldaS_R*nConvEig+nConvEig;
    PHIST_CHK_IERR(SUBR( GenSchurDecomp ) (T_raw+offT, ldaT, T_A_raw+offT_A, ldaT_A, 
                                           S_L_raw+offS_L, ldaS_L, S_R_raw+offS_R, ldaS_R,
                                           nV-nConvEig, nSelect, nSort, which, tol, ev_H+nConvEig, iflag), *iflag);
    // compute Rayleigh quotients w.r.t. the harmonic Ritz values
    for (int i=nConvEig; i<nV; i++)
    {
      rq_H[i] = st::conj(T_raw[i*ldaT+i])*T_A_raw[i*ldaT_A+i];
    }

    // we still need to add the missing parts of T, T_A, S_L, S_R
    if( nConvEig > 0 )
    {
      // upper left part of S_L, S_R, T (identity matrices)
      for(int i = 0; i < nConvEig; i++)
      {
        S_L_raw[ldaS_L*i+i] = st::one();
        S_R_raw[ldaS_R*i+i] = st::one();
        T_raw[ldaT*i+i] = st::one();
      }

      // upper left part of T_A = R
      PHIST_CHK_IERR(SUBR( sdMat_view_block ) (T_A_, &T_A, 0, nConvEig-1, 0, nConvEig-1, iflag), *iflag);
      PHIST_CHK_IERR(SUBR( sdMat_get_block  ) (R,    T_A,  0, nConvEig-1, 0, nConvEig-1, iflag), *iflag);

      // upper right part of T_A and T
      PHIST_CHK_IERR(SUBR( sdMat_view_block  ) (S_R_, &S_R,  nConvEig, nV-1,       nConvEig, nEig_-1, iflag), *iflag);
      PHIST_CHK_IERR(SUBR( sdMat_view_block  ) (H_,   &Hq,   0,        nConvEig-1, nConvEig, nEig_-1, iflag), *iflag);
      PHIST_CHK_IERR(SUBR( sdMat_view_block  ) (T_,   &Tq,   0,        nConvEig-1, nConvEig, nEig_-1, iflag), *iflag);
      PHIST_CHK_IERR(SUBR( sdMat_view_block  ) (H_A_, &H_Aq, 0,        nConvEig-1, nConvEig, nEig_-1,    iflag), *iflag);
      PHIST_CHK_IERR(SUBR( sdMat_view_block  ) (T_A_, &T_Aq, 0,        nConvEig-1, nConvEig, nEig_-1,    iflag), *iflag);

      PHIST_CHK_IERR(SUBR( sdMat_times_sdMat ) (st::one(), Hq, S_R, st::zero(), Tq, iflag), *iflag);
      PHIST_CHK_IERR(SUBR( sdMat_times_sdMat ) (st::one(), H_Aq, S_R, st::zero(), T_Aq, iflag), *iflag);
    }
    
#ifdef TESTING
{
  // check that H S_R = S_L T
  PHIST_CHK_IERR(SUBR( sdMat_view_block ) (S_L_,&S_L, 0, nV-1,    0, nEig_-1, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_view_block ) (S_R_,&S_R, 0, nV-1,    0, nEig_-1, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_view_block ) (T_,&T, 0, nEig_-1, 0, nEig_-1, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(sdMat_view_block)(Htmp_, &Htmp, 0, nV-1,    0, nEig_-1, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_view_block ) (T_A_,&T_A, 0, nEig_-1, 0, nEig_-1, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(sdMat_view_block)(H_Atmp_, &H_Atmp, 0, nV-1,    0, nEig_-1, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(sdMat_times_sdMat)(st::one(), Hful, S_R, st::zero(), Htmp, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(sdMat_times_sdMat)(-st::one(), S_L, T, st::one(), Htmp, iflag), *iflag);

  // check that H S_R = S_L T
  PHIST_CHK_IERR(SUBR( sdMat_view_block ) (T_A_,&T_A, 0, nEig_-1, 0, nEig_-1, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(sdMat_view_block)(H_Atmp_, &H_Atmp, 0, nV-1,    0, nEig_-1, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(sdMat_times_sdMat)(st::one(), H_Aful, S_R, st::zero(), H_Atmp, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(sdMat_times_sdMat)(-st::one(), S_L, T_A, st::one(), H_Atmp, iflag), *iflag);

  // get max norms
  _MT_ absErrT = mt::zero(), absErrT_A=mt::zero();
  for(int i = 0; i < nEig_; i++)
    for(int j = 0; j < nV; j++)
    {
      absErrT   = mt::max(absErrT,   st::abs(Htmp_raw[i*ldaHtmp+j]));
      absErrT_A = mt::max(absErrT_A, st::abs(H_Atmp_raw[i*ldaH_Atmp+j]));
    }
  PHIST_SOUT(PHIST_INFO, "H*S_R   - S_L*T  : %8.4e\n"
                         "H_A*S_R - S_L*T_A: %8.4e\n", absErrT,absErrT_A);
  //PHIST_CHK_IERR(SUBR(sdMat_print)(H, iflag), *iflag);
  //PHIST_CHK_IERR(SUBR(sdMat_print)(Q_H, iflag), *iflag);
  //PHIST_CHK_IERR(SUBR(sdMat_print)(R_H, iflag), *iflag);
  //PHIST_CHK_IERR(SUBR(sdMat_print)(Htmp, iflag), *iflag);

}
#endif

    // update approximate Schur form of A (keeping the already computed part locked)
    // additional vectors q are present in Q so that Q has a total of nEig_ columns,
    // of which the first nConvEig are `hard locked', i.e. we do not rebuild them   
    // anymore.
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (Q_,   &Qq,                   nConvEig, nEig_-1,   iflag), *iflag);
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BQ_,  &BQq,                  nConvEig, nEig_-1,   iflag), *iflag);

    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (S_R_,&S_R,nConvEig,nV-1,      nConvEig, nEig_-1,   iflag), *iflag);
    PHIST_CHK_IERR(SUBR( mvec_times_sdMat ) (st::one(), V,    S_R, st::zero(), Qq,  iflag), *iflag);

    if( B_op != NULL )
    {
      PHIST_CHK_IERR(SUBR( mvec_times_sdMat ) (st::one(), BV,   S_R, st::zero(), BQq,  iflag), *iflag);
    }

    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (T_A_,&T_Aq, 0, nEig_-1, nConvEig, nEig_-1, iflag), *iflag);
    PHIST_CHK_IERR(SUBR( sdMat_set_block  ) (R,    T_Aq, 0, nEig_-1, nConvEig, nEig_-1, iflag), *iflag);

    ////////////////////////////////////////////////////////////////////////////////////
    // compute the residual for all the unconverged part q,r of the QR decomposition: //
    //, -res = -(Aq - B*q*r) = Bqr-Aq
    ////////////////////////////////////////////////////////////////////////////
    
    // first set t_res = Aq = AV S_R = W H_A S_R
    int ncolsR=nEig_-nConvEig;
    int nqv=nV-nConvEig;
    PHIST_CHK_IERR(SUBR( mvec_view_block  ) (res, &t_res, 0, ncolsR-1,   iflag), *iflag);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Htmp_,&Htmp, 0, nqv-1, 0, ncolsR-1, iflag), *iflag);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_A_,&H_A,   nConvEig, nV-1, nConvEig, nV-1, iflag), *iflag);
    PHIST_CHK_IERR(SUBR( sdMat_times_sdMat) (st::one(), H_A, S_R, st::zero(), Htmp, iflag), *iflag);
    PHIST_CHK_IERR(SUBR( mvec_times_sdMat ) (st::one(), W,   Htmp, st::zero(), t_res, iflag), *iflag);
    // overwrite t_res with the residual: -res = -(Aq - Bqr) = + BQq*r - res
    PHIST_CHK_IERR(SUBR( mvec_times_sdMat ) (st::one(), BQ,   T_Aq, -st::one(), t_res, iflag), *iflag);
    // calculate norm of the residual
    PHIST_CHK_IERR(SUBR( mvec_norm2 ) (t_res, resNorm+nConvEig, iflag), *iflag);

    // reset views
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (S_L_, &S_L, 0,     nV-1,      0,             nV-1,      iflag), *iflag);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (S_R_, &S_R, 0,     nV-1,      0,             nV-1,      iflag), *iflag);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (T_  , &T,   0,     nV-1,      0,             nV-1,      iflag), *iflag);
    PHIST_CHK_IERR(SUBR( sdMat_view_block ) (T_A_, &T_A, 0,     nV-1,      0,             nV-1,      iflag), *iflag);

#if PHIST_OUTLEV>=PHIST_DEBUG
  PHIST_SOUT(PHIST_DEBUG,"nV=%d; nConvEig=%d, nEig_=%d\n",nV,nConvEig,nEig_);
  PHIST_SOUT(PHIST_DEBUG,"H=\n");
  PHIST_CHK_IERR(SUBR(sdMat_print)(H,iflag),*iflag);
  PHIST_SOUT(PHIST_DEBUG,"H_A=\n");
  PHIST_CHK_IERR(SUBR(sdMat_print)(H_A,iflag),*iflag);
  PHIST_SOUT(PHIST_DEBUG,"S_L=\n");
  PHIST_CHK_IERR(SUBR(sdMat_print)(S_L,iflag),*iflag);
  PHIST_SOUT(PHIST_DEBUG,"S_R=\n");
  PHIST_CHK_IERR(SUBR(sdMat_print)(S_R,iflag),*iflag);
  PHIST_SOUT(PHIST_DEBUG,"T=\n");
  PHIST_CHK_IERR(SUBR(sdMat_print)(T,iflag),*iflag);
  PHIST_SOUT(PHIST_DEBUG,"T_A=\n");
  PHIST_CHK_IERR(SUBR(sdMat_print)(T_A,iflag),*iflag);
  PHIST_SOUT(PHIST_DEBUG,"Vful=\n");
  PHIST_CHK_IERR(SUBR(mvec_print)(Vful,iflag),*iflag);
  PHIST_SOUT(PHIST_DEBUG,"Wful=\n");
  PHIST_CHK_IERR(SUBR(mvec_print)(Wful,iflag),*iflag);
  PHIST_SOUT(PHIST_DEBUG,"Q_=\n");
  PHIST_CHK_IERR(SUBR(mvec_print)(Q_,iflag),*iflag);
#endif                                           

    
#ifdef TESTING
{
  // check that the residual is orthogonal to Q (should be by construction!)
  PHIST_CHK_IERR( SUBR( sdMat_view_block ) (Htmp_, &Htmp, 0, nEig_-1, nConvEig, nEig_-1, iflag), *iflag);
  PHIST_CHK_IERR( SUBR( mvecT_times_mvec ) (st::one(), Q, t_res, st::zero(), Htmp, iflag), *iflag);
  int nQ;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(Q,&nQ,iflag),*iflag);
  _MT_ equalEps = st::abs(Htmp_raw[0]);
  for(int i = 0; i < nEig_; i++)
    for(int j = nConvEig; j < nEig_; j++)
      equalEps = mt::max(equalEps, st::abs(Htmp_raw[i*ldaHtmp+j]));
  PHIST_OUT(PHIST_INFO, "Res orthogonality wrt. Q: %e (nQ=%d)\n", equalEps,nQ);
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


    // TODO - subspacejada has a (currently disabled) feature to
    // reorder multiple eigenvalues in the Schur form by residual norm,
    // this is not implemented here, yet, until we fix the implementation
    // there. The function ReorderPartialGenSchurDecomp is already available,
    // but not tested, in phist_schur_decomp.h
    std::vector<int> resPermutation(nEig_);
    for(int i = 0; i < nEig_; i++)
      resPermutation[i] = i;

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
      // determine which part of S_R we need to apply
      int imin, imax, jmin, jmax;
      if (symmetric)
      {
        imin=nConvEig;
        jmin=nConvEig;
      }
      else
      {
        imin=0;
        jmin=0;
      }

      // to avoid unnecessary subspace transformations, shrink the searchspace one iteration earlier...
      if( nV + 2*blockDim > maxBase )
      {
        PHIST_SOUT(PHIST_INFO,"Shrinking search space (one iteration earlier) from %d to %d\n", nV, minBase);
        imax=nV-1;
        jmax=minBase-1;
      }
      else
      {
          imax=nV-1;
        jmax=nV-1;
      }

      PHIST_CHK_IERR(SUBR( sdMat_view_block ) (S_L_,  &S_L,  imin, imax, jmin, jmax,    iflag), *iflag);
      PHIST_CHK_IERR(SUBR( sdMat_view_block ) (S_R_,  &S_R,  imin, imax, jmin, jmax,    iflag), *iflag);

      // copy Q into the first columns of W_ as well so that we orthogonalize W against Q when calling orthog below
      PHIST_CHK_IERR(SUBR(mvec_view_block) (Q_,   &Qq, nConvEig, nConvEig+nNewConvEig-1,   iflag), *iflag);
      PHIST_CHK_IERR(SUBR(mvec_set_block ) (W_,    Qq, nConvEig, nConvEig+nNewConvEig-1,iflag),*iflag);

      // reorder V and H TODO: can we exploit symmetry?
      if( symmetric && false)
      {
        PHIST_CHK_IERR(SUBR( transform_searchSpaceHarmonic ) (V, W, BV, H, H_A, S_L, S_R, B_op != NULL, iflag), *iflag);
      }
      else
      {
        PHIST_CHK_IERR(SUBR( transform_searchSpaceHarmonic ) (Vful, Wful, BVful, Hful, H_Aful, S_L, S_R, B_op != NULL, iflag), *iflag);
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
          sigma[k] = ct::real(rq_H[i]);
#else
          sigma[k] = rq_H[i];
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
        int imin=0, imax=nV-1, jmin=0, jmax=minBase-1;
        if( symmetric )
        {
          imin=nConvEig;
          jmin=nConvEig;
        }
        
        PHIST_CHK_IERR(SUBR( sdMat_view_block ) (S_L_,  &S_L,  imin, imax, jmin, jmax, iflag), *iflag);
        PHIST_CHK_IERR(SUBR( sdMat_view_block ) (S_R_,  &S_R,  imin, imax, jmin, jmax, iflag), *iflag);
        
        if (symmetric && false)
        {
          PHIST_CHK_IERR(SUBR( transform_searchSpaceHarmonic ) (V, W, BV, H, H_A, S_L, S_R, B_op != NULL, iflag), *iflag);
        }
        else
        {
          PHIST_CHK_IERR(SUBR( transform_searchSpaceHarmonic ) (Vful, Wful, BVful, Hful, H_Aful, S_L, S_R, B_op != NULL, iflag), *iflag);
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
      PHIST_CHK_IERR(SUBR( mvec_view_block  ) (At_,&Wv, 0, k-1, iflag), *iflag);

      // orthogonalize t against [Q V]
      PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Htmp_,  &HVv, 0,     nV-1,      nV,    nV+k-1,    iflag), *iflag);
      PHIST_CHK_IERR(SUBR( sdMat_view_block ) (Htmp_,  &Hvv, nV,    nV+k-1,    nV,    nV+k-1,    iflag), *iflag);
      int rankV;
      *iflag=PHIST_ROBUST_REDUCTIONS;
      PHIST_CHK_NEG_IERR(SUBR( orthog ) (Vful, Vv, B_op, Hvv, HVv, 5, &rankV, iflag), *iflag);

      // calculate AVv, BVv
      PHIST_CHK_IERR( A_op->apply(st::one(), A_op->A, Vv, st::zero(), Wv, iflag), *iflag);
      if( B_op != NULL )
      {
        //TODO: make orthog use and update BVv
        PHIST_CHK_IERR(SUBR( mvec_view_block  ) (BV_, &BVv,                   nV,    nV+k-1,    iflag), *iflag);
        PHIST_CHK_IERR( B_op->apply(st::one(), B_op->A, Vv, st::zero(), BVv, iflag), *iflag);
      }

      // enlarge search space
      // first update views
      PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_A_,  &HVv, 0,     nV-1,      nV,    nV+k-1,    iflag), *iflag);
      PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_A_,  &HvV, nV,    nV+k-1,    0,     nV-1,      iflag), *iflag);
      PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_A_,  &Hvv, nV,    nV+k-1,    nV,    nV+k-1,    iflag), *iflag);
      // orthogonalize At against [Q W] to get Wv and fill new column of H_A
      PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_A_,&H_Aq, 0,     nV-1,    nV,    nV+k-1,    iflag), *iflag);
      PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_A_,&H_A, nV,   nV+k-1,    nV,    nV+k-1,    iflag), *iflag);
      int rankW;
      *iflag=PHIST_ROBUST_REDUCTIONS;
      PHIST_CHK_NEG_IERR(SUBR( orthog ) (Wful, Wv, B_op, H_A, H_Aq, 5, &rankW, iflag), *iflag);
      // TODO: only take non-random vector if *iflag > 0

      // update H=W'V
      PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,  &HVv, 0,     nV-1,      nV,    nV+k-1,    iflag), *iflag);
      PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,  &HvV, nV,    nV+k-1,    0,     nV-1,      iflag), *iflag);
      PHIST_CHK_IERR(SUBR( sdMat_view_block ) (H_,  &Hvv, nV,    nV+k-1,    nV,    nV+k-1,    iflag), *iflag);

      // For A=A', H=W'V is not symmetric here, in contrast to subspacejada. H'H_A is, but I'm not sure
      // if we can exploit that here to save a reduction
      PHIST_CHK_IERR(SUBR( mvecT_times_mvec ) (st::one(), Wful,  Vv, st::zero(), HVv, iflag), *iflag);      
      PHIST_CHK_IERR(SUBR( mvecT_times_mvec ) (st::one(), Wv, Vful,  st::zero(), HvV, iflag), *iflag);
      PHIST_CHK_IERR(SUBR( mvecT_times_mvec ) (st::one(), Wv, Vv, st::zero(), Hvv, iflag), *iflag);
      
      // use set block to put Vv and Wv really into V and W
      PHIST_CHK_IERR(SUBR( mvec_set_block ) (V_, Vv,  nV, nV+k-1, iflag), *iflag);
      PHIST_CHK_IERR(SUBR( mvec_set_block ) (W_, Wv, nV, nV+k-1, iflag), *iflag);
      // increase nV
      nV = nV + k;
    } // k > 0

    UPDATE_SUBSPACE_VIEWS;
TESTING_CHECK_SUBSPACE_INVARIANTS;
  }

  // copy result to Q_
  PHIST_CHK_IERR(SUBR(mvec_set_block)(Q__, Q, 0, nEig_-1, iflag), *iflag);
  // copy resulting eigenvalues to ev
  for(int i = 0; i < nEig_; i++) ev[i] = ev_H[i];


  //------------------------------- delete vectors and matrices --------------------
  PHIST_CHK_IERR(SUBR(jadaCorrectionSolver_delete)(innerSolv, iflag), *iflag);

  // delete views
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (H, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (H_A, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Hq, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (H_Aq, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (T, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (T_A, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Tq, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (T_Aq, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (S_L,iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (S_R, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Hvv, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (HvV, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (HVv, iflag), *iflag);

  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Hful,iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (H_Aful,iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Htmp,iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (R,   iflag), *iflag);

  PHIST_CHK_IERR(SUBR( mvec_delete  ) (Qq,  iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (BQq, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (t,   iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (t_res,iflag),*iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (BVv, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (BV,  iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (Wv, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (W,  iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (Vv,  iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (Vtmp,iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (V,   iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (Vful,iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (Wful,iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (BVful,iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (Qtil,iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (BQtil,iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (Q,   iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (BQ,  iflag), *iflag);

  delete[] ev_H;
  delete[] rq_H;

  // delete mvecs and sdMats
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (S_L_,iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (S_R_,iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (T_,iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (T_A_,iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (H_,  iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (H_A_,  iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Htmp_,iflag), *iflag);
  if( B_op != NULL )
  {
    PHIST_CHK_IERR(SUBR( mvec_delete )(BV_, iflag), *iflag);
    PHIST_CHK_IERR(SUBR( mvec_delete )(BQ_, iflag), *iflag);
  }
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (Q_,  iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (t_,  iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (At_, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (res, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (W_, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (Vtmp_,iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_delete  ) (V_,  iflag), *iflag);
}

