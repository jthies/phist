#include "phist_blockedgmres_helper_def.hpp"

// small C++ wrapper that deletes temporary vectors at the end of a scope
class TYPE(mvec_GC)
{
  public:

  TYPE(mvec_GC)()
  {
    x_=NULL;
    y_=NULL;
    z_=NULL;
  }
  
  ~TYPE(mvec_GC)()
  {
    int iflag;
    if (z_==x_ || z_==y_) z_=NULL;
    if (z_) SUBR(mvec_delete)(z_,&iflag);
    if (x_) SUBR(mvec_delete)(x_,&iflag);
    if (y_) SUBR(mvec_delete)(y_,&iflag);
    x_=NULL;
    y_=NULL;
    z_=NULL;
  }
  
  //! the actual objects
  TYPE(mvec_ptr) x_,y_,z_;

};

// create new state objects. We just get an array of (NULL-)pointers
void SUBR(blockedGMRESstates_create)(TYPE(blockedGMRESstate_ptr) state[], int numSys, phist_const_map_ptr map, int maxBas,int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;

  if (numSys <= 0)
    return;

  phist_const_comm_ptr comm;
  PHIST_CHK_IERR(phist_map_get_comm(map,&comm,iflag),*iflag);

  // setup buffer of mvecs to be used later
#ifdef PHIST_HAVE_TEUCHOS
  Teuchos::RCP<TYPE(MvecRingBuffer)> mvecBuff(new TYPE(MvecRingBuffer)(maxBas+1));
#else
  TYPE(MvecRingBuffer)* mvecBuff = new TYPE(MvecRingBuffer)(maxBas+1);
#endif
  PHIST_CHK_IERR( mvecBuff->create_mvecs(map, numSys, iflag), *iflag);

  // set up individual states
  for(int i = 0; i < numSys; i++)
  {
    // initialization data for the next state
    TYPE(blockedGMRESstate) tmp = {i,(_MT_)0.5,-2,0,-1,0,NULL,NULL,NULL,NULL,NULL,-mt::one(),-mt::one(),st::zero(),NULL};

    // create state
    state[i] = new TYPE(blockedGMRESstate)(tmp);

    // allocate members
    PHIST_CHK_IERR(SUBR( sdMat_create )(&state[i]->H_, maxBas+1, maxBas, comm, iflag), *iflag);
    PHIST_CHK_IERR(SUBR( mvec_create  )(&state[i]->b_, map,      1,            iflag), *iflag);
    state[i]->cs_ = new ST[maxBas];
    state[i]->sn_ = new ST[maxBas];
    state[i]->rs_ = new ST[maxBas+1];

#ifdef PHIST_HAVE_TEUCHOS
    // assign MvecRingBuffer (with reference counting)
    state[i]->Vbuff = (void*) new Teuchos::RCP<TYPE(MvecRingBuffer)>(mvecBuff);
#else
    // TODO - check memory management, this used to be an RCP
    state[i]->Vbuff = (void*)mvecBuff;
#endif
  }
}


// delete blockedGMRESstate object
void SUBR(blockedGMRESstates_delete)(TYPE(blockedGMRESstate_ptr) state[], int numSys, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
#ifndef PHIST_HAVE_TEUCHOS
  if (numSys==0) return;
  
    PHIST_CAST_PTR_FROM_VOID(TYPE(MvecRingBuffer), mvecBuff, state[0]->Vbuff, *iflag);    
    PHIST_CHK_IERR(mvecBuff->delete_mvecs(iflag), *iflag);
    
    delete mvecBuff;
#endif

  for(int i = 0; i < numSys; i++)
  {
    PHIST_CHK_IERR(SUBR( sdMat_delete ) (state[i]->H_, iflag), *iflag);
    PHIST_CHK_IERR(SUBR( mvec_delete  ) (state[i]->b_, iflag), *iflag);
    delete [] state[i]->cs_;
    delete [] state[i]->sn_;
    delete [] state[i]->rs_;
#ifdef PHIST_HAVE_TEUCHOS
    PHIST_CAST_PTR_FROM_VOID(Teuchos::RCP<TYPE(MvecRingBuffer)>, mvecBuff, state[i]->Vbuff, *iflag);
    PHIST_CHK_IERR((*mvecBuff)->delete_mvecs(iflag), *iflag);
    delete mvecBuff;
#endif
    delete state[i];
  }
}


// reset blockedGMRES state.
void SUBR(blockedGMRESstate_reset)(TYPE(blockedGMRESstate_ptr) S, TYPE(const_mvec_ptr) b, TYPE(const_mvec_ptr) x0, int *iflag)
{
#include "phist_std_typedefs.hpp"  
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
  
  int previous_status = S->status;
  S->status=-2; // not initialized, if this function fails in some way to set status the next _iterate call will complain

  // get mvecBuff
#ifdef PHIST_HAVE_TEUCHOS
  PHIST_CAST_PTR_FROM_VOID(Teuchos::RCP<TYPE(MvecRingBuffer)>, mvecBuffPtr, S->Vbuff, *iflag);
  Teuchos::RCP<TYPE(MvecRingBuffer)> mvecBuff = *mvecBuffPtr;
#else
  PHIST_CAST_PTR_FROM_VOID(TYPE(MvecRingBuffer), mvecBuff, S->Vbuff, *iflag);
#endif


  // release mvecs currently marked used by this state
  for(int j = 0; j < S->curDimV_; j++)
  {
    int Vind = mvecBuff->prevIndex(S->lastVind_,j);
    mvecBuff->decRef(Vind);
  }

  S->curDimV_ = 0;

  // set H to zero
  PHIST_CHK_IERR(SUBR(sdMat_put_value)(S->H_, st::zero(), iflag), *iflag);


  // only freed resources, object still not initialized (status -2)
  if( b == NULL && x0 == NULL )
  {
    S->status=-2;
    return;
  }

  if( b == NULL && (S->normR0_ == -mt::one() || previous_status == -2) )
  {
    PHIST_OUT(PHIST_ERROR,"on the first call to blockedGMRESstate_reset you *must* provide the RHS vector");
    *iflag=-1;
    return;
  }

  if( b != NULL )
  {
    // new rhs -> need to recompute ||b-A*x0||
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(), b, st::zero(), S->b_, iflag), *iflag);
    S->status = -1; // indicates that the residual has to be calculated
    S->totalIter = 0;
    S->normR0_ = -mt::one(); // indicates the initial residual has to be set (first start)
  }

  if( x0 == NULL )
  {
    // great, we can directly apply a first gmres step as we don't need to compute A*x0

    PHIST_CHK_IERR(SUBR( mvec_norm2 ) (S->b_, &S->normR_, iflag), *iflag);
    if( S->normR0_ < mt::zero() )
      S->normR0_ = S->normR_;
    S->rs_[0] = S->normR_;
    S->prevBeta_ = S->normR_;

    S->lastVind_ = mvecBuff->lastIndex();
    S->curDimV_ = 1;
    TYPE(mvec_ptr) r0 = NULL;
    PHIST_CHK_IERR(SUBR( mvec_view_block )(mvecBuff->at(S->lastVind_), &r0, S->id, S->id, iflag), *iflag);
    mvecBuff->incRef(S->lastVind_);
    if( S->normR_ != mt::zero() )
    {
      _ST_ scale = st::one() / S->normR_;
      PHIST_CHK_IERR(SUBR( mvec_add_mvec ) (scale, S->b_, st::zero(), r0, iflag), *iflag);
    }
    PHIST_CHK_IERR(SUBR( mvec_delete ) (r0, iflag), *iflag);
    S->status=1; // 1: iterating/unconverged
  }
  else // x != NULL
  {
    // initialize everything to calculate b-A*x0 in the next call to iterate
    S->lastVind_ = mvecBuff->lastIndex();
    PHIST_CHK_IERR(SUBR( mvec_set_block ) (mvecBuff->at(S->lastVind_), x0, S->id, S->id, iflag), *iflag);
    mvecBuff->incRef(S->lastVind_);
    S->prevBeta_ = st::zero();
    S->status=-1; // restart, need to compute residual
  }

}


// calculate approximate solution.
// first solves triangular system s=R\y for all states and then updates the solution of all the given
// states in a vectorized way, x+=Vs
void SUBR(blockedGMRESstates_updateSol)(TYPE(blockedGMRESstate_ptr) S[], int numSys, 
        TYPE(const_linearOp_ptr) rightPrecon,
        TYPE(mvec_ptr) x, _MT_* resNorm, bool scaleSolutionToOne, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag = 0;

  if( numSys <= 0 )
    return;

  // if there are multiple systems, get the maximal id
  int maxId = 0;
  for(int i = 0; i < numSys; i++)
    maxId = std::max(maxId,S[i]->id);
    
  // if there is no (right) preconditioner, update x+=Vs directly, otherwise
  // compute z=Vs first and then update x+=P\z
  TYPE(mvec_ptr) z=x;
  if (rightPrecon!=NULL)
  {
    z=NULL;
    PHIST_CHK_IERR(SUBR(mvec_create)(&z,rightPrecon->range_map,numSys,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_put_value)(z,st::zero(),iflag),*iflag);
  }

  bool ordered = true;
  for(int i = 0; i < numSys; i++)
    ordered = ordered && (S[i]->id == i);

  // make sure all lastVind_ are the same
  int lastVind = S[0]->lastVind_;
  for(int i = 0; i < numSys; i++)
    PHIST_CHK_IERR(*iflag = (S[i]->lastVind_ != lastVind) ? -1 : 0, *iflag);

#ifdef PHIST_HAVE_TEUCHOS
  PHIST_CAST_PTR_FROM_VOID(Teuchos::RCP<TYPE(MvecRingBuffer)>, mvecBuffPtr, S[0]->Vbuff, *iflag);
  Teuchos::RCP<TYPE(MvecRingBuffer)> mvecBuff = *mvecBuffPtr;
#else
  PHIST_CAST_PTR_FROM_VOID(TYPE(MvecRingBuffer), mvecBuff, S[0]->Vbuff, *iflag);
#endif

  // make sure all systems use the same mvecBuff
  for(int i = 0; i < numSys; i++)
  {
#ifdef PHIST_HAVE_TEUCHOS
    PHIST_CAST_PTR_FROM_VOID(Teuchos::RCP<TYPE(MvecRingBuffer)>, mvecBuffPtr_i, S[i]->Vbuff, *iflag);
    PHIST_CHK_IERR(*iflag = (*mvecBuffPtr_i != *mvecBuffPtr) ? -1 : 0, *iflag);
#else
    PHIST_CAST_PTR_FROM_VOID(TYPE(MvecRingBuffer), mvecBuffPtr_i, S[i]->Vbuff, *iflag);
    PHIST_CHK_IERR(*iflag = (mvecBuffPtr_i != mvecBuff) ? -1 : 0, *iflag);
#endif
  }

  // determine maximal and shared (minimal) dimensions of subspaces
  int maxCurDimV = 0;
  int sharedCurDimV = mvecBuff->size();
  for(int i = 0; i < numSys; i++)
  {
    maxCurDimV = std::max(maxCurDimV, S[i]->curDimV_);
    sharedCurDimV = std::min(sharedCurDimV, S[i]->curDimV_);
  }

  // calculate resNorm
  for(int i = 0; i < numSys; i++)
  {
    if( S[i]->curDimV_ <= 0 )
      resNorm[i] = -mt::one();
    else if( S[i]->curDimV_ == 1 )
      resNorm[i] = mt::one();
    else
      resNorm[i] = S[i]->normR_/S[i]->normR0_;
  }

  // no iteration done yet?
  if( maxCurDimV <= 1 ) return;


  // allocate space for y
  _ST_ *yglob = new _ST_[(maxId+1)*(maxCurDimV-1)];
  for(int i = 0; i < (maxId+1)*(maxCurDimV-1); i++)
    yglob[i] = st::zero();
  phist_lidx ldy = (maxId+1);

  // calculate y by solving the triangular systems
  for(int i = 0; i < numSys; i++)
  {
    // nothing to do here?
    if( S[i]->curDimV_ <= 1 )
      continue;

    // helpful variables
    int m = S[i]->curDimV_-1;
    ST *H_raw=NULL;
    phist_lidx ldH;
    PHIST_CHK_IERR(SUBR(sdMat_extract_view)(S[i]->H_,&H_raw,&ldH,iflag),*iflag);
    _ST_ *y = &yglob[S[i]->id+ldy*(maxCurDimV-S[i]->curDimV_)];


#if 0
//PHIST_OUTLEV>=PHIST_DEBUG
    PHIST_SOUT(PHIST_DEBUG,"blockedGMRES_updateSol[%d], curDimV=%d, H=\n",i,S[i]->curDimV_);
    {
      TYPE(sdMat_ptr) H = NULL;
      PHIST_CHK_IERR(SUBR(sdMat_view_block)(S[i]->H_, &H, 0, m-1, 0, m-1, iflag), *iflag);
      PHIST_CHK_IERR(SUBR(sdMat_print)(H,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_delete)(H,iflag),*iflag);
    }
    PHIST_SOUT(PHIST_DEBUG,"rs=\n");
    for (int k=0;k<m;k++)
    {
      PHIST_SOUT(PHIST_DEBUG,"%16.8f+%16.8fi\n",st::real(S[i]->rs_[k]),st::imag(S[i]->rs_[k]));      
    }
#endif

    // set y to rs
    for(int j = 0; j < m; j++)
      y[ldy*j] = S[i]->rs_[j];


#if defined(TESTING) && (PHIST_OUTLEV>=PHIST_DEBUG)
{
  // setup e-vector
  _ST_ e[m+1];
  e[0] = st::one();
  for(int j = 1; j < m+1; j++)
    e[j] = st::zero();
  // apply givens rotations
  for(int j = 0; j < m; j++)
  {
    _ST_ tmp = st::conj(S[i]->cs_[j]) * e[j]  +  st::conj(S[i]->sn_[j]) * e[j+1];
    e[j+1] = -S[i]->sn_[j] * e[j]  +  S[i]->cs_[j] * e[j+1];
    e[j] = tmp;
  }
  // check that y is givens_rotations applied to e
  PHIST_SOUT(PHIST_INFO, "rs/norm0:");
  for(int j = 0; j < m; j++)
  {
    PHIST_SOUT(PHIST_INFO, "\t%8.4e + i%8.4e", st::real(y[ldy*j])/S[i]->normR0_, st::imag(y[ldy*j])/S[i]->normR0_);
  }
  PHIST_SOUT(PHIST_INFO, "\nabs(rs(j)/norm0)):%8.4e", st::abs(y[ldy*(m-1)])/S[i]->normR0_);
  PHIST_SOUT(PHIST_INFO, "\nrot(e_1):");
  for(int j = 0; j <= m; j++)
  {
    PHIST_SOUT(PHIST_INFO, "\t%8.4e + i%8.4e", st::real(e[j]), st::imag(e[j]));
  }
  PHIST_SOUT(PHIST_INFO, "\nabs(rot(e_1)):%8.4e\n", st::abs(e[m]));
}
#endif


    // solve triangular system
    phist_blas_idx ildH=static_cast<phist_blas_idx>(ldH);
    phist_blas_idx ildy=static_cast<phist_blas_idx>(ldy);
    if (ildH<0 || ildy<0)
    {
      *iflag=PHIST_INTEGER_OVERFLOW;
      return;
    }
    PHIST_CHK_IERR(SUBR(sdMat_from_device)(S[i]->H_,iflag),*iflag);
    PHIST_CHK_IERR(PHIST_TG_PREFIX(TRSV)("U","N","N",&m,(st::blas_scalar_t*)H_raw,&ildH,(st::blas_scalar_t*)y, &ildy),*iflag);


    // if we are only interested in the directions Vi*yi and appropriate AVi*yi,
    // then this scaling may help to improve the conditioning of a following orthogonalization step!
    if( scaleSolutionToOne )
    {
      // scale y to one
      _MT_ scale = mt::zero();
      for(int j = 0; j < m; j++)
        scale += st::real(st::conj(y[ldy*j])*y[ldy*j]);
      scale = mt::one()/sqrt(scale);
      PHIST_DEB("scaling solution with: %8.4e\n", scale);
      for(int j = 0; j < m; j++)
        y[ldy*j] *= scale;
    }
  }


// put all iterations in one big compute task; this speeds up the tests with ghost (significantly)
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
  // add up solution
  TYPE(mvec_ptr) Vj = NULL, x_i = NULL;
  PHIST_DEB("blockedGMRESstates_updateSol maxCurDimV = %d, sharedCurDimV = %d, ordered = %d\n", maxCurDimV, sharedCurDimV, (int)ordered);
  for(int j = 0; j < maxCurDimV-1; j++)
  {
    int Vind = mvecBuff->prevIndex(lastVind,maxCurDimV-1-j);
    _ST_ *yj = yglob + ldy*j;
    PHIST_DEB("blockedGMRESstates_updateSol j = %d, Vind = %d, yj = ", j, Vind);
    for(int i = 0; i < numSys; i++)
    {
#ifdef IS_COMPLEX
      PHIST_DEB("\t%e%+ei (%d)", st::real(yj[S[i]->id]),
                                 st::imag(yj[S[i]->id]), S[i]->id);
#else
      PHIST_DEB("\t%e (%d)", yj[S[i]->id], S[i]->id);
#endif
    }
    PHIST_DEB("\n");

    // compute z=V*y (if right precond is present) or
    // x=x+V*y if there is no right preconditioning (z points to x in that case)
    if( j >= maxCurDimV-sharedCurDimV && ordered )
    {
      // update solution of all systems at once
      PHIST_CHK_IERR(SUBR(mvec_view_block)(mvecBuff->at(Vind), &Vj, 0, maxId, iflag), *iflag);
      PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(yj, Vj, st::one(), z, iflag), *iflag);
    }
    else
    {
      // update solution of single systems
      for(int i = 0; i < numSys; i++)
      {
        if( j >= maxCurDimV-S[i]->curDimV_ )
        {
          PHIST_CHK_IERR(SUBR(mvec_view_block)(mvecBuff->at(Vind), &Vj, S[i]->id, S[i]->id, iflag), *iflag);
          PHIST_CHK_IERR(SUBR(mvec_view_block)(x, &x_i, i, i, iflag), *iflag);
          PHIST_CHK_IERR(SUBR(mvec_add_mvec)(yj[S[i]->id], Vj, st::one(), z, iflag), *iflag);
        }
      }
    }
  }
 
  if (rightPrecon!=NULL)
  {
    PHIST_CHK_IERR(rightPrecon->apply(st::one(),rightPrecon->A,z,st::one(),x,iflag),*iflag);
  }

  PHIST_CHK_IERR(SUBR(mvec_delete)(x_i, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_delete)(Vj, iflag), *iflag);
  if (z!=x) PHIST_CHK_IERR(SUBR(mvec_delete)(z,iflag),*iflag);
PHIST_TASK_END(iflag)
  delete[] yglob;
}


// implementation of gmres on several systems simultaneously
void SUBR(blockedGMRESstates_iterate)(TYPE(const_linearOp_ptr) Aop, TYPE(const_linearOp_ptr) Pop,
        TYPE(blockedGMRESstate_ptr) S[], int numSys, int* nIter, bool useIMGS, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag = 0;

  int maxIter = (*nIter)>0 ? *nIter: 9999999;
  *nIter=0;

#if PHIST_OUTLEV>=PHIST_DEBUG
  PHIST_SOUT(PHIST_DEBUG,"starting function iterate() with %d systems\n curDimVs: ",numSys);
  for (int i=0;i<numSys;i++)
  {
    PHIST_SOUT(PHIST_DEBUG,"%d ",S[i]->curDimV_);
  }
  PHIST_SOUT(PHIST_DEBUG,"\n");
#endif

  if( numSys <= 0 )
    return;

  // if there are multiple systems, get the maximal id
  int maxId = 0;
  for(int i = 0; i < numSys; i++) maxId = std::max(maxId,S[i]->id);
  int minId = maxId;
  for(int i = 0; i < numSys; i++) minId = std::min(minId,S[i]->id);

#ifdef PHIST_HAVE_TEUCHOS
  PHIST_CAST_PTR_FROM_VOID(Teuchos::RCP<TYPE(MvecRingBuffer)>, mvecBuffPtr, S[0]->Vbuff, *iflag);
  Teuchos::RCP<TYPE(MvecRingBuffer)> mvecBuff = *mvecBuffPtr;
#else
  PHIST_CAST_PTR_FROM_VOID(TYPE(MvecRingBuffer), mvecBuff, S[0]->Vbuff, *iflag);
#endif

  // make sure all systems use the same mvecBuff
  for(int i = 0; i < numSys; i++)
  {
#ifdef PHIST_HAVE_TEUCHOS
    PHIST_CAST_PTR_FROM_VOID(Teuchos::RCP<TYPE(MvecRingBuffer)>, mvecBuffPtr_i, S[i]->Vbuff, *iflag);
    PHIST_CHK_IERR(*iflag = (*mvecBuffPtr_i != *mvecBuffPtr) ? -1 : 0, *iflag);
#else
    PHIST_CAST_PTR_FROM_VOID(TYPE(MvecRingBuffer), mvecBuffPtr_i, S[i]->Vbuff, *iflag);
    PHIST_CHK_IERR(*iflag = (mvecBuffPtr_i != mvecBuff) ? -1 : 0, *iflag);
#endif
  }

  // determine maximal and shared (minimal) dimensions of subspaces
  int maxCurDimV = 0;
  int sharedCurDimV = mvecBuff->size();
  for(int i = 0; i < numSys; i++)
  {
    maxCurDimV = std::max(maxCurDimV, S[i]->curDimV_);
    sharedCurDimV = std::min(sharedCurDimV, S[i]->curDimV_);
  }

  // work vectors for x, y = Aop(x) and z=P\r
  TYPE(mvec_GC) work;
  
  // z is Precond\x, create it as a temporary vector
  if (Pop!=NULL)
  {
    PHIST_CHK_IERR(SUBR(mvec_create)(&work.z_,Pop->range_map,numSys,iflag),*iflag);
  }

  {
    // make sure all lastVind_ are the same
    int lastVind = S[0]->lastVind_;
    for(int i = 0; i < numSys; i++)
      PHIST_CHK_IERR(*iflag = (S[i]->lastVind_ != lastVind) ? -1 : 0, *iflag);

    // x0 / last element of krylov subspace
    PHIST_CHK_IERR(SUBR( mvec_view_block )( mvecBuff->at(lastVind), &work.x_, minId, maxId, iflag), *iflag);

#ifdef TESTING
// print a visualization of the current state
std::vector< std::vector<int> > mvecUsedBy(maxId+1);
for(int i = 0; i < maxId+1; i++)
  mvecUsedBy[i].resize(mvecBuff->size(),false);
std::vector<int> idUsed(maxId+1,false);
for(int i = 0; i < numSys; i++)
{
  idUsed[S[i]->id] = true;
  for(int j = 0; j < S[i]->curDimV_; j++)
  {
    int Vind = mvecBuff->prevIndex(lastVind,j);
    mvecUsedBy[S[i]->id][Vind] = true;
  }
}
PHIST_SOUT(PHIST_INFO,"Blocked GMRES status:\n");
for(int j = 0; j < mvecBuff->size(); j++)
{
  PHIST_SOUT(PHIST_INFO,"--");
}
PHIST_SOUT(PHIST_INFO,"\n");
for(int i = 0; i < maxId+1; i++)
{
  for(int j = 0; j < mvecBuff->size(); j++)
  {
    if( j == lastVind && idUsed[i] )
    {
      PHIST_SOUT(PHIST_INFO," |");
    }
    else if( mvecUsedBy[i][j] )
    {
      PHIST_SOUT(PHIST_INFO," +");
    }
    else
    {
      PHIST_SOUT(PHIST_INFO,"  ");
    }
  }
  PHIST_SOUT(PHIST_INFO,"\n");
}
for(int j = 0; j < mvecBuff->size(); j++)
{
  PHIST_SOUT(PHIST_INFO,"--");
}
PHIST_SOUT(PHIST_INFO,"\n");
for(int j = 0; j < mvecBuff->size(); j++)
{
  PHIST_SOUT(PHIST_INFO," %1d",mvecBuff->refCount(j));
}
PHIST_SOUT(PHIST_INFO,"\n");
#endif
  }


  // views into V_
  TYPE(mvec_ptr) Vj = NULL, Vk = NULL;
  // views into H_
  TYPE(sdMat_ptr) R1 = NULL, R2 = NULL;


  // we return as soon as one system converges or reaches its
  // maximum permitted number of iterations. The decision about what to do
  // next is then left to the caller.
  int anyConverged = 0;
  int anyFull      = 0;
  int anyFailed    = 0;



  PHIST_SOUT(PHIST_VERBOSE,"GMRES iteration started\n");
  PHIST_SOUT(PHIST_VERBOSE,"=======================\n");

  // check if there are already converged/failed systems
  for(int i = 0; i < numSys; i++)
  {
    // check convergence
    MT relres = S[i]->normR_ / S[i]->normR0_;
    MT absres = S[i]->normR_;
    if( S[i]->normR0_ != -mt::one() && ( absres < 1000*st::eps() || relres < S[i]->tol ) )
    {
      S[i]->status = 0; // mark as converged
      anyConverged++;
    }
    else if( S[i]->totalIter>=maxIter )
    {
      S[i]->status = 3;
      anyFailed++;
    }
    else if( S[i]->curDimV_ >= mvecBuff->size() )
    {
      S[i]->status = 2; // mark as failed/restart needed
      anyFull++;
    }
  }
  for (int i=0;i<numSys;i++)
  {
    if (S[i]->curDimV_>0)
    {
      PHIST_SOUT(PHIST_VERBOSE,"[%d]: %d\t%8.4e\t(%8.4e)\n", i, 
        S[i]->curDimV_-1,S[i]->normR_/S[i]->normR0_,S[i]->normR_);
    }
    else
    {
      PHIST_SOUT(PHIST_VERBOSE,"[%d]: restarted\n",i);
    }
  }

// put all iterations in one big compute task; this speeds up the tests with ghost (significantly)
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
  while( anyConverged==0 && anyFailed==0 && anyFull==0 )
  {
    //    % get new vector for y
    int nextIndex;
    PHIST_CHK_IERR( mvecBuff->getNextUnused(nextIndex,iflag), *iflag);
    PHIST_CHK_IERR(SUBR( mvec_view_block ) (mvecBuff->at(nextIndex), &work.y_, minId, maxId, iflag), *iflag);
    
    // apply (right) preconditioning *if this is not a restart*
    if (Pop!=NULL)
    {
      int numRestarted=0;
      for (int i=0; i<numSys; i++) 
      {
        numRestarted+=(S[i]->curDimV_==0)?1:0;
      }
      // if all systems were restarted, we don't need to apply the preconditioner in the first iteration at all
      // (in this 'iteration' only the residual is computed)
      if (numRestarted==numSys)
      {
        PHIST_CHK_IERR( SUBR(mvec_add_mvec)(st::one(),work.x_,st::zero(),work.z_,iflag),*iflag);
      }
      else
      {
        // apply preconditioner to all vectors for simplicity and overwrite the ones that should not be preconditioned
        // because the corresponding system was restarted
        PHIST_CHK_IERR( Pop->apply (st::one(), Pop->A, work.x_, st::zero(), work.z_, iflag), *iflag);
        if (numRestarted>0)
        {
          for (int i=0; i<numSys; i++)
          {
            int j=S[i]->curDimV_;
            if( j == 0 )
            {
              // (re-)start: r_0 = b - A*x_0, so throw away the P\x and 
              // replace by x for this column
              TYPE(mvec_ptr) tmpX=NULL, tmpZ=NULL;
              PHIST_CHK_IERR(SUBR(mvec_view_block)(work.x_,&tmpX,S[i]->id-minId,S[i]->id-minId,iflag),*iflag);
              PHIST_CHK_IERR(SUBR(mvec_view_block)(work.z_,&tmpZ,S[i]->id-minId,S[i]->id-minId,iflag),*iflag);
              PHIST_CHK_IERR( SUBR(mvec_add_mvec)(st::one(),tmpX,st::zero(),tmpZ,iflag),*iflag);
            }
          }
        }//numRestarted>0
      }//numRestarted<numSys
    }//Pop!=NULL
    else
    {
      work.z_=work.x_;
    }

    //    % apply operator A
    PHIST_CHK_IERR( Aop->apply (st::one(), Aop->A, work.z_, st::zero(), work.y_, iflag), *iflag);


    //    % initialize GMRES for (re-)started systems
    for(int i = 0; i < numSys; i++)
    {
      int j = S[i]->curDimV_;
      if( j == 0 )
      {
        // (re-)start: r_0 = b - A*x_0
        PHIST_CHK_IERR( SUBR(mvec_view_block) (work.y_, &Vj, S[i]->id-minId, S[i]->id-minId, iflag), *iflag);
        PHIST_CHK_IERR( SUBR(mvec_add_mvec) (st::one(), S[i]->b_, -st::one(), Vj, iflag), *iflag);
      }
    }
    // increment ref counters in mvecBuff and set lastVind_
    for(int i = 0; i < numSys; i++)
    {
      if( S[i]->curDimV_ == 0 )
      {
        // we have used one index for x0, release it
        mvecBuff->decRef(S[i]->lastVind_);
      }
      S[i]->lastVind_ = nextIndex;
      mvecBuff->incRef(nextIndex);
    }

    for(int i = 0; i < numSys; i++)
    {
      // get H daat from GPU (if it resides there)
      PHIST_CHK_IERR(SUBR(sdMat_from_device)(S[i]->H_,iflag),*iflag);
    }

    //    % Arnoldi update with iterated modified Gram Schmidt
    {
      std::vector<_MT_> ynorm(maxId+1-minId,-mt::one());
      std::vector<_MT_> prev_ynorm(maxId+1-minId,-mt::one());
      for(int mgsIter = 0; mgsIter < 3; mgsIter++)
      {
        // calculate norm
        prev_ynorm = ynorm;
        if( useIMGS || mgsIter == 1 )
        {
          PHIST_CHK_IERR(SUBR(mvec_norm2)(work.y_, &ynorm[0], iflag), *iflag);
        }

        bool needAnotherIteration = (mgsIter == 0);
        PHIST_SOUT(PHIST_DEBUG,"reduction in norm in IMGS:");
        for(int i = 0; i < numSys; i++)
        {
          PHIST_SOUT(PHIST_DEBUG,"\t%8.4e", ynorm[S[i]->id-minId]/prev_ynorm[S[i]->id-minId]);
          if( ynorm[S[i]->id-minId] < 0.75 * prev_ynorm[S[i]->id-minId] )
            needAnotherIteration = true;
        }
        PHIST_SOUT(PHIST_DEBUG,"\n");
        if( mgsIter > 0)
        {
          if( !(needAnotherIteration && useIMGS) )
            break;
          PHIST_SOUT(PHIST_DEBUG, "Additional MGS iteration in blockedGMRES!\n");
        }


        for(int j = 0; j < maxCurDimV; j++)
        {
          int Vind = mvecBuff->prevIndex(S[0]->lastVind_,maxCurDimV-j);
          _ST_ tmp[maxId+1-minId];

          bool calculatedDot = false;
          if( j >= maxCurDimV-sharedCurDimV && maxId+1-minId == numSys )
          {
            // MGS step for all systems at once
            PHIST_CHK_IERR(SUBR( mvec_view_block ) (mvecBuff->at(Vind), &Vk, minId, maxId, iflag), *iflag);
            PHIST_CHK_IERR(SUBR( mvec_dot_mvec   ) (Vk, work.y_, tmp, iflag), *iflag);
            for(int i = 0; i < maxId+1-minId; i++)
              tmp[i] = -tmp[i];
            PHIST_CHK_IERR(SUBR( mvec_vadd_mvec  ) (tmp, Vk, st::one(), work.y_, iflag), *iflag);
            calculatedDot = true;
          }

          // store in H (and do MGS steps for single systems)
          for(int i = 0; i < numSys; i++)
          {
            int j_ = j - (maxCurDimV-S[i]->curDimV_);
            if( j_ >= 0 )
            {
    //std::cout << "In blockedGMRES arnoldi: j " << j << " Vind " << Vind << " i " << i << " j_ " << j_ << " curDimV " << S[i]->curDimV_ << std::endl;
              if( !calculatedDot )
              {
                // MGS step for single system
                PHIST_CHK_IERR(SUBR( mvec_view_block ) (mvecBuff->at(Vind), &Vk, S[i]->id,       S[i]->id,       iflag), *iflag);
                PHIST_CHK_IERR(SUBR( mvec_view_block ) (work.y_,             &Vj, S[i]->id-minId, S[i]->id-minId, iflag), *iflag);
                PHIST_CHK_IERR(SUBR( mvec_dot_mvec   ) (Vk, Vj, &tmp[S[i]->id-minId], iflag), *iflag);
                tmp[S[i]->id-minId] = -tmp[S[i]->id-minId];
                PHIST_CHK_IERR(SUBR( mvec_add_mvec   ) (tmp[S[i]->id-minId], Vk, st::one(), Vj, iflag), *iflag);
              }

              // store in H
              ST *Hj=NULL;
              phist_lidx ldH;
              PHIST_CHK_IERR(SUBR(sdMat_extract_view)(S[i]->H_,&Hj,&ldH,iflag),*iflag); 
              Hj += (S[i]->curDimV_-1)*ldH;
              Hj[j_] += -tmp[S[i]->id-minId];
            }
          }
        }
      }

      for(int i = 0; i < numSys; i++)
      {
        PHIST_CHK_IERR(SUBR(sdMat_to_device)(S[i]->H_,iflag),*iflag);
      }

      // normalize resulting vector
      // we have already calculated the norm (stored in ynorm)
      for(int i = 0; i < numSys; i++)
      {
        int j = S[i]->curDimV_;
        if( j == 0 )
        {
          // initilize rs_ and normR, normR0
          S[i]->rs_[0] = ynorm[S[i]->id-minId];
          S[i]->normR_ = ynorm[S[i]->id-minId];
          if( S[i]->normR0_ == -mt::one() ) S[i]->normR0_ = S[i]->normR_;
        }
        else
        {
          // raw view of H
          ST *Hj=NULL;
          phist_lidx ldH;
          PHIST_CHK_IERR(SUBR(sdMat_extract_view)(S[i]->H_,&Hj,&ldH,iflag),*iflag); 
          Hj += (S[i]->curDimV_-1)*ldH;
          Hj[j] = ynorm[S[i]->id-minId];
        }
      }
      _ST_ scale[maxId+1-minId];
      for(int i = 0; i < maxId+1-minId; i++)
        scale[i] = st::one() / ynorm[i];
      PHIST_CHK_IERR(SUBR(mvec_vscale)(work.y_, scale, iflag), *iflag);
    }

    maxCurDimV++;
    sharedCurDimV++;
#if defined(TESTING) && (PHIST_OUTLEV>=PHIST_DEBUG)
// check subspace orthogonality
{
  for(int i = 0; i < numSys; i++)
  {
    int nj = S[i]->curDimV_+1;
    _ST_ orth[nj][nj];
    _MT_ maxOrthErr = mt::zero();
    for(int j = 0; j < nj; j++)
    {
      for(int k = 0; k < nj; k++)
      {
        int Vjind = mvecBuff->prevIndex(nextIndex,j);
        int Vkind = mvecBuff->prevIndex(nextIndex,k);
        PHIST_CHK_IERR(SUBR(mvec_view_block)(mvecBuff->at(Vjind), &Vj, S[i]->id, S[i]->id, iflag), *iflag);
        PHIST_CHK_IERR(SUBR(mvec_view_block)(mvecBuff->at(Vkind), &Vk, S[i]->id, S[i]->id, iflag), *iflag);
        PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(Vj,Vk,&orth[j][k],iflag), *iflag);
        if( j == k )
          maxOrthErr = std::max(maxOrthErr, st::abs(orth[j][k]-st::one()));
        else
          maxOrthErr = std::max(maxOrthErr, st::abs(orth[j][k]));
      }
    }
    if( maxOrthErr > 100*mt::eps() )
    {
      PHIST_SOUT(PHIST_INFO,"subspace orthogonality of subspace %d: %8.4e\n", i, maxOrthErr);
    }
  }
}
// check arnoldi/krylov property for last (untransformed row of H): AV_k = V_(k+1) * H_(k+1,k)
{
  TYPE(mvec_ptr) tmpVec = NULL, tmpVec_ = NULL;
  phist_const_map_ptr map;
  PHIST_CHK_IERR(SUBR(mvec_get_map)(work.y_, &map, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_create)(&tmpVec_, map, maxId+1, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_view_block)(tmpVec_, &tmpVec, minId, maxId, iflag), *iflag);
  PHIST_CHK_IERR( Aop->apply (st::one(), Aop->A, work.x_, st::zero(), tmpVec, iflag), *iflag);
  for(int i = 0; i < numSys; i++)
  {
    PHIST_CHK_IERR(SUBR(mvec_view_block)(tmpVec_, &tmpVec, S[i]->id, S[i]->id, iflag), *iflag);
    ST *Hj=NULL;
    phist_lidx ldH; 
    PHIST_CHK_IERR(SUBR(sdMat_extract_view)(S[i]->H_,&Hj,&ldH,iflag),*iflag); 
    Hj += (S[i]->curDimV_-1)*ldH;
    _MT_ maxHerr = mt::zero();
    for(int j = 0; j < S[i]->curDimV_; j++)
    {
      int Vjind = mvecBuff->prevIndex(nextIndex,S[i]->curDimV_-j);
      PHIST_CHK_IERR(SUBR(mvec_view_block)(mvecBuff->at(Vjind), &Vj, S[i]->id, S[i]->id, iflag), *iflag);
      _ST_ Hnj;
      PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(Vj, tmpVec, &Hnj, iflag), *iflag);
      maxHerr = std::max(maxHerr, st::abs(Hj[j]-Hnj));
    }
    if( maxHerr > 100*mt::eps() )
    {
      PHIST_SOUT(PHIST_INFO,"accuracy of last column of H of system %d: %8.4e\n", i, maxHerr);
    }
  }
  PHIST_CHK_IERR(SUBR(mvec_delete)(tmpVec, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_delete)(tmpVec_, iflag), *iflag);
}
#endif

    //    % update QR factorization of H
    for(int i = 0; i < numSys; i++)
    {
      int j = S[i]->curDimV_;
      if( j == 0 )
        continue;

      // raw view of H
      ST *Hj=NULL;
      phist_lidx ldH; 
      PHIST_CHK_IERR(SUBR(sdMat_extract_view)(S[i]->H_,&Hj,&ldH,iflag),*iflag); 
      Hj += (j-1)*ldH;
      // apply previous Gives rotations to column j
      _ST_ tmp;
      for(int k = 0; k < j-1; k++)
      {
        tmp = st::conj(S[i]->cs_[k])*Hj[k] + st::conj(S[i]->sn_[k])*Hj[k+1];
        Hj[k+1] = -S[i]->sn_[k]*Hj[k] + S[i]->cs_[k]*Hj[k+1];
        Hj[k] = tmp;
      }
      // new Givens rotation to eliminate H(j+1,j)
#ifdef IS_COMPLEX
      _MT_ cs;
      PHIST_TG_PREFIX(LARTG)((blas_cmplx*)&Hj[j-1],(blas_cmplx*)&Hj[j],&cs,(blas_cmplx*)&S[i]->sn_[j-1],(blas_cmplx*)&tmp);
      S[i]->cs_[j-1] = (_ST_) cs;
      S[i]->sn_[j-1] = st::conj(S[i]->sn_[j-1]);
#else
      PHIST_TG_PREFIX(LARTG)(&Hj[j-1],&Hj[j],&S[i]->cs_[j-1],&S[i]->sn_[j-1],&tmp);
      //{
        //_MT_ len = mt::sqrt(st::real(st::conj(Hj[j-1])*Hj[j-1])+st::real(st::conj(Hj[j])*Hj[j]));
        //S[i]->cs_[j-1] = Hj[j-1]/len;
        //S[i]->sn_[j-1] = Hj[j]/len;
        //// and apply it
        //tmp = st::conj(S[i]->cs_[j-1])*Hj[j-1] + st::conj(S[i]->sn_[j-1])*Hj[j];
      //}
#endif
#ifdef TESTING
{
  PHIST_OUT(PHIST_DEBUG,"(Hj[j-1],Hj[j]) = (%8.4e+i%8.4e, %8.4e + i%8.4e)\n", st::real(Hj[j-1]), st::imag(Hj[j-1]),st::real(Hj[j]),st::imag(Hj[j]));
  PHIST_OUT(PHIST_DEBUG,"(c,s) = (%8.4e+i%8.4e, %8.4e+i%8.4e)\n", st::real(S[i]->cs_[j-1]),st::imag(S[i]->cs_[j-1]),st::real(S[i]->sn_[j-1]),st::imag(S[i]->sn_[j-1]));
  PHIST_OUT(PHIST_DEBUG,"r = %8.4e + i%8.4e\n", st::real(tmp),st::imag(tmp));
  _ST_ r_ = st::conj(S[i]->cs_[j-1])*Hj[j-1] + st::conj(S[i]->sn_[j-1])*Hj[j];
  _ST_ zero_ = -S[i]->sn_[j-1]*Hj[j-1] + S[i]->cs_[j-1]*Hj[j];
  PHIST_OUT(PHIST_DEBUG,"(r, 0) = (%8.4e + i%8.4e, %8.4e+i%8.4e)\n", st::real(r_), st::imag(r_), st::real(zero_), st::imag(zero_));
  PHIST_CHK_IERR(*iflag = (st::abs(r_-tmp) < 1.e-5) ? 0 : -1, *iflag);
  PHIST_CHK_IERR(*iflag = (st::abs(zero_) < 1.e-5) ? 0 : -1, *iflag);
}
#endif
      // eliminate Hj[j]
      Hj[j-1] = tmp;
      Hj[j] = st::zero();
      // apply to RHS
      tmp = st::conj(S[i]->cs_[j-1])*S[i]->rs_[j-1];
      S[i]->rs_[j] = -S[i]->sn_[j-1]*S[i]->rs_[j-1];
      S[i]->rs_[j-1] = tmp;

      // update current residual norm
      S[i]->normR_=st::abs(S[i]->rs_[j]);
    }

    for(int i = 0; i < numSys; i++)
    {
      // push H to GPU if necessary
      PHIST_CHK_IERR(SUBR(sdMat_to_device)(S[i]->H_,iflag),*iflag);
    }

    //    % check convergence, update subspace dimension etc
    for(int i = 0; i < numSys; i++)
    {
      int j = S[i]->curDimV_;

      if (S[i]->curDimV_) S[i]->totalIter++; // do not count first round after restart twice
      S[i]->curDimV_++;

      // check convergence
      MT relres = S[i]->normR_ / S[i]->normR0_;
      MT absres = S[i]->normR_;
      if( absres < 1000*st::eps() || relres < S[i]->tol )
      {
        S[i]->status = 0; // mark as converged
        anyConverged++;
      }
      else if( S[i]->curDimV_ >= mvecBuff->size() )
      {
        S[i]->status = 2; // mark as failed/restart needed
        anyFull++;
      }
      else
      {
        S[i]->status = 1; // iterating, not converged yet
      }
    }

    // use work.y_ as input in the next iteration
    std::swap(work.x_,work.y_);



    for (int i=0;i<numSys;i++)
    {
      PHIST_SOUT(PHIST_VERBOSE,"[%d]: %d\t%8.4e\t(%8.4e)\n", i, S[i]->curDimV_-1,S[i]->normR_/S[i]->normR0_,S[i]->normR_);
    }
    (*nIter)++;
  }
PHIST_TASK_END(iflag)

  PHIST_SOUT(PHIST_VERBOSE,"%d converged, %d need restart, %d failed.\n",anyConverged,anyFull,anyFailed);
  PHIST_SOUT(PHIST_VERBOSE,"---------------------------------------  \n");

  // delete remaining views (note that our mvec_GC object "work" takes care of x,y and z)
  PHIST_CHK_IERR(SUBR(mvec_delete)(Vj,     iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_delete)(Vk,     iflag), *iflag);

  *iflag=99;
  for (int i=0; i<numSys; i++) *iflag=std::min(*iflag,S[i]->status);
}

