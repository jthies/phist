#include "phist_blockedgmres_helper_def.hpp"

#ifdef PHIST_OUTLEV
#undef PHIST_OUTLEV
#endif
#define PHIST_OUTLEV 0

#ifdef TESTING
#undef TESTING
#endif

// create new state objects. We just get an array of (NULL-)pointers
void SUBR(blockedGMRESstates_create)(TYPE(blockedGMRESstate_ptr) state[], int numSys, const_map_ptr_t map, int maxBas,int* ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
  *ierr=0;

  if (numSys <= 0)
    return;

  const_comm_ptr_t comm;
  PHIST_CHK_IERR(phist_map_get_comm(map,&comm,ierr),*ierr);

  // setup buffer of mvecs to be used later
#ifdef PHIST_HAVE_BELOS
  Teuchos::RCP<TYPE(MvecRingBuffer)> mvecBuff(new TYPE(MvecRingBuffer)(maxBas+1));
#else
  TYPE(MvecRingBuffer)* mvecBuff = new TYPE(MvecRingBuffer)(maxBas+1);
#endif
  PHIST_CHK_IERR( mvecBuff->create_mvecs(map, numSys, ierr), *ierr);

  // set up individual states
  for(int i = 0; i < numSys; i++)
  {
    // initialization data for the next state
    TYPE(blockedGMRESstate) tmp = {i,(_MT_)0.5,-2,0,-1,0,NULL,NULL,NULL,NULL,NULL,-mt::one(),-mt::one(),st::zero(),NULL};

    // create state
    state[i] = new TYPE(blockedGMRESstate)(tmp);

    // allocate members
    PHIST_CHK_IERR(SUBR( sdMat_create )(&state[i]->H_, maxBas+1, maxBas, comm, ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_create  )(&state[i]->b_, map,      1,            ierr), *ierr);
    state[i]->cs_ = new ST[maxBas];
    state[i]->sn_ = new ST[maxBas];
    state[i]->rs_ = new ST[maxBas+1];

#ifdef PHIST_HAVE_BELOS
    // assign MvecRingBuffer (with reference counting)
    state[i]->Vbuff = (void*) new Teuchos::RCP<TYPE(MvecRingBuffer)>(mvecBuff);
#else
    // TODO - check memory management, this used to be an RCP
    state[i]->Vbuff = (void*)mvecBuff;
#endif
  }
}


// delete blockedGMRESstate object
void SUBR(blockedGMRESstates_delete)(TYPE(blockedGMRESstate_ptr) state[], int numSys, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
#ifndef PHIST_HAVE_BELOS
  if (numSys==0) return;
  
    CAST_PTR_FROM_VOID(TYPE(MvecRingBuffer), mvecBuff, state[0]->Vbuff, *ierr);    
    PHIST_CHK_IERR(mvecBuff->delete_mvecs(ierr), *ierr);
    
    delete mvecBuff;
#endif

  for(int i = 0; i < numSys; i++)
  {
    PHIST_CHK_IERR(SUBR( sdMat_delete ) (state[i]->H_, ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_delete  ) (state[i]->b_, ierr), *ierr);
    delete [] state[i]->cs_;
    delete [] state[i]->sn_;
    delete [] state[i]->rs_;
#ifdef PHIST_HAVE_BELOS
    CAST_PTR_FROM_VOID(Teuchos::RCP<TYPE(MvecRingBuffer)>, mvecBuff, state[i]->Vbuff, *ierr);
    PHIST_CHK_IERR((*mvecBuff)->delete_mvecs(ierr), *ierr);
    delete mvecBuff;
#endif
    delete state[i];
  }
}


// reset blockedGMRES state.
void SUBR(blockedGMRESstate_reset)(TYPE(blockedGMRESstate_ptr) S, TYPE(const_mvec_ptr) b, TYPE(const_mvec_ptr) x0, int *ierr)
{
#include "phist_std_typedefs.hpp"  
  ENTER_FCN(__FUNCTION__);
  *ierr=0;

  // get mvecBuff
#ifdef PHIST_HAVE_BELOS
  CAST_PTR_FROM_VOID(Teuchos::RCP<TYPE(MvecRingBuffer)>, mvecBuffPtr, S->Vbuff, *ierr);
  Teuchos::RCP<TYPE(MvecRingBuffer)> mvecBuff = *mvecBuffPtr;
#else
  CAST_PTR_FROM_VOID(TYPE(MvecRingBuffer), mvecBuff, S->Vbuff, *ierr);
#endif

  // release mvecs currently marked used by this state
  for(int j = 0; j < S->curDimV_; j++)
  {
    int Vind = mvecBuff->prevIndex(S->lastVind_,j);
    mvecBuff->decRef(Vind);
  }
  S->curDimV_ = 0;

  // only freed resources
  if( b == NULL && x0 == NULL )
  {
    S->status = -2;
    return;
  }

  if( b == NULL && (S->normR0_ == -mt::one() || S->status == -2) )
  {
    PHIST_OUT(PHIST_ERROR,"on the first call to blockedGMRESstate_reset you *must* provide the RHS vector");
    *ierr=-1;
    return;
  }

  if( b != NULL )
  {
    // new rhs -> need to recompute ||b-A*x0||
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(), b, st::zero(), S->b_, ierr), *ierr);
    S->status = -1;
    S->totalIter = 0;
    S->normR0_ = -mt::one();
  }

  // set H to zero
  PHIST_CHK_IERR(SUBR(sdMat_put_value)(S->H_, st::zero(), ierr), *ierr);

  if( x0 == NULL )
  {
    // great, we can directly apply a first gmres step as we don't need to compute A*x0

    PHIST_CHK_IERR(SUBR( mvec_norm2 ) (S->b_, &S->normR_, ierr), *ierr);
    if( S->normR0_ < mt::zero() )
      S->normR0_ = S->normR_;
    S->rs_[0] = S->normR_;
    S->prevBeta_ = S->normR_;

    S->lastVind_ = mvecBuff->lastIndex();
    S->curDimV_ = 1;
    TYPE(mvec_ptr) r0 = NULL;
    PHIST_CHK_IERR(SUBR( mvec_view_block )(mvecBuff->at(S->lastVind_), &r0, S->id, S->id, ierr), *ierr);
    mvecBuff->incRef(S->lastVind_);
    if( S->normR_ != mt::zero() )
    {
      _ST_ scale = st::one() / S->normR_;
      PHIST_CHK_IERR(SUBR( mvec_add_mvec ) (scale, S->b_, st::zero(), r0, ierr), *ierr);
    }
    PHIST_CHK_IERR(SUBR( mvec_delete ) (r0, ierr), *ierr);
  }
  else // x != NULL
  {
    // initialize everything to calculate b-A*x0 in the next call to iterate
    S->lastVind_ = mvecBuff->lastIndex();
    PHIST_CHK_IERR(SUBR( mvec_set_block ) (mvecBuff->at(S->lastVind_), x0, S->id, S->id, ierr), *ierr);
    mvecBuff->incRef(S->lastVind_);
    S->prevBeta_ = st::zero();
  }

  // update status
  if( S->status >= 0 )
    S->status = 1;
}


// calculate approximate solution
void SUBR(blockedGMRESstates_updateSol)(TYPE(blockedGMRESstate_ptr) S[], int numSys, TYPE(mvec_ptr) x, _MT_* resNorm, bool scaleSolutionToOne, int* ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
  *ierr = 0;

  if( numSys <= 0 )
    return;

  // if there are multiple systems, get the maximal id
  int maxId = 0;
  for(int i = 0; i < numSys; i++)
    maxId = std::max(maxId,S[i]->id);

  bool ordered = true;
  for(int i = 0; i < numSys; i++)
    ordered = ordered && (S[i]->id == i);

  // make sure all lastVind_ are the same
  int lastVind = S[0]->lastVind_;
  for(int i = 0; i < numSys; i++)
    PHIST_CHK_IERR(*ierr = (S[i]->lastVind_ != lastVind) ? -1 : 0, *ierr);

#ifdef PHIST_HAVE_BELOS
  CAST_PTR_FROM_VOID(Teuchos::RCP<TYPE(MvecRingBuffer)>, mvecBuffPtr, S[0]->Vbuff, *ierr);
  Teuchos::RCP<TYPE(MvecRingBuffer)> mvecBuff = *mvecBuffPtr;
#else
  CAST_PTR_FROM_VOID(TYPE(MvecRingBuffer), mvecBuff, S[0]->Vbuff, *ierr);
#endif

  // make sure all systems use the same mvecBuff
  for(int i = 0; i < numSys; i++)
  {
#ifdef PHIST_HAVE_BELOS
    CAST_PTR_FROM_VOID(Teuchos::RCP<TYPE(MvecRingBuffer)>, mvecBuffPtr_i, S[i]->Vbuff, *ierr);
    PHIST_CHK_IERR(*ierr = (*mvecBuffPtr_i != *mvecBuffPtr) ? -1 : 0, *ierr);
#else
    CAST_PTR_FROM_VOID(TYPE(MvecRingBuffer), mvecBuffPtr_i, S[i]->Vbuff, *ierr);
    PHIST_CHK_IERR(*ierr = (mvecBuffPtr_i != mvecBuff) ? -1 : 0, *ierr);
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
  if( maxCurDimV <= 1 )
    return;


  // allocate space for y
  _ST_ *yglob = new _ST_[(maxId+1)*(maxCurDimV-1)];
  for(int i = 0; i < (maxId+1)*(maxCurDimV-1); i++)
    yglob[i] = st::zero();
  lidx_t ldy = (maxId+1);

  // calculate y by solving the triangular systems
  for(int i = 0; i < numSys; i++)
  {
    // nothing to do here?
    if( S[i]->curDimV_ <= 1 )
      continue;

    // helpful variables
    int m = S[i]->curDimV_-1;
    ST *H_raw=NULL;
    lidx_t ldH;
    PHIST_CHK_IERR(SUBR(sdMat_extract_view)(S[i]->H_,&H_raw,&ldH,ierr),*ierr);
    _ST_ *y = &yglob[S[i]->id+ldy*(maxCurDimV-S[i]->curDimV_)];


#if 0
//PHIST_OUTLEV>=PHIST_DEBUG
    PHIST_SOUT(PHIST_DEBUG,"blockedGMRES_updateSol[%d], curDimV=%d, H=\n",i,S[i]->curDimV_);
    {
      TYPE(sdMat_ptr) H = NULL;
      PHIST_CHK_IERR(SUBR(sdMat_view_block)(S[i]->H_, &H, 0, m-1, 0, m-1, ierr), *ierr);
      PHIST_CHK_IERR(SUBR(sdMat_print)(H,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(sdMat_delete)(H,ierr),*ierr);
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


#ifdef TESTING
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
    blas_idx_t ildH=static_cast<blas_idx_t>(ldH);
    blas_idx_t ildy=static_cast<blas_idx_t>(ldy);
    if (ildH<0 || ildy<0)
    {
      *ierr=PHIST_INTEGER_OVERFLOW;
      return;
    }
    PHIST_CHK_IERR(PREFIX(TRSV)("U","N","N",&m,(st::blas_scalar_t*)H_raw,&ildH,(st::blas_scalar_t*)y, &ildy, ierr),*ierr);


    // if we are only interested in the directions Vi*yi and appropriate AVi*yi,
    // then this scaling may help to improve the conditioning of a following orthogonalization step!
    if( scaleSolutionToOne )
    {
      // scale y to one
      _MT_ scale = mt::zero();
      for(int j = 0; j < m; j++)
        scale += st::real(st::conj(y[ldy*j])*y[ldy*j]);
      scale = mt::one()/sqrt(scale);
      PHIST_SOUT(PHIST_DEBUG,"scaling solution with: %8.4e\n", scale);
      for(int j = 0; j < m; j++)
        y[ldy*j] *= scale;
    }
  }


  // add up solution
  TYPE(mvec_ptr) Vj = NULL, x_i = NULL;
  for(int j = 0; j < maxCurDimV-1; j++)
  {
    int Vind = mvecBuff->prevIndex(lastVind,maxCurDimV-1-j);
    _ST_ *yj = yglob + ldy*j;
//std::cout << "j " << j << " yj " << *yj << " maxCurDimV " << maxCurDimV << " sharedCurDimV " << sharedCurDimV << " Vind " << Vind << std::endl;

    if( j >= maxCurDimV-sharedCurDimV && ordered )
    {
      // update solution of all systems at once
      PHIST_CHK_IERR(SUBR(mvec_view_block)(mvecBuff->at(Vind), &Vj, 0, maxId, ierr), *ierr);
      PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(yj, Vj, st::one(), x, ierr), *ierr);
    }
    else
    {
      // update solution of single systems
      for(int i = 0; i < numSys; i++)
      {
        if( j >= maxCurDimV-S[i]->curDimV_ )
        {
          PHIST_CHK_IERR(SUBR(mvec_view_block)(mvecBuff->at(Vind), &Vj, S[i]->id, S[i]->id, ierr), *ierr);
          PHIST_CHK_IERR(SUBR(mvec_view_block)(x, &x_i, i, i, ierr), *ierr);
          PHIST_CHK_IERR(SUBR(mvec_add_mvec)(yj[S[i]->id], Vj, st::one(), x, ierr), *ierr);
        }
      }
    }
  }

  PHIST_CHK_IERR(SUBR(mvec_delete)(x_i, ierr), *ierr);
  PHIST_CHK_IERR(SUBR(mvec_delete)(Vj, ierr), *ierr);
}


// implementation of gmres on several systems simultaneously
void SUBR(blockedGMRESstates_iterate)(TYPE(const_op_ptr) Aop, TYPE(blockedGMRESstate_ptr) S[], int numSys, int* nIter, bool useIMGS, int* ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
  *ierr = 0;

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
  for(int i = 0; i < numSys; i++)
    maxId = std::max(maxId,S[i]->id);
  int minId = maxId;
  for(int i = 0; i < numSys; i++)
    minId = std::min(minId,S[i]->id);

#ifdef PHIST_HAVE_BELOS
  CAST_PTR_FROM_VOID(Teuchos::RCP<TYPE(MvecRingBuffer)>, mvecBuffPtr, S[0]->Vbuff, *ierr);
  Teuchos::RCP<TYPE(MvecRingBuffer)> mvecBuff = *mvecBuffPtr;
#else
  CAST_PTR_FROM_VOID(TYPE(MvecRingBuffer), mvecBuff, S[0]->Vbuff, *ierr);
#endif

  // make sure all systems use the same mvecBuff
  for(int i = 0; i < numSys; i++)
  {
#ifdef PHIST_HAVE_BELOS
    CAST_PTR_FROM_VOID(Teuchos::RCP<TYPE(MvecRingBuffer)>, mvecBuffPtr_i, S[i]->Vbuff, *ierr);
    PHIST_CHK_IERR(*ierr = (*mvecBuffPtr_i != *mvecBuffPtr) ? -1 : 0, *ierr);
#else
    CAST_PTR_FROM_VOID(TYPE(MvecRingBuffer), mvecBuffPtr_i, S[i]->Vbuff, *ierr);
    PHIST_CHK_IERR(*ierr = (mvecBuffPtr_i != mvecBuff) ? -1 : 0, *ierr);
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

  // work vector for x and y = Aop(x)
  TYPE(mvec_ptr) work_x = NULL;
  TYPE(mvec_ptr) work_y = NULL;

  {
    // make sure all lastVind_ are the same
    int lastVind = S[0]->lastVind_;
    for(int i = 0; i < numSys; i++)
      PHIST_CHK_IERR(*ierr = (S[i]->lastVind_ != lastVind) ? -1 : 0, *ierr);

    // x0 / last element of krylov subspace
    PHIST_CHK_IERR(SUBR( mvec_view_block )( mvecBuff->at(lastVind), &work_x, minId, maxId, ierr), *ierr);

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
PHIST_SOUT(PHIST_INFO,"Pipelined GMRES status:\n");
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
  int anyFailed = 0;



  PHIST_SOUT(PHIST_VERBOSE,"GMRES iteration started\n");
  PHIST_SOUT(PHIST_VERBOSE,"=======================\n");

  // check if there are already converged/failed systems
  for(int i = 0; i < numSys; i++)
  {
    // check convergence
    MT relres = S[i]->normR_ / S[i]->normR0_;
    MT absres = S[i]->normR_;
    if( S[i]->normR0_ != -mt::one() && ( absres < 100*st::eps() || relres < S[i]->tol ) )
    {
      S[i]->status = 0; // mark as converged
      anyConverged++;
    }
    else if( S[i]->curDimV_ >= mvecBuff->size() )
    {
      S[i]->status = 2; // mark as failed/restart needed
      anyFailed++;
    }
  }

  while( anyConverged == 0 && anyFailed == 0 )
  {
    //    % get new vector for y
    int nextIndex;
    PHIST_CHK_IERR( mvecBuff->getNextUnused(nextIndex,ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_view_block ) (mvecBuff->at(nextIndex), &work_y, minId, maxId, ierr), *ierr);

    //    % apply the operator of the matrix A
    PHIST_CHK_IERR( Aop->apply (st::one(), Aop->A, work_x, st::zero(), work_y, ierr), *ierr);


    //    % initialize GMRES for (re-)started systems
    for(int i = 0; i < numSys; i++)
    {
      int j = S[i]->curDimV_;
      if( j == 0 )
      {
        // (re-)start: r_0 = b - A*x_0
        PHIST_CHK_IERR( SUBR(mvec_view_block) (work_y, &Vj, S[i]->id-minId, S[i]->id-minId, ierr), *ierr);
        PHIST_CHK_IERR( SUBR(mvec_add_mvec) (st::one(), S[i]->b_, -st::one(), Vj, ierr), *ierr);
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


    //    % arnoldi update with iterated modified gram schmidt
    {
      std::vector<_MT_> ynorm(maxId+1-minId,-mt::one());
      std::vector<_MT_> prev_ynorm(maxId+1-minId,-mt::one());
      for(int mgsIter = 0; mgsIter < 3; mgsIter++)
      {
        // calculate norm
        prev_ynorm = ynorm;
        if( useIMGS || mgsIter == 1 )
        {
          PHIST_CHK_IERR(SUBR(mvec_norm2)(work_y, &ynorm[0], ierr), *ierr);
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
          PHIST_SOUT(PHIST_INFO, "Additional MGS iteration in blockedGMRES!\n");
        }


        for(int j = 0; j < maxCurDimV; j++)
        {
          int Vind = mvecBuff->prevIndex(S[0]->lastVind_,maxCurDimV-j);
          _ST_ tmp[maxId+1-minId];

          bool calculatedDot = false;
          if( j >= maxCurDimV-sharedCurDimV && maxId+1-minId == numSys )
          {
            // MGS step for all systems at once
            PHIST_CHK_IERR(SUBR( mvec_view_block ) (mvecBuff->at(Vind), &Vk, minId, maxId, ierr), *ierr);
            PHIST_CHK_IERR(SUBR( mvec_dot_mvec   ) (Vk, work_y, tmp, ierr), *ierr);
            for(int i = 0; i < maxId+1-minId; i++)
              tmp[i] = -tmp[i];
            PHIST_CHK_IERR(SUBR( mvec_vadd_mvec  ) (tmp, Vk, st::one(), work_y, ierr), *ierr);
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
                PHIST_CHK_IERR(SUBR( mvec_view_block ) (mvecBuff->at(Vind), &Vk, S[i]->id,       S[i]->id,       ierr), *ierr);
                PHIST_CHK_IERR(SUBR( mvec_view_block ) (work_y,             &Vj, S[i]->id-minId, S[i]->id-minId, ierr), *ierr);
                PHIST_CHK_IERR(SUBR( mvec_dot_mvec   ) (Vk, Vj, &tmp[S[i]->id-minId], ierr), *ierr);
                tmp[S[i]->id-minId] = -tmp[S[i]->id-minId];
                PHIST_CHK_IERR(SUBR( mvec_add_mvec   ) (tmp[S[i]->id-minId], Vk, st::one(), Vj, ierr), *ierr);
              }

              // store in H
              ST *Hj=NULL;
              lidx_t ldH; 
              PHIST_CHK_IERR(SUBR(sdMat_extract_view)(S[i]->H_,&Hj,&ldH,ierr),*ierr); 
              Hj += (S[i]->curDimV_-1)*ldH;
              Hj[j_] += -tmp[S[i]->id-minId];
            }
          }
        }
      }

      // normalize resulting vector
      // we have already calculated the norm (stored in ynorm)
      for(int i = 0; i < numSys; i++)
      {
        int j = S[i]->curDimV_;
        if( j == 0 )
        {
          // initilize rs_
          S[i]->rs_[0] = ynorm[S[i]->id-minId];
          S[i]->normR_ = ynorm[S[i]->id-minId];
          if( S[i]->normR0_ == -mt::one() )
            S[i]->normR0_ = S[i]->normR_;
        }
        else
        {
          // raw view of H
          ST *Hj=NULL;
          lidx_t ldH; 
          PHIST_CHK_IERR(SUBR(sdMat_extract_view)(S[i]->H_,&Hj,&ldH,ierr),*ierr); 
          Hj += (S[i]->curDimV_-1)*ldH;
          Hj[j] = ynorm[S[i]->id-minId];
        }
      }
      _ST_ scale[maxId+1-minId];
      for(int i = 0; i < maxId+1-minId; i++)
        scale[i] = st::one() / ynorm[i];
      PHIST_CHK_IERR(SUBR(mvec_vscale)(work_y, scale, ierr), *ierr);
    }
    maxCurDimV++;
    sharedCurDimV++;
#ifdef TESTING
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
        PHIST_CHK_IERR(SUBR(mvec_view_block)(mvecBuff->at(Vjind), &Vj, S[i]->id, S[i]->id, ierr), *ierr);
        PHIST_CHK_IERR(SUBR(mvec_view_block)(mvecBuff->at(Vkind), &Vk, S[i]->id, S[i]->id, ierr), *ierr);
        PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(Vj,Vk,&orth[j][k],ierr), *ierr);
        if( j == k )
          maxOrthErr = std::max(maxOrthErr, st::abs(orth[j][k]-st::one()));
        else
          maxOrthErr = std::max(maxOrthErr, st::abs(orth[j][k]));
      }
    }
    PHIST_SOUT(PHIST_INFO,"subspace orthogonality of subspace %d: %8.4e\n", i, maxOrthErr);
    //for(int j = 0; j < nj; j++)
    //{
      //for(int k = 0; k < nj; k++)
      //{
        //PHIST_SOUT(PHIST_INFO,"\t%8.4e", st::abs(orth[j][k]-( (j==k) ? st::one() : st::zero() )) );
      //}
      //PHIST_SOUT(PHIST_INFO,"\n");
    //}
  }
}
// check arnoldi/krylov property for last (untransformed row of H): AV_k = V_(k+1) * H_(k+1,k)
{
  TYPE(mvec_ptr) tmpVec = NULL, tmpVec_ = NULL;
  const_map_ptr_t map;
  PHIST_CHK_IERR(SUBR(mvec_get_map)(work_y, &map, ierr), *ierr);
  PHIST_CHK_IERR(SUBR(mvec_create)(&tmpVec_, map, maxId+1, ierr), *ierr);
  PHIST_CHK_IERR(SUBR(mvec_view_block)(tmpVec_, &tmpVec, minId, maxId, ierr), *ierr);
  PHIST_CHK_IERR( Aop->apply (st::one(), Aop->A, work_x, st::zero(), tmpVec, ierr), *ierr);
  for(int i = 0; i < numSys; i++)
  {
    PHIST_CHK_IERR(SUBR(mvec_view_block)(tmpVec_, &tmpVec, S[i]->id, S[i]->id, ierr), *ierr);
    ST *Hj=NULL;
    lidx_t ldH; 
    PHIST_CHK_IERR(SUBR(sdMat_extract_view)(S[i]->H_,&Hj,&ldH,ierr),*ierr); 
    Hj += (S[i]->curDimV_-1)*ldH;
    PHIST_SOUT(PHIST_INFO,"accuracy of last column of H of system %d:\n", i);
    for(int j = 0; j < S[i]->curDimV_; j++)
    {
      int Vjind = mvecBuff->prevIndex(nextIndex,S[i]->curDimV_-j);
      PHIST_CHK_IERR(SUBR(mvec_view_block)(mvecBuff->at(Vjind), &Vj, S[i]->id, S[i]->id, ierr), *ierr);
      _ST_ Hnj;
      PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(Vj, tmpVec, &Hnj, ierr), *ierr);
      PHIST_SOUT(PHIST_INFO,"\t%8.4e", st::abs(Hj[j]-Hnj));
    }
    PHIST_SOUT(PHIST_INFO,"\n");
  }
  PHIST_CHK_IERR(SUBR(mvec_delete)(tmpVec, ierr), *ierr);
  PHIST_CHK_IERR(SUBR(mvec_delete)(tmpVec_, ierr), *ierr);
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
      lidx_t ldH; 
      PHIST_CHK_IERR(SUBR(sdMat_extract_view)(S[i]->H_,&Hj,&ldH,ierr),*ierr); 
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
      PREFIX(LARTG)((blas_cmplx_t*)&Hj[j-1],(blas_cmplx_t*)&Hj[j],&cs,(blas_cmplx_t*)&S[i]->sn_[j-1],(blas_cmplx_t*)&tmp);
      S[i]->cs_[j-1] = (_ST_) cs;
      S[i]->sn_[j-1] = st::conj(S[i]->sn_[j-1]);
#else
      PREFIX(LARTG)(&Hj[j-1],&Hj[j],&S[i]->cs_[j-1],&S[i]->sn_[j-1],&tmp);
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
  PHIST_OUT(PHIST_VERBOSE,"(Hj[j-1],Hj[j]) = (%8.4e+i%8.4e, %8.4e + i%8.4e)\n", st::real(Hj[j-1]), st::imag(Hj[j-1]),st::real(Hj[j]),st::imag(Hj[j]));
  PHIST_OUT(PHIST_VERBOSE,"(c,s) = (%8.4e+i%8.4e, %8.4e+i%8.4e)\n", st::real(S[i]->cs_[j-1]),st::imag(S[i]->cs_[j-1]),st::real(S[i]->sn_[j-1]),st::imag(S[i]->sn_[j-1]));
  PHIST_OUT(PHIST_VERBOSE,"r = %8.4e + i%8.4e\n", st::real(tmp),st::imag(tmp));
  _ST_ r_ = st::conj(S[i]->cs_[j-1])*Hj[j-1] + st::conj(S[i]->sn_[j-1])*Hj[j];
  _ST_ zero_ = -S[i]->sn_[j-1]*Hj[j-1] + S[i]->cs_[j-1]*Hj[j];
  PHIST_OUT(PHIST_VERBOSE,"(r, 0) = (%8.4e + i%8.4e, %8.4e+i%8.4e)\n", st::real(r_), st::imag(r_), st::real(zero_), st::imag(zero_));
  PHIST_CHK_IERR(*ierr = (st::abs(r_-tmp) < 1.e-5) ? 0 : -1, *ierr);
  PHIST_CHK_IERR(*ierr = (st::abs(zero_) < 1.e-5) ? 0 : -1, *ierr);
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


    //    % check convergence, update subspace dimension etc
    for(int i = 0; i < numSys; i++)
    {
      int j = S[i]->curDimV_;

      S[i]->curDimV_++;
      S[i]->totalIter++;

      // check convergence
      MT relres = S[i]->normR_ / S[i]->normR0_;
      MT absres = S[i]->normR_;
      if( absres < 100*st::eps() || relres < S[i]->tol )
      {
        S[i]->status = 0; // mark as converged
        anyConverged++;
      }
      else if( S[i]->curDimV_ >= mvecBuff->size() )
      {
        S[i]->status = 2; // mark as failed/restart needed
        anyFailed++;
      }
      else
      {
        S[i]->status = 1; // iterating, not converged yet
      }
    }

    // use work_y as input in the next iteration
    std::swap(work_x,work_y);



#if PHIST_OUTLEV>=PHIST_VERBOSE
    for (int i=0;i<numSys;i++)
    {
    PHIST_SOUT(PHIST_VERBOSE,"[%d]: %d\t%8.4e\t(%8.4e)\n",i,
          S[i]->curDimV_-1,S[i]->normR_/S[i]->normR0_,S[i]->normR_);
    }
#endif

    (*nIter)++;
  }

  PHIST_SOUT(PHIST_VERBOSE,"%d converged, %d failed.\n",anyConverged,anyFailed);
  PHIST_SOUT(PHIST_VERBOSE,"-----------------------\n");

  // delete views
  PHIST_CHK_IERR(SUBR(mvec_delete)(work_x, ierr), *ierr);
  PHIST_CHK_IERR(SUBR(mvec_delete)(work_y, ierr), *ierr);
  PHIST_CHK_IERR(SUBR(mvec_delete)(Vj,     ierr), *ierr);
  PHIST_CHK_IERR(SUBR(mvec_delete)(Vk,     ierr), *ierr);

  if (anyConverged > 0)
    *ierr=0;
      
  if (anyFailed > 0)
    *ierr=1;
}
