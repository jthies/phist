// small class for a ring buffer for the subspaces
// TODO: use a real template instead of this "TYPE"-stuff?
class TYPE(MvecRingBuffer)
{
  public:
    TYPE(MvecRingBuffer)(int size) : mvecs_(size,NULL), mvecs_used_(size,0), lastIndex_(0) {}

    // we can handle failures probably more cleanly if not done in the constructor
    void create_mvecs(const_map_ptr_t map, int nvecs, int* ierr)
    {
      for(int i = 0; i < mvecs_.size(); i++)
      {
        PHIST_CHK_IERR(*ierr = (mvecs_[i] != NULL) ? -1 : 0, *ierr);
        PHIST_CHK_IERR(SUBR( mvec_create ) (&mvecs_[i], map, nvecs, ierr), *ierr);
      }
    }

    // must be called before the destructor
    void delete_mvecs(int *ierr)
    {
      for(int i = 0; i < mvecs_.size(); i++)
      {
        PHIST_CHK_IERR(SUBR( mvec_delete ) (mvecs_[i], ierr), *ierr);
        mvecs_[i] = NULL;
      }
    }

    ~TYPE(MvecRingBuffer)()
    {
      int ierr = 0;
      // print an error if there are still allocated mvecs
      // as this could use up all memory quite fast if used multiple times!
      for(int i = 0; i < mvecs_.size(); i++)
      {
        PHIST_CHK_IERR(ierr = (mvecs_[i] != NULL) ? -1 : 0, ierr);
      }
    }

    // get vector at index i
    TYPE(mvec_ptr) at(int i) {return mvecs_.at(i);}

    int size() {return mvecs_.size();}

    // just to make sure no element is used twice
    void incRef(int i) {mvecs_used_.at(i)++;}
    void decRef(int i) {mvecs_used_.at(i)--;}

    // get last used index
    int lastIndex() {return lastIndex_;}

    // get preceeding index
    int prevIndex(int i, int n = 1) {return (i+size()-n)%size();}

    // get next index and make sure it is unused
    void getNextUnused(int &nextIndex, int *ierr)
    {
      *ierr = 0;
      nextIndex = (lastIndex_+1)%mvecs_.size();
      if( mvecs_used_[nextIndex] != 0 )
      {
        *ierr = -1;
        return;
      }
      lastIndex_ = nextIndex;
    }

  private:
    std::vector<TYPE(mvec_ptr)> mvecs_;
    std::vector<int> mvecs_used_;
    int lastIndex_;

    // hide copy constructor etc
    TYPE(MvecRingBuffer)(const TYPE(MvecRingBuffer)&);
    const TYPE(MvecRingBuffer)& operator=(const TYPE(MvecRingBuffer)&);
};


// create new state objects. We just get an array of (NULL-)pointers
void SUBR(pgmresStates_create)(TYPE(pgmresState_ptr) state[], int numSys, const_map_ptr_t map, int maxBas,int* ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
  *ierr=0;

  if (numSys <= 0)
    return;

  const_comm_ptr_t comm;
  PHIST_CHK_IERR(phist_map_get_comm(map,&comm,ierr),*ierr);

  // setup buffer of mvecs to be used later
  Teuchos::RCP<TYPE(MvecRingBuffer)> mvecBuff(new TYPE(MvecRingBuffer)(maxBas+1));
  PHIST_CHK_IERR( mvecBuff->create_mvecs(map, numSys, ierr), *ierr);

  // set up individual states
  for(int i = 0; i < numSys; i++)
  {
    // initialization data for the next state
    TYPE(pgmresState) tmp = {i,(_MT_)0.5,-2,0,-1,0,NULL,NULL,NULL,NULL,NULL,-mt::one(),-mt::one(),NULL};

    // create state
    state[i] = new TYPE(pgmresState)(tmp);

    // allocate members
    PHIST_CHK_IERR(SUBR( sdMat_create )(&state[i]->H_, maxBas+1, maxBas, comm, ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_create  )(&state[i]->b_, map,      1,            ierr), *ierr);
    state[i]->cs_ = new MT[maxBas];
    state[i]->sn_ = new ST[maxBas];
    state[i]->rs_ = new ST[maxBas];

    // assign MvecRingBuffer (with reference counting)
    state[i]->Vbuff = (void*) new Teuchos::RCP<TYPE(MvecRingBuffer)>(mvecBuff);
  }
}


// delete pgmresState object
void SUBR(pgmresStates_delete)(TYPE(pgmresState_ptr) state[], int numSys, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  *ierr=0;

  for(int i = 0; i < numSys; i++)
  {
    PHIST_CHK_IERR(SUBR( sdMat_delete ) (state[i]->H_, ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_delete  ) (state[i]->b_, ierr), *ierr);
    delete [] state[i]->cs_;
    delete [] state[i]->sn_;
    delete [] state[i]->rs_;
    CAST_PTR_FROM_VOID(Teuchos::RCP<TYPE(MvecRingBuffer)>, mvecBuff, state[i]->Vbuff, *ierr);
    PHIST_CHK_IERR((*mvecBuff)->delete_mvecs(ierr), *ierr);
    delete mvecBuff;
    delete state[i];
  }
}


// reset pgmres state.
void SUBR(pgmresState_reset)(TYPE(pgmresState_ptr) S, TYPE(const_mvec_ptr) b, TYPE(const_mvec_ptr) x0, int *ierr)
{
#include "phist_std_typedefs.hpp"  
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  
  if( b == NULL && S->normR0_ == -mt::one() )
  {
    PHIST_OUT(PHIST_ERROR,"on the first call to pgmresState_reset you *must* provide the RHS vector");
    *ierr=-1;
    return;
  }

  if( b != NULL )
  {
    // new rhs -> need to recompute ||b-A*x0||
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(), b, st::zero(), S->b_, ierr), *ierr);
    S->status = -1;
    S->totalIter = 0;
  }


  // get mvecBuff
  CAST_PTR_FROM_VOID(Teuchos::RCP<TYPE(MvecRingBuffer)>, mvecBuffPtr, S->Vbuff, *ierr);
  Teuchos::RCP<TYPE(MvecRingBuffer)> mvecBuff = *mvecBuffPtr;

  if( x0 == NULL )
  {
    // great, we can directly apply a first gmres step as we don't need to compute A*x0

    PHIST_CHK_IERR(SUBR( mvec_norm2 ) (S->b_, &S->normR_, ierr), *ierr);
    if( S->normR0_ < mt::zero() )
      S->normR0_ = S->normR_;
    S->rs_[0] = S->normR_;

    S->lastVind_ = mvecBuff->lastIndex();
    S->curDimV_ = 1;
    TYPE(mvec_ptr) r0 = NULL;
    PHIST_CHK_IERR(SUBR( mvec_view_block )(mvecBuff->at(S->lastVind_), &r0, S->id, S->id, ierr), *ierr);
    mvecBuff->incRef(S->lastVind_);
    if( S->normR_ != 0 )
    {
      _ST_ scale = st::one() / S->normR_;
      PHIST_CHK_IERR(SUBR( mvec_add_mvec ) (scale, S->b_, st::zero(), r0, ierr), *ierr);
    }
    PHIST_CHK_IERR(SUBR( mvec_delete ) (r0, ierr), *ierr);
  }
  else // x != NULL
  {
    // initialize everything to calculate b-A*x0 in the next call to iterate
    S->curDimV_ = 0;
    S->normR0_  = -mt::one(); // needs to be recomputed
    S->lastVind_ = mvecBuff->lastIndex();
    PHIST_CHK_IERR(SUBR( mvec_set_block ) (mvecBuff->at(S->lastVind_), x0, S->id, S->id, ierr), *ierr);
    mvecBuff->incRef(S->lastVind_);
  }
}


// calculate approximate solution
void SUBR(pgmresStates_updateSol)(TYPE(pgmresState_ptr) S_array[], int numSys, TYPE(mvec_ptr) x, _MT_* resNorm, bool scaleSolutionToOne, int* ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
  *ierr = 0;

  if( numSys <= 0 )
    return;

  // if there are multiple systems, make sure their id is correct!
  for(int i = 0; i < numSys; i++)
    PHIST_CHK_IERR(*ierr = (S_array[i]->id != i) ? -1 : 0, *ierr);

  // make sure all lastVind_ are the same
  int lastVind = S_array[0]->lastVind_;
  for(int i = 0; i < numSys; i++)
    PHIST_CHK_IERR(*ierr = (S_array[i]->lastVind_ != lastVind) ? -1 : 0, *ierr);

  CAST_PTR_FROM_VOID(Teuchos::RCP<TYPE(MvecRingBuffer)>, mvecBuffPtr, S_array[0]->Vbuff, *ierr);
  Teuchos::RCP<TYPE(MvecRingBuffer)> mvecBuff = *mvecBuffPtr;

  // make sure all systems use the same mvecBuff
  for(int i = 0; i < numSys; i++)
  {
    CAST_PTR_FROM_VOID(Teuchos::RCP<TYPE(MvecRingBuffer)>, mvecBuffPtr_i, S_array[i]->Vbuff, *ierr);
    PHIST_CHK_IERR(*ierr = (*mvecBuffPtr_i != *mvecBuffPtr) ? -1 : 0, *ierr);
  }

  // determine maximal and shared (minimal) dimensions of subspaces
  int maxCurDimV = 0;
  int sharedCurDimV = mvecBuff->size();
  for(int i = 0; i < numSys; i++)
  {
    maxCurDimV = std::max(maxCurDimV, S_array[i]->curDimV_);
    sharedCurDimV = std::min(sharedCurDimV, S_array[i]->curDimV_);
  }

  // calculate resNorm
  for(int i = 0; i < numSys; i++)
  {
    if( S_array[i]->curDimV_ <= 0 )
      resNorm[i] = -mt::one();
    else if( S_array[i]->curDimV_ == 1 )
      resNorm[i] = mt::one();
    else
      resNorm[i] = S_array[i]->normR_/S_array[i]->normR0_;
  }

  // no iteration done yet?
  if( maxCurDimV <= 1 )
    return;


  // allocate space for y
  _ST_ *yglob = new _ST_[numSys*maxCurDimV];
  for(int i = 0; i < numSys*maxCurDimV; i++)
    yglob[i] = st::zero();
  int ldy = numSys;

  // calculate y by solving the triangular systems
  for(int i = 0; i < numSys; i++)
  {
    // nothing to do here?
    if( S_array[i]->curDimV_ <= 1 )
      continue;

    // helpful variables
    TYPE(const_pgmresState_ptr) S = S_array[i];
    int m = S->curDimV_ - 1;
    ST *H_raw=NULL;
    lidx_t ldH;
    PHIST_CHK_IERR(SUBR(sdMat_extract_view)(S->H_,&H_raw,&ldH,ierr),*ierr);
    _ST_ *y = &yglob[i+ldy*(maxCurDimV-S->curDimV_)];


#if PHIST_OUTLEV>=PHIST_DEBUG
    PHIST_SOUT(PHIST_DEBUG,"pgmres_updateSol[%d], curDimV=%d, H=\n",i,S->curDimV_);
    {
      TYPE(sdMat_ptr) H = NULL;
      PHIST_CHK_IERR(SUBR(sdMat_view_block)(S->H_, &H, 0, m-1, 0, m-1, ierr), *ierr);
      PHIST_CHK_IERR(SUBR(sdMat_print)(H,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(sdMat_delete)(H,ierr),*ierr);
    }
    PHIST_SOUT(PHIST_DEBUG,"rs=\n");
    for (int k=0;k<m;k++)
    {
      PHIST_SOUT(PHIST_DEBUG,"%16.8f+%16.8fi\n",st::real(S->rs_[k]),st::imag(S->rs_[k]));      
    }
#endif

    // set y to rs
    for(int j = 0; j < m; j++)
      y[ldy*j] = S->rs_[j];

#ifdef TESTING
{
  // setup e-vector
  _ST_ e[m];
  e[0] = st::one();
  for(int j = 1; j < m+1; j++)
    e[j] = st::zero();
  // apply givens rotations
  for(int j = 0; j < m-1; j++)
  {
    _ST_ tmp = S->cs_[j] * e[j]  +  S->sn_[j] * e[j+1];
    e[j+1] = -st::conj(S->sn_[j]) * e[j]  +  S->cs_[j] * e[j+1];
    e[j] = tmp;
  }
  // check that y is givens_rotations applied to e
  PHIST_SOUT(PHIST_INFO, "rs/norm0:");
  for(int j = 0; j < m; j++)
    PHIST_SOUT(PHIST_INFO, "\t%8.4e + i%8.4e", st::real(y[ldy*j])/S->normR0_, st::imag(y[ldy*j])/S->normR0_);
  PHIST_SOUT(PHIST_INFO, "\nabs(rs(j)/norm0)):%8.4e", st::abs(y[ldy*(m-1)])/S->normR0_);
  PHIST_SOUT(PHIST_INFO, "\nrot(e_1):");
  for(int j = 0; j < m; j++)
    PHIST_SOUT(PHIST_INFO, "\t%8.4e + i%8.4e", st::real(e[j]), st::imag(e[j]));
  PHIST_SOUT(PHIST_INFO, "\nabs(rot(e_1)):%8.4e\n", st::abs(e[m-1]));
}
#endif

    // solve triangular system
    PHIST_CHK_IERR(PREFIX(TRSV)("U","N","N",&m,(st::blas_scalar_t*)H_raw,&ldH,(st::blas_scalar_t*)y, &ldy, ierr),*ierr);

    // if we are only interested in the directions Vi*yi and appropriate AVi*yi,
    // then this scaling may help to improve the conditioning of a following orthogonalization step!
    if( scaleSolutionToOne )
    {
      // scale y to one
      _MT_ scale = mt::zero();
      for(int j = 0; j < m; j++)
        scale += st::real(st::conj(y[ldy*j])*y[ldy*j]);
      scale = mt::one()/sqrt(scale);
      for(int j = 0; j < m; j++)
        y[ldy*j] *= scale;
    }
  }


  // add up solution
  TYPE(mvec_ptr) Vj = NULL, x_i = NULL;
  for(int j = 1; j <= maxCurDimV; j++)
  {
    int Vind = mvecBuff->prevIndex(S_array[0]->lastVind_,maxCurDimV-j);
    _ST_ *yj = yglob + ldy*(maxCurDimV-j);

    if( j > maxCurDimV-sharedCurDimV )
    {
      // update solution of all systems at once
      PHIST_CHK_IERR(SUBR(mvec_view_block)(mvecBuff->at(Vind), &Vj, 0, numSys-1, ierr), *ierr);
      PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(yj, Vj, st::one(), x, ierr), *ierr);
    }
    else
    {
      // update solution of single systems
      for(int i = 0; i < numSys; i++)
      {
        if( j > maxCurDimV-S_array[i]->curDimV_ )
        {
          PHIST_CHK_IERR(SUBR(mvec_view_block)(mvecBuff->at(Vind), &Vj, S_array[i]->id, S_array[i]->id, ierr), *ierr);
          PHIST_CHK_IERR(SUBR(mvec_view_block)(x, &x_i, i, i, ierr), *ierr);
          PHIST_CHK_IERR(SUBR(mvec_add_mvec)(yj[i], Vj, st::one(), x, ierr), *ierr);
        }
      }
    }
  }

  PHIST_CHK_IERR(SUBR(mvec_delete)(x_i, ierr), *ierr);
  PHIST_CHK_IERR(SUBR(mvec_delete)(Vj, ierr), *ierr);
}


// implementation of gmres on several systems simultaneously
void SUBR(pgmresStates_iterate)(TYPE(const_op_ptr) Aop, TYPE(pgmresState_ptr) S[], int numSys, int* nIter, int* ierr)
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

  // make sure the id of the systems here is correct
// TODO: we could allow to continue the calculation of a single system even for i != id
  for(int i = 0; i < numSys; i++)
    PHIST_CHK_IERR(*ierr = (S[i]->id != i) ? -1 : 0, *ierr);

  CAST_PTR_FROM_VOID(Teuchos::RCP<TYPE(MvecRingBuffer)>, mvecBuffPtr, S[0]->Vbuff, *ierr);
  Teuchos::RCP<TYPE(MvecRingBuffer)> mvecBuff = *mvecBuffPtr;

  // make sure all systems use the same mvecBuff
  for(int i = 0; i < numSys; i++)
  {
    CAST_PTR_FROM_VOID(Teuchos::RCP<TYPE(MvecRingBuffer)>, mvecBuffPtr_i, S[i]->Vbuff, *ierr);
    PHIST_CHK_IERR(*ierr = (*mvecBuffPtr_i != *mvecBuffPtr) ? -1 : 0, *ierr);
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
    PHIST_CHK_IERR(SUBR( mvec_view_block )( mvecBuff->at(lastVind), &work_x, 0, numSys-1, ierr), *ierr);
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

//TODO: check dimensions of all systems!


  PHIST_SOUT(PHIST_VERBOSE,"GMRES iteration started\n");
  PHIST_SOUT(PHIST_VERBOSE,"=======================\n");

  while( anyConverged == 0 && anyFailed == 0 )
  {
    //    % get new vector for y
    int nextIndex;
    PHIST_CHK_IERR( mvecBuff->getNextUnused(nextIndex,ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_view_block ) (mvecBuff->at(nextIndex), &work_y, 0, numSys-1, ierr), *ierr);


    //    % apply the operator of the matrix A
    PHIST_CHK_IERR( Aop->apply (st::one(), Aop->A, work_x, st::zero(), work_y, ierr), *ierr);


    //    % initialize GMRES for (re-)started systems
    for(int i = 0; i < numSys; i++)
    {
      int j = S[i]->curDimV_;
      if( j == 0 )
      {
        // (re-)start: r_0 = b - A*x_0
        PHIST_CHK_IERR( SUBR(mvec_view_block) (work_y, &Vj, i, i, ierr), *ierr);
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


    //    % arnoldi update with modified gram schmidt
    for(int j = 1; j < maxCurDimV; j++)
    {
      int Vind = mvecBuff->prevIndex(S[0]->lastVind_,maxCurDimV-j);
      _ST_ tmp[numSys];

      bool calculatedDot = false;
      if( j > maxCurDimV-sharedCurDimV )
      {
        // MGS step for all systems at once
        PHIST_CHK_IERR(SUBR( mvec_view_block ) (mvecBuff->at(Vind), &Vk, 0, numSys-1, ierr), *ierr);
        PHIST_CHK_IERR(SUBR( mvec_dot_mvec   ) (work_y, Vk, tmp, ierr), *ierr);
        for(int i = 0; i < numSys; i++)
          tmp[i] = -tmp[i];
        PHIST_CHK_IERR(SUBR( mvec_vadd_mvec  ) (tmp, Vk, st::one(), work_y, ierr), *ierr);
        calculatedDot = true;
      }

      // store in H (and do MGS steps for single systems)
      for(int i = 0; i < numSys; i++)
      {
        int j_ = j - (maxCurDimV-S[i]->curDimV_+1);
        if( j_ >= 0 )
        {
          if( !calculatedDot )
          {
            // MGS step for single system
            PHIST_CHK_IERR(SUBR( mvec_view_block ) (mvecBuff->at(Vind), &Vk, i, i, ierr), *ierr);
            PHIST_CHK_IERR(SUBR( mvec_view_block ) (work_y,             &Vj, i, i, ierr), *ierr);
            PHIST_CHK_IERR(SUBR( mvec_dot_mvec   ) (Vj, Vk, &tmp[i], ierr), *ierr);
            tmp[i] = -tmp[i];
            PHIST_CHK_IERR(SUBR( mvec_add_mvec   ) (tmp[i], Vk, st::one(), Vj, ierr), *ierr);
          }

          // store in H
          ST *Hj=NULL;
          lidx_t ldH; 
          PHIST_CHK_IERR(SUBR(sdMat_extract_view)(S[i]->H_,&Hj,&ldH,ierr),*ierr); 
          Hj += (S[i]->curDimV_-1)*ldH;
          Hj[j_] = -tmp[i];
        }
      }
    }
    {
      // normalize resulting vector
      _MT_ tmp[numSys];
      PHIST_CHK_IERR(SUBR(mvec_normalize)(work_y, tmp, ierr), *ierr);
      for(int i = 0; i < numSys; i++)
      {
        int j = S[i]->curDimV_;
        // raw view of H
        ST *Hj=NULL;
        lidx_t ldH; 
        PHIST_CHK_IERR(SUBR(sdMat_extract_view)(S[i]->H_,&Hj,&ldH,ierr),*ierr); 
        Hj += (S[i]->curDimV_-1)*ldH;
        Hj[j] = tmp[i];
      }
    }
    maxCurDimV++;
    sharedCurDimV++;


    //    % update QR factorization of H
    for(int i = 0; i < numSys; i++)
    {
      int j = S[i]->curDimV_;
      if( j == 0 )
      {
        S[i]->rs_[0] = S[i]->normR_;
        continue;
      }

      // raw view of H
      ST *Hj=NULL;
      lidx_t ldH; 
      PHIST_CHK_IERR(SUBR(sdMat_extract_view)(S[i]->H_,&Hj,&ldH,ierr),*ierr); 
      Hj += (j-1)*ldH;
      // apply previous Gives rotations to column j
      _ST_ tmp;
      for(int k = 0; k < j-1; k++)
      {
        tmp = S[i]->cs_[k]*Hj[k] + S[i]->sn_[k]*Hj[k+1];
        Hj[k+1] = -st::conj(S[i]->sn_[k])*Hj[k] + S[i]->cs_[k]*Hj[k+1];
        Hj[k] = tmp;
      }
      // new Givens rotation to eliminate H(j+1,j)
#ifdef IS_COMPLEX
      PREFIX(LARTG)((blas_cmplx_t*)&Hj[j-1],(blas_cmplx_t*)&Hj[j],&S[i]->cs_[j-1],(blas_cmplx_t*)&S[i]->sn_[j-1],(blas_cmplx_t*)&tmp);
#else
      PREFIX(LARTG)(&Hj[j-1],&Hj[j],&S[i]->cs_[j-1],&S[i]->sn_[j-1],&tmp);
#endif
#ifdef TESTING
{
  PHIST_OUT(PHIST_VERBOSE,"(Hj[j-1],Hj[j]) = (%8.4e+i%8.4e, %8.4e + i%8.4e)\n", st::real(Hj[j-1]), st::imag(Hj[j-1]),st::real(Hj[j]),st::imag(Hj[j]));
  PHIST_OUT(PHIST_VERBOSE,"(c,s) = (%8.4e, %8.4e+i%8.4e)\n", S[i]->cs_[j-1],st::real(S[i]->sn_[j-1]),st::imag(S[i]->sn_[j-1]));
  PHIST_OUT(PHIST_VERBOSE,"r = %8.4e + i%8.4e\n", st::real(tmp),st::imag(tmp));
  _ST_ r_ = S[i]->cs_[j-1]*Hj[j-1] + S[i]->sn_[j-1]*Hj[j];
  _ST_ zero_ = -st::conj(S[i]->sn_[j-1])*Hj[j-1] + S[i]->cs_[j-1]*Hj[j];
  PHIST_OUT(PHIST_VERBOSE,"(r, 0) = (%8.4e + i%8.4e, %8.4e+i%8.4e)\n", st::real(r_), st::imag(r_), st::real(zero_), st::imag(zero_));
  PHIST_CHK_IERR(*ierr = (st::abs(r_-tmp) < 1.e-5) ? 0 : -1, *ierr);
  PHIST_CHK_IERR(*ierr = (st::abs(zero_) < 1.e-5) ? 0 : -1, *ierr);
}
#endif
      // eliminate Hj[j]
      Hj[j-1] = tmp;
      Hj[j] = st::zero();
      // apply to RHS
      tmp = S[i]->cs_[j-1]*S[i]->rs_[j-1];
      S[i]->rs_[j] = -st::conj(S[i]->sn_[j-1])*S[i]->rs_[j-1];
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
      else if( S[i]->curDimV_ >= mvecBuff->size()-1 )
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

