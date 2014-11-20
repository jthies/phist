#include "phist_blockedgmres_helper_def.hpp"

// implementation of minres on several systems simultaneously
void SUBR(blockedMINRESstates_iterate)(TYPE(const_op_ptr) Aop, TYPE(blockedGMRESstate_ptr) S[], int numSys, int* nIter, int* ierr)
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
PHIST_SOUT(PHIST_INFO,"Pipelined MINRES status:\n");
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



  PHIST_SOUT(PHIST_VERBOSE,"MINRES iteration started\n");
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

  for (int i=0;i<numSys;i++)
  {
    PHIST_SOUT(PHIST_VERBOSE,"[%d]: %d\t%8.4e\t(%8.4e)\n", i, S[i]->curDimV_-1,S[i]->normR_/S[i]->normR0_,S[i]->normR_);
  }

  while( anyConverged == 0 && anyFailed == 0 )
  {
    //    % get new vector for y
    int nextIndex;
    PHIST_CHK_IERR( mvecBuff->getNextUnused(nextIndex,ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_view_block ) (mvecBuff->at(nextIndex), &work_y, minId, maxId, ierr), *ierr);

    //    % apply the operator of the matrix A
    PHIST_CHK_IERR( Aop->apply (st::one(), Aop->A, work_x, st::zero(), work_y, ierr), *ierr);


    //    % initialize MINRES for (re-)started systems
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


    //    % lanczos update
    {
      _ST_ prevBeta[maxId+1-minId];
      //_ST_ prevBeta_[maxId+1-minId];
      std::vector<_MT_> beta(maxId+1-minId, -mt::one());
      _ST_ alpha[maxId+1-minId];

      // alpha = work_x' * work_y
      PHIST_CHK_IERR( SUBR(mvec_dot_mvec) (work_x, work_y, alpha, ierr), *ierr);
      for(int i = 0; i < maxId+1-minId; i++)
        alpha[i] = -alpha[i];
      // obtain prevBeta from state
      for(int i = 0; i < numSys; i++)
      {
        if( S[i]->curDimV_ > 1 )
          prevBeta[S[i]->id-minId] = -S[i]->prevBeta_;
        else
          prevBeta[S[i]->id-minId] = st::zero();
      }
      // lanczos: work_y = work_y - beta*v_(k-1) - alpha*work_x
      int prevVind = mvecBuff->prevIndex(S[0]->lastVind_,2);
      PHIST_CHK_IERR( SUBR(mvec_view_block) (mvecBuff->at(prevVind), &Vk, minId, maxId, ierr), *ierr);
      //PHIST_CHK_IERR( SUBR(mvec_dot_mvec) (Vk, work_y, prevBeta_, ierr), *ierr);
      //PHIST_SOUT(PHIST_INFO, "prevBeta (correct val)");
      //for(int i = 0; i < maxId+1-minId; i++)
      //{
        //PHIST_SOUT(PHIST_INFO, "\t%8.4e (%8.4e)", prevBeta[i], prevBeta_[i]);
      //}
      //PHIST_SOUT(PHIST_INFO, "\n");
      PHIST_CHK_IERR( SUBR(mvec_vadd_mvec) (alpha,    work_x, st::one(), work_y, ierr), *ierr);
      PHIST_CHK_IERR( SUBR(mvec_vadd_mvec) (prevBeta, Vk,     st::one(), work_y, ierr), *ierr);

      // calculate new beta
      PHIST_CHK_IERR( SUBR(mvec_norm2) (work_y, &beta[0], ierr), *ierr);
      for(int i = 0; i < numSys; i++)
      {
        S[i]->prevBeta_ = beta[S[i]->id-minId];
        int j = S[i]->curDimV_-1;
        if( j >= 0 )
        {
          // store in H
          ST *Hj=NULL;
          lidx_t ldH;
          PHIST_CHK_IERR(SUBR(sdMat_extract_view)(S[i]->H_,&Hj,&ldH,ierr),*ierr);
          Hj += (S[i]->curDimV_-1)*ldH;
          Hj[j] = -alpha[S[i]->id-minId];

          if( j-1 >= 0 )
            Hj[j-1] = -prevBeta[S[i]->id-minId];
        }
      }

      // normalize
      // we have already calculated the norm (stored in beta)
      for(int i = 0; i < numSys; i++)
      {
        int j = S[i]->curDimV_-1;
        if( j == -1 )
        {
          // initilize rs_
          S[i]->rs_[0] = beta[S[i]->id-minId];
          S[i]->normR_ = beta[S[i]->id-minId];
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
          Hj[j+1] = beta[S[i]->id-minId];
        }
      }
      _ST_ scale[maxId+1-minId];
      for(int i = 0; i < maxId+1-minId; i++)
        scale[i] = st::one() / beta[i];
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



    for (int i=0;i<numSys;i++)
    {
      PHIST_SOUT(PHIST_VERBOSE,"[%d]: %d\t%8.4e\t(%8.4e)\n", i, S[i]->curDimV_-1, S[i]->normR_/S[i]->normR0_, S[i]->normR_);
    }

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

