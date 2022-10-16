/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_blockedgmres_helper_def.hpp"

// implementation of minres on several systems simultaneously
void SUBR(blockedMINRESstates_iterate)(TYPE(const_linearOp_ptr) Aop, 
                                       TYPE(const_linearOp_ptr) rightPrecon,
                                       TYPE(blockedGMRESstate_ptr) S[], int numSys, int *nIter, int* iflag)
{
#include "phist_std_typedefs.hpp"
  *iflag = 0;
  if (numSys==0) return; // do not appear in timing stats
  PHIST_ENTER_FCN(__FUNCTION__);

  int maxIter=(*nIter)>0?*nIter:9999999;
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

  if (rightPrecon!=NULL)
  {
    PHIST_SOUT(PHIST_WARNING,"preconditioning not implemented in %s\n",__FUNCTION__);
    
  }

  // if there are multiple systems, get the maximal id
  int maxId = 0;
  for(int i = 0; i < numSys; i++)
    maxId = std::max(maxId,S[i]->id);
  int minId = maxId;
  for(int i = 0; i < numSys; i++)
    minId = std::min(minId,S[i]->id);

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

  // work vector for x and y = Aop(x)
  TYPE(mvec_ptr) work_x = NULL;
  TYPE(mvec_ptr) work_y = NULL;

  {
    // make sure all lastVind_ are the same
    int lastVind = S[0]->lastVind_;
    for(int i = 0; i < numSys; i++)
      PHIST_CHK_IERR(*iflag = (S[i]->lastVind_ != lastVind) ? -1 : 0, *iflag);

    // x0 / last element of krylov subspace
    PHIST_CHK_IERR(SUBR( mvec_view_block )( mvecBuff->at(lastVind), &work_x, minId, maxId, iflag), *iflag);

#ifdef PHIST_TESTING
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
  int anyFull      = 0;
  int anyFailed    = 0;



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
    else if( S[i]->totalIter >=maxIter )
    {
      S[i]->status = 3; // mark as failed
      anyFull++;
    }
    else if( S[i]->curDimV_ >= mvecBuff->size() )
    {
      S[i]->status = 2; // mark as restart needed
      anyFull++;
    }
  }

  for (int i=0;i<numSys;i++)
  {
    PHIST_SOUT(PHIST_VERBOSE,"[%d]: %d\t%8.4e\t(%8.4e)\n", i, S[i]->curDimV_-1,S[i]->normR_/S[i]->normR0_,S[i]->normR_);
  }

// put all iterations in one big compute task; this speeds up the tests with ghost (significantly)
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
  while( anyConverged==0 && anyFull==0 && anyFailed==0 )
  {
    //    % get new vector for y
    int nextIndex;
    PHIST_CHK_IERR( mvecBuff->getNextUnused(nextIndex,iflag), *iflag);
    PHIST_CHK_IERR(SUBR( mvec_view_block ) (mvecBuff->at(nextIndex), &work_y, minId, maxId, iflag), *iflag);

    //    % apply the operator of the matrix A
    PHIST_CHK_IERR( Aop->apply (st::one(), Aop->A, work_x, st::zero(), work_y, iflag), *iflag);


    //    % initialize MINRES for (re-)started systems
    for(int i = 0; i < numSys; i++)
    {
      int j = S[i]->curDimV_;
      if( j == 0 )
      {
        // (re-)start: r_0 = b - A*x_0
        PHIST_CHK_IERR( SUBR(mvec_view_block) (work_y, &Vj, S[i]->id-minId, S[i]->id-minId, iflag), *iflag);
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
      PHIST_CHK_IERR(SUBR(sdMat_from_device)(S[i]->H_,iflag),*iflag);
    }

    //    % lanczos update
    {
      _ST_ prevBeta[maxId+1-minId];
      std::vector<_MT_> beta(maxId+1-minId, -mt::one());
      _ST_ alpha[maxId+1-minId];

      // alpha = work_x' * work_y
      if( sharedCurDimV > 0 )
      {
        PHIST_CHK_IERR( SUBR(mvec_dot_mvec) (work_x, work_y, alpha, iflag), *iflag);
      }
      for(int i = 0; i < numSys; i++)
      {
        if( S[i]->curDimV_ > 0 )
          alpha[S[i]->id-minId] = -alpha[S[i]->id-minId];
        else
          alpha[S[i]->id-minId] = st::zero();
      }

      // obtain prevBeta from state
      for(int i = 0; i < numSys; i++)
      {
        if( S[i]->curDimV_ > 1 )
          prevBeta[S[i]->id-minId] = -S[i]->prevBeta_;
        else
          prevBeta[S[i]->id-minId] = st::zero();
      }

      // lanczos: work_y = work_y - beta*v_(k-1) - alpha*work_x
      if(sharedCurDimV > 0)
      {
        PHIST_CHK_IERR( SUBR(mvec_vadd_mvec) (alpha,    work_x, st::one(), work_y, iflag), *iflag);
      }
      if( sharedCurDimV > 1)
      {
        int prevVind = mvecBuff->prevIndex(S[0]->lastVind_,2);
        PHIST_CHK_IERR( SUBR(mvec_view_block) (mvecBuff->at(prevVind), &Vk, minId, maxId, iflag), *iflag);
        PHIST_CHK_IERR( SUBR(mvec_vadd_mvec) (prevBeta, Vk,     st::one(), work_y, iflag), *iflag);
      }

      // calculate new beta
      PHIST_CHK_IERR( SUBR(mvec_norm2) (work_y, &beta[0], iflag), *iflag);
      for(int i = 0; i < numSys; i++)
      {
        S[i]->prevBeta_ = beta[S[i]->id-minId];
        int j = S[i]->curDimV_-1;
        if( j >= 0 )
        {
          // store in H
          ST *Hj=NULL;
          phist_lidx ldH;
          PHIST_CHK_IERR(SUBR(sdMat_extract_view)(S[i]->H_,&Hj,&ldH,iflag),*iflag);
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
          phist_lidx ldH; 
          PHIST_CHK_IERR(SUBR(sdMat_extract_view)(S[i]->H_,&Hj,&ldH,iflag),*iflag); 
          Hj += (S[i]->curDimV_-1)*ldH;
          Hj[j+1] = beta[S[i]->id-minId];
        }
      }
      _ST_ scale[maxId+1-minId];
      for(int i = 0; i < maxId+1-minId; i++)
        scale[i] = st::one() / beta[i];
      PHIST_CHK_IERR(SUBR(mvec_vscale)(work_y, scale, iflag), *iflag);
    }
    maxCurDimV++;
    sharedCurDimV++;
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_DEBUG)
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
  PHIST_CHK_IERR(SUBR(mvec_get_map)(work_y, &map, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_create)(&tmpVec_, map, maxId+1, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_view_block)(tmpVec_, &tmpVec, minId, maxId, iflag), *iflag);
  PHIST_CHK_IERR( Aop->apply (st::one(), Aop->A, work_x, st::zero(), tmpVec, iflag), *iflag);
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
      maxHerr = std::max(maxHerr,st::abs(Hj[j]-Hnj));
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
      PHIST_TG_PREFIX(LARTGP)(Hj[j-1],Hj[j],&S[i]->cs_[j-1],&S[i]->sn_[j-1],&tmp);
      //{
        //_MT_ len = mt::sqrt(st::real(st::conj(Hj[j-1])*Hj[j-1])+st::real(st::conj(Hj[j])*Hj[j]));
        //S[i]->cs_[j-1] = Hj[j-1]/len;
        //S[i]->sn_[j-1] = Hj[j]/len;
        //// and apply it
        //tmp = st::conj(S[i]->cs_[j-1])*Hj[j-1] + st::conj(S[i]->sn_[j-1])*Hj[j];
      //}
#endif
#ifdef PHIST_TESTING
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
      PHIST_CHK_IERR(SUBR(sdMat_to_device)(S[i]->H_,iflag),*iflag);
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
      else if(S[i]->totalIter >= maxIter )
      {
        S[i]->status=3;
        anyFailed++;
      }
      else if( S[i]->curDimV_ >= mvecBuff->size() )
      {
        S[i]->status = 2; // mark asrestart needed
        anyFull++;
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
PHIST_TASK_END(iflag)

#if PHIST_OUTLEV>=PHIST_VERBOSE
    PHIST_SOUT(PHIST_VERBOSE,"%d converged, %d need restart", anyConverged,anyFull);
    if (anyFailed)
    {
      PHIST_SOUT(PHIST_VERBOSE,", %d exceeded max iter.\n",anyFailed);
    }
    else
    {
      PHIST_SOUT(PHIST_VERBOSE,"\n");
    }
    PHIST_SOUT(PHIST_VERBOSE,"---------------------------------------  \n");
#endif

  // delete views
  PHIST_CHK_IERR(SUBR(mvec_delete)(work_x, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_delete)(work_y, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_delete)(Vj,     iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_delete)(Vk,     iflag), *iflag);

  *iflag=99;
  for (int i=0; i<numSys; i++) *iflag=std::min(*iflag,S[i]->status);
}

