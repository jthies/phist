// create new state objects. We just get an array of (NULL-)pointers
void SUBR(jadaInnerGmresStates_create)(TYPE(jadaInnerGmresState_ptr) state[], int numSys,
        const_map_ptr_t map, int maxBas,int* ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  if (numSys==0) return;
  const_comm_ptr_t comm;
  PHIST_CHK_IERR(phist_map_get_comm(map,&comm,ierr),*ierr);

  // setup a "queue" of mvecs to use later
  int totally_needed_mvecs = 2*(maxBas+2)*numSys;
  std::vector<TYPE(mvec_ptr)> *unused_mvecs = new std::vector<TYPE(mvec_ptr)>();
  for(int i = 0; i < totally_needed_mvecs; i++)
  {
    TYPE(mvec_ptr) tmp;
    PHIST_CHK_IERR(SUBR(mvec_create)(&tmp,map,numSys,ierr),*ierr);
    unused_mvecs->push_back(tmp);
  }
  std::vector<TYPE(mvec_ptr)> *used_mvecs = new std::vector<TYPE(mvec_ptr)>();
  
  for (int i=0;i<numSys;i++)
  {
    state[i] = new TYPE(jadaInnerGmresState);
    state[i]->id=i;
    // set some default options
    state[i]->tol=0.5; // typical starting tol for JaDa inner iterations...
    state[i]->ierr=-2;// not initialized
    state[i]->maxBas_=maxBas;
    // we allow one additional vector to be stored in the basis so that
    // we can have V(:,i+1) = A*V(:,i) temporarily
    state[i]->V_ = new TYPE(mvec_ptr)[maxBas+1];
    state[i]->AV_ = new TYPE(mvec_ptr)[maxBas+1];
    state[i]->glob_unused_mvecs_ = (void*)unused_mvecs;
    state[i]->glob_used_mvecs_ = (void*)used_mvecs;
    PHIST_CHK_IERR(SUBR(sdMat_create)(&state[i]->H_, maxBas+1, maxBas, comm,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_create)(&state[i]->x0_,map,1,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_create)(&state[i]->b_,map,1,ierr),*ierr);
    state[i]->cs_ = new MT[maxBas];
    state[i]->sn_ = new ST[maxBas];
    state[i]->rs_ = new ST[maxBas];
    state[i]->curDimV_=0;
    state[i]->normR0_=-mt::one(); // not initialized
    state[i]->normR_=-mt::one();
  }
}


//! delete jadaInnerGmresState object
void SUBR(jadaInnerGmresStates_delete)(TYPE(jadaInnerGmresState_ptr) state[], int numSys, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  if( numSys > 0 )
  {
    std::vector<TYPE(mvec_ptr)> *unused_mvecs = (std::vector<TYPE(mvec_ptr)>*) state[0]->glob_unused_mvecs_;
    std::vector<TYPE(mvec_ptr)> *used_mvecs = (std::vector<TYPE(mvec_ptr)>*) state[0]->glob_used_mvecs_;
    for(int i = 0; i < unused_mvecs->size(); i++)
    {
      PHIST_CHK_IERR(SUBR(mvec_delete)(unused_mvecs->at(i),ierr), *ierr);
    }
    for(int i = 0; i < used_mvecs->size(); i++)
    {
      PHIST_CHK_IERR(SUBR(mvec_delete)(used_mvecs->at(i),ierr), *ierr);
    }
    delete unused_mvecs;
    delete used_mvecs;
  }
  for (int i=0;i<numSys;i++)
  {
    PHIST_CHK_IERR(SUBR(mvec_delete)(state[i]->x0_,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_delete)(state[i]->b_,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_delete)(state[i]->H_,ierr),*ierr);
    delete [] state[i]->AV_;
    delete [] state[i]->V_;
    delete [] state[i]->cs_;
    delete [] state[i]->sn_;
    delete [] state[i]->rs_;
    delete state[i];
  }
}

// reset jadaInnerGmres state.
void SUBR(jadaInnerGmresState_reset)(TYPE(jadaInnerGmresState_ptr) S, TYPE(const_mvec_ptr) b,
        TYPE(const_mvec_ptr) x0,int *ierr)
{
#include "phist_std_typedefs.hpp"  
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  
  if (b==NULL && S->normR0_ == -mt::one())
  {
    PHIST_OUT(PHIST_ERROR,"on the first call to jadaInnerGmresState_reset you *must* provide the RHS vector");
    *ierr=-1;
    return;
  }
  else if (b!=NULL)
  {
    // new rhs -> need to recompute ||b-A*x0||
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(), b, st::zero(), S->b_, ierr), *ierr);
    S->ierr = -1;
    S->totalIter = 0;
  }
  S->curDimV_=0;
  S->normR0_=-mt::one(); // needs to be computed in next iterate call
  // set V_0=X_0. iterate() will have to compute the residual and normalize it,
  // because the actual V_0 we want is r/||r||_2, but we can't easily apply the
  // operator to a single vector.
  PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),x0,st::zero(), S->x0_, ierr),*ierr);
  for (int i=0;i<S->maxBas_;i++)
    S->rs_[i]=st::zero();
  PHIST_CHK_IERR(SUBR(sdMat_put_value)(S->H_,st::zero(),ierr),*ierr);
}


void SUBR(jadaInnerGmresStates_updateSol)(TYPE(jadaInnerGmresState_ptr) S_array[], int numSys, TYPE(mvec_ptr) x, TYPE(mvec_ptr) Ax, _MT_* resNorm, bool scaleSolutionToOne, int* ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);

  const_comm_ptr_t comm=NULL;
  PHIST_CHK_IERR(SUBR(mvec_get_comm)(x,&comm,ierr),*ierr);
  TYPE(mvec_ptr) x_i=NULL;
  TYPE(mvec_ptr) Ax_i=NULL;
  TYPE(mvec_ptr) V=NULL;
  *ierr = 0;

  if( numSys <= 0 )
    return;

  int sharedCurDimV = S_array[0]->curDimV_;
  for(int i = 0; i < numSys; i++)
  {
    PHIST_CHK_IERR(*ierr = (sharedCurDimV == S_array[i]->curDimV_) ? 0 : -99, *ierr);
  }

  // no iteration done yet?
  if( sharedCurDimV <= 0 )
  {
    for(int i = 0; i < numSys; i++)
      resNorm[i] = -mt::one();
    return;
  }
  if( sharedCurDimV == 1 )
  {
    for(int i = 0; i < numSys; i++)
      resNorm[i] = S_array[i]->normR0_;
    return;
  }

  _ST_ *y = new _ST_[numSys*S_array[0]->maxBas_];
  int ldy = numSys;

  // calculate y by solving the triangular systems
  for (int i=0;i<numSys;i++)
  {
    TYPE(const_jadaInnerGmresState_ptr) S = S_array[i];
    ST *H_raw=NULL;
    lidx_t ldH;

    int m=S->curDimV_-1;
    // no iteration done yet?
    if( m < 0 )
    {
      resNorm[i] = -mt::one();
      continue;
    }
    if( m == 0 )
    {
      resNorm[i] = S->normR0_;
      continue;
    }
    resNorm[i] = S->normR_/S->normR0_;

    PHIST_CHK_IERR(SUBR(sdMat_extract_view)(S->H_,&H_raw,&ldH,ierr),*ierr);

#if PHIST_OUTLEV>=PHIST_DEBUG
    PHIST_SOUT(PHIST_DEBUG,"jadaInnerGmres_updateSol[%d], curDimV=%d, H=\n",i,S->curDimV_);
    {
      TYPE(sdMat_ptr) H = NULL;
      PHIST_CHK_IERR(SUBR(sdMat_view_block)(S->H_, &H, 0, m+1, 0, m, ierr), *ierr);
      PHIST_CHK_IERR(SUBR(sdMat_print)(H,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(sdMat_delete)(H,ierr),*ierr);
    }
    PHIST_SOUT(PHIST_DEBUG,"rs=\n");
    for (int k=0;k<m;k++)
    {
      PHIST_SOUT(PHIST_DEBUG,"%16.8f+%16.8fi\n",st::real(S->rs_[k]),st::imag(S->rs_[k]));      
    }
#endif

    // y = H\rs, H upper triangular
    const char* uplo="U";
    const char* trans="N";
    const char* diag="N";
    // set y to rs
    for(int j = 0; j < m; j++)
      y[i+ldy*j] = S->rs_[j];
    PHIST_CHK_IERR(PREFIX(TRSV)(uplo,trans,diag,&m,
                                        (st::blas_scalar_t*)H_raw,&ldH,
                                        (st::blas_scalar_t*)&y[i], &ldy, ierr),*ierr);

    // if we are only interested in the directions Vi*yi and appropriate AVi*yi,
    // then this scaling may help to improve the conditioning of a following orthogonlization step!
    if( scaleSolutionToOne )
    {
      // scale y to one
      _MT_ scale = mt::zero();
      for(int j = 0; j < m; j++)
        scale += st::real(st::conj(y[i+ldy*j])*y[i+ldy*j]);
      scale = mt::one()/sqrt(scale);
      for(int j = 0; j < m; j++)
        y[i+ldy*j] *= scale;
    }
  }


  // add up solution
  for(int j = 0; j < sharedCurDimV-1; j++)
  {
    PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(&y[ldy*j], S_array[0]->V_[j], st::one(), x, ierr), *ierr);
    if( Ax != NULL )
    {
      PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(&y[ldy*j], S_array[0]->AV_[j], st::one(), Ax, ierr), *ierr);
    }
  }
}


// implementation of gmres on several systems simultaneously
void SUBR(jadaInnerGmresStates_iterate)(TYPE(const_op_ptr) jdOp,
        TYPE(jadaInnerGmresState_ptr) S[], int numSys,
        int* nIter, int* ierr)
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

  // get map
  const_map_ptr_t map = NULL;
  PHIST_CHK_IERR( SUBR(mvec_get_map) (S[0]->V_, &map, ierr), *ierr);
  
  // work vector for x and y = jdOp(x)
  TYPE(mvec_ptr) work_x = NULL;
  TYPE(mvec_ptr) work_y = NULL;
  // get pointers to mvec-stacks
  std::vector<TYPE(mvec_ptr)> *unused_mvecs = (std::vector<TYPE(mvec_ptr)>*) S[0]->glob_unused_mvecs_;
  std::vector<TYPE(mvec_ptr)> *used_mvecs = (std::vector<TYPE(mvec_ptr)>*) S[0]->glob_used_mvecs_;
#warning "assuming full restart, partial restart not supported yet"
  unused_mvecs->insert(unused_mvecs->end(),used_mvecs->begin(),used_mvecs->end());
  used_mvecs->clear();
  // borrow work_x, it is only used for one iteration
  PHIST_CHK_IERR( *ierr = (unused_mvecs->size() > 2) ? 0 : -1, *ierr);
  work_x = unused_mvecs->at(0);
  {
    // check dimensions
    int nvec;
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(work_x,&nvec,ierr),*ierr);
    PHIST_CHK_IERR( *ierr = (nvec == numSys) ? 0 : -99, *ierr);
  }

  // views into V_
  TYPE(mvec_ptr) Vj = NULL, AVj = NULL, Vprev = NULL;
  // views into H_
  TYPE(sdMat_ptr) R1 = NULL, R2 = NULL;
  TYPE(mvec_ptr) view_Ax = NULL;

  // we return as soon as one system converges or reaches its
  // maximum permitted number of iterations. The decision about what to do
  // next is then left to the caller.
  int anyConverged = 0;
  int anyFailed = 0;
  // check wether one of the systems cannot iterate
  for(int i = 0; i < numSys; i++)
  {
    if( S[i]->curDimV_ >= S[i]->maxBas_ )
      anyFailed++;
  }


  // gather work_x
#ifndef PHIST_KERNEL_LIB_FORTRAN
  for(int i = 0; i < numSys; i++)
  {
    int jprev = S[i]->curDimV_-1;
    if( jprev < 0 )
    {
      // (re-)start
      PHIST_CHK_IERR( SUBR(mvec_set_block ) (work_x, S[i]->x0_, i, i, ierr), *ierr);
    }
    else
    {
      PHIST_CHK_IERR( SUBR(mvec_view_block) (S[i]->V_, &Vj, jprev, jprev, ierr), *ierr);
      PHIST_CHK_IERR( SUBR(mvec_set_block ) (work_x, Vj, i, i, ierr), *ierr);
    }
  }
#else
  {
    TYPE(mvec_ptr) work_xi[numSys];
    for(int i = 0; i < numSys; i++)
    {
      work_xi[i] = NULL;
      int jprev = S[i]->curDimV_-1;
      if( jprev < 0 )
      {
        // (re-)start
        PHIST_CHK_IERR( SUBR(mvec_view_block ) (S[i]->x0_, &work_xi[i], 0, 0, ierr), *ierr);
      }
      else
      {
        PHIST_CHK_IERR( SUBR(mvec_view_block) (S[i]->V_[jprev], &work_xi[i], i, i, ierr), *ierr);
      }
    }
    PHIST_CHK_IERR( SUBR(mvec_gather_mvecs) (work_x, (TYPE(const_mvec_ptr)*)work_xi, numSys, ierr), *ierr);
    for(int i = 0; i < numSys; i++)
    {
      PHIST_CHK_IERR( SUBR(mvec_delete) (work_xi[i], ierr), *ierr);
    }
  }
#endif

  while( anyConverged == 0 && anyFailed == 0 )
  {
    // we need a new mvec for work_y
    work_y = unused_mvecs->back();
    used_mvecs->push_back(work_y);
    unused_mvecs->pop_back();

    //    % apply the jadaOp
    PHIST_CHK_IERR( jdOp->apply (st::one(), jdOp->A, work_x, st::zero(), work_y, ierr), *ierr);
    PHIST_CHK_IERR( SUBR(jadaOp_view_AX) (jdOp->A, &view_Ax, ierr), *ierr);

    //    % scatter view_Ax
    // we need storage for Ax, reuse work_x here
    work_x = unused_mvecs->back();
    used_mvecs->push_back(work_x);
    unused_mvecs->pop_back();
    // copy data
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(), view_Ax, st::zero(), work_x, ierr), *ierr);
    // append it to the AV lists of the states
    for(int i = 0; i < numSys; i++)
    {
      int jprev = std::max(S[i]->curDimV_-1,0);
      S[i]->AV_[jprev] = work_x;
    }

    // use work_y for Vj to avoid unnecessary copies and scatter afterwords

    //    % initialize GMRES for (re-)started systems
    for(int i = 0; i < numSys; i++)
    {
      int j = S[i]->curDimV_;
      PHIST_CHK_IERR( SUBR(mvec_view_block) (work_y, &Vj, i, i, ierr), *ierr);

      if( j == 0 )
      {
        //    % (re-)start: normalize r_0 = b - A*x_0
        PHIST_SOUT(PHIST_VERBOSE,"jadaInnerGmres state %d (re-)starts\n",i);

        // we need b-Ax0 in the first step
        PHIST_CHK_IERR( SUBR(mvec_add_mvec) (st::one(), S[i]->b_, -st::one(), Vj, ierr), *ierr);

        // normalize
        _MT_ norm;
        PHIST_CHK_IERR( SUBR(mvec_normalize) (Vj, &norm, ierr), *ierr);
        S[i]->rs_[0] = norm;
        S[i]->normR_ = norm;

        // if this is no restart, store normR0
        if( S[i]->normR0_ < mt::zero() )
          S[i]->normR0_ = norm;
      }
    }

    // TODO: % arnoldi update with modified gram schmidt
    {
      int sharedCurDimV = S[0]->curDimV_;
      for(int i = 0; i < numSys; i++)
      {
        PHIST_CHK_IERR(*ierr = (sharedCurDimV == S[i]->curDimV_) ? 0 : -99, *ierr);
      }
      int j = sharedCurDimV;
      if( j > 0 )
      {
        _ST_*tmp = new _ST_[numSys];
        for(int k = 0; k < sharedCurDimV; k++)
        {
          // simply gram schmidt applied to all systems at once
          PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(work_y,S[0]->V_[k],tmp,ierr),*ierr);
          // store in H
          for(int i = 0; i < numSys; i++)
          {
            // raw view of H
            ST *Hj=NULL;
            lidx_t ldH; 
            PHIST_CHK_IERR(SUBR(sdMat_extract_view)(S[i]->H_,&Hj,&ldH,ierr),*ierr); 
            Hj += (j-1)*ldH;
            Hj[k] = tmp[i];
          }
          for(int i = 0; i < numSys; i++)
            tmp[i] = -tmp[i];
          PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(tmp,S[0]->V_[k],st::one(),work_y,ierr),*ierr);
        }
        delete[] tmp;

        _MT_ *tmp_ = new _MT_[numSys];
        PHIST_CHK_IERR(SUBR(mvec_normalize)(work_y,tmp_,ierr),*ierr);
        for(int i = 0; i < numSys; i++)
        {
          // raw view of H
          ST *Hj=NULL;
          lidx_t ldH; 
          PHIST_CHK_IERR(SUBR(sdMat_extract_view)(S[i]->H_,&Hj,&ldH,ierr),*ierr); 
          Hj += (j-1)*ldH;
          Hj[j] = tmp_[i];
        }
        delete[] tmp_;
      }
    }


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
        tmp = S[i]->cs_[k]*Hj[k] + S[i]->sn_[k]*Hj[k+1];
        Hj[k+1] = -st::conj(S[i]->sn_[k])*Hj[k] + S[i]->cs_[k]*Hj[k+1];
      }
      // new Givens rotation to eliminate H(j+1,j)
#ifdef IS_COMPLEX
      PREFIX(LARTG)((blas_cmplx_t*)&Hj[j-1],(blas_cmplx_t*)&Hj[j],&S[i]->cs_[j-1],(blas_cmplx_t*)&S[i]->sn_[j-1],(blas_cmplx_t*)&tmp);
#else
      PREFIX(LARTG)(&Hj[j-1],&Hj[j],&S[i]->cs_[j-1],&S[i]->sn_[j-1],&tmp);
#endif
#ifdef TESTING
{
  PHIST_OUT(PHIST_VERBOSE,"(Hj[j-1],Hj[j]) = (%8.4e+%8.4ei, %8.4e+%8.4ei)\n", st::real(Hj[j-1]), st::imag(Hj[j-1]),st::real(Hj[j]),st::imag(Hj[j]));
  PHIST_OUT(PHIST_VERBOSE,"(c,s) = (%8.4e, %8.4e+%8.4ei)\n", S[i]->cs_[j-1],st::real(S[i]->sn_[j-1]),st::imag(S[i]->sn_[j-1]));
  PHIST_OUT(PHIST_VERBOSE,"r = %8.4e+%8.4ei\n", st::real(tmp),st::imag(tmp));
  _ST_ r_ = S[i]->cs_[j-1]*Hj[j-1] + S[i]->sn_[j-1]*Hj[j];
  _ST_ zero_ = -st::conj(S[i]->sn_[j-1])*Hj[j-1] + S[i]->cs_[j-1]*Hj[j];
  PHIST_OUT(PHIST_VERBOSE,"(r, 0) = (%8.4e+8.4%ei, %8.4e+8.4%ei)\n", st::real(r_), st::imag(r_), st::real(zero_), st::imag(zero_));
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
      if( relres < S[i]->tol || absres < S[i]->tol )
      {
        S[i]->ierr = 0; // mark as converged
        anyConverged++;
      }
      else if( S[i]->curDimV_ >= S[i]->maxBas_ )
      {
        S[i]->ierr = 2; // mark as failed/restart needed
        anyFailed++;
      }
      else
      {
        S[i]->ierr = 1; // iterating, not converged yet
      }
    }

    // scatter work_y
    // append work_y to the V_ lists of the states
    for(int i = 0; i < numSys; i++)
    {
      int j = S[i]->curDimV_-1;
      S[i]->V_[j] = work_y;
    }
    // use work_y as input in the next iteration
    work_x = work_y;



    PHIST_SOUT(PHIST_VERBOSE,"GMRES iteration states\n");
    PHIST_SOUT(PHIST_VERBOSE,"======================\n");
#if PHIST_OUTLEV>=PHIST_VERBOSE
    for (int i=0;i<numSys;i++)
    {
    PHIST_SOUT(PHIST_VERBOSE,"[%d]: %d\t%8.4e\t(%8.4e)\n",i,
          S[i]->curDimV_-1,S[i]->normR_/S[i]->normR0_,S[i]->normR_);
    }
#endif
    PHIST_SOUT(PHIST_VERBOSE,"%d converged, %d failed.\n",anyConverged,anyFailed);
    PHIST_SOUT(PHIST_VERBOSE,"----------------------\n");

    (*nIter)++;
  }

  if (anyConverged > 0)
    *ierr=0;
      
  if (anyFailed > 0)
    *ierr=1;
}

