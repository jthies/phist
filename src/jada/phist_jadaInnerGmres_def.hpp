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

  // memory for the bases V is allocated in one big chunk
  int tot_nV = (maxBas+2)*numSys;
//TYPE(mvec_ptr) Vglob=NULL;
//PHIST_CHK_IERR(SUBR(mvec_create)(&Vglob, map,tot_nV, ierr),*ierr);
  
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
    PHIST_CHK_IERR(SUBR(mvec_create)(&state[i]->V_,map,maxBas+1,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_create)(&state[i]->AV_,map,maxBas+1,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_create)(&state[i]->H_, maxBas+1, maxBas, comm,ierr),*ierr);
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
  for (int i=0;i<numSys;i++)
  {
    PHIST_CHK_IERR(SUBR(mvec_delete)(state[i]->b_,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_delete)(state[i]->H_,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_delete)(state[i]->AV_,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_delete)(state[i]->V_,ierr),*ierr);
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
  PHIST_CHK_IERR(SUBR(mvec_set_block)(S->V_,x0,0,0,ierr),*ierr);
  for (int i=0;i<S->maxBas_;i++)
    S->rs_[i]=st::zero();
  PHIST_CHK_IERR(SUBR(sdMat_put_value)(S->H_,st::zero(),ierr),*ierr);
}


void SUBR(jadaInnerGmresStates_updateSol)(TYPE(jadaInnerGmresState_ptr) S_array[], int numSys, TYPE(mvec_ptr) x, TYPE(mvec_ptr) Ax, _MT_* resNorm, int* ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);

  const_comm_ptr_t comm=NULL;
  PHIST_CHK_IERR(SUBR(mvec_get_comm)(x,&comm,ierr),*ierr);
  TYPE(mvec_ptr) x_i=NULL;
  TYPE(mvec_ptr) Ax_i=NULL;

  for (int i=0;i<numSys;i++)
  {
    PHIST_CHK_IERR(SUBR(mvec_view_block)(x,&x_i,i,i,ierr),*ierr);
    if( Ax != NULL )
    {
      PHIST_CHK_IERR(SUBR(mvec_view_block)(Ax,&Ax_i,i,i,ierr),*ierr);
    }
    // compute y by solving the triangular system
    TYPE(const_jadaInnerGmresState_ptr) S = S_array[i];
    TYPE(sdMat_ptr) y=NULL;
    ST *H_raw=NULL, *y_raw=NULL;
    lidx_t ldH,ldy;

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

    PHIST_CHK_IERR(SUBR(sdMat_create)(&y,m,1,comm,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_extract_view)(y,&y_raw,&ldy,ierr),*ierr);
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
    for (int i=0;i<m;i++)
    {
      PHIST_SOUT(PHIST_DEBUG,"%16.8f+%16.8fi\n",st::real(S->rs_[i]),st::imag(S->rs_[i]));      
    }
#endif

    for (int i=0;i<m;i++)
      y_raw[i]=S->rs_[i];

    // y = H\rs, H upper triangular
    const char* uplo="U";
    const char* trans="N";
    const char* diag="N";
    int nrhs=1;
    // set y to rs
    for(int j = 0; j < m; j++)
      y_raw[j] = S->rs_[j];
    PHIST_CHK_IERR(PREFIX(TRTRS)(uplo,trans,diag,&m,&nrhs,
                                        (st::blas_scalar_t*)H_raw,&ldH,
                                        (st::blas_scalar_t*)y_raw, &ldy, ierr),*ierr);


#if PHIST_OUTLEV>=PHIST_DEBUG
    PHIST_OUT(PHIST_DEBUG,"y=\n");
    PHIST_CHK_IERR(SUBR(sdMat_print)(y,ierr),*ierr);
#endif

    TYPE(mvec_ptr) V=NULL;
    PHIST_CHK_IERR(SUBR(mvec_view_block)(S->V_,&V,0,m-1,ierr),*ierr);

    // X = X + M\V*y. TODO: with right preconditioning, split this into two parts and apply
    // preconditioner to all systems simultaneously outside the loop (need tmp vector)
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(st::one(),V,y,st::one(),x_i,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_delete)(V,ierr),*ierr);
    if( Ax != NULL )
    {
      TYPE(mvec_ptr) AV=NULL;
      PHIST_CHK_IERR(SUBR(mvec_view_block)(S->AV_,&AV,0,m-1,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(st::one(),AV,y,st::one(),Ax_i,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_delete)(AV,ierr),*ierr);
    }
    PHIST_CHK_IERR(SUBR(sdMat_delete)(y,ierr),*ierr);
  }
  PHIST_CHK_IERR(SUBR(mvec_delete)(x_i,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_delete)(Ax_i,ierr),*ierr);
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
  const_map_ptr_t map;
  PHIST_CHK_IERR( SUBR(mvec_get_map) (S[0]->V_, &map, ierr), *ierr);
  
  // work vector for x and y = jdOp(x)
  TYPE(mvec_ptr) work_x;
  PHIST_CHK_IERR( SUBR(mvec_create ) (&work_x, map, numSys, ierr), *ierr);
  TYPE(mvec_ptr) work_y;
  PHIST_CHK_IERR( SUBR(mvec_create ) (&work_y, map, numSys, ierr), *ierr);

  // views into V_
  TYPE(mvec_ptr) Vj = NULL, AVj = NULL, Vprev = NULL;
  // views into H_
  TYPE(sdMat_ptr) R1 = NULL, R2 = NULL;
  TYPE(mvec_ptr) work_Ax = NULL;

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

  while( anyConverged == 0 && anyFailed == 0 )
  {
    // gather work_x
    for(int i = 0; i < numSys; i++)
    {
      int jprev = std::max(S[i]->curDimV_-1,0);
      PHIST_CHK_IERR( SUBR(mvec_view_block) (S[i]->V_, &Vj, jprev, jprev, ierr), *ierr);
      PHIST_CHK_IERR( SUBR(mvec_set_block ) (work_x, Vj, i, i, ierr), *ierr);
    }

    // apply the jadaOp
    PHIST_CHK_IERR( jdOp->apply (st::one(), jdOp->A, work_x, st::zero(), work_y, ierr), *ierr);
    PHIST_CHK_IERR( SUBR(jadaOp_view_AX) (jdOp->A, &work_Ax, ierr), *ierr);

    // distribute work_y
    for(int i = 0; i < numSys; i++)
    {
      int j = S[i]->curDimV_;
      PHIST_CHK_IERR( SUBR(mvec_view_block) (S[i]->V_, &Vj, j, j, ierr), *ierr);
      PHIST_CHK_IERR( SUBR(mvec_get_block ) (work_y, Vj, i, i, ierr), *ierr);

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
      else
      {
        // set block in AV
        PHIST_CHK_IERR( SUBR(mvec_view_block) (S[i]->AV_, &AVj, j-1, j-1, ierr), *ierr);
        PHIST_CHK_IERR( SUBR(mvec_get_block ) (work_Ax, AVj, i, i, ierr), *ierr);

        //    % arnoldi update:
        // view appropriate blocks in V_ and H_
        PHIST_CHK_IERR( SUBR(mvec_view_block ) (S[i]->V_, &Vprev, 0, j-1, ierr), *ierr);
        // view H(j,j) as R1
        PHIST_CHK_IERR(SUBR(sdMat_view_block)(S[i]->H_,&R1,j,j,j-1,j-1,ierr),*ierr);
        // view H(1:j,j) as R2
        PHIST_CHK_IERR(SUBR(sdMat_view_block)(S[i]->H_,&R2,0,j-1,j-1,j-1,ierr),*ierr);
        //orthogonalize
        PHIST_CHK_IERR(SUBR(orthog)(Vprev,Vj,R1,R2,2,ierr),*ierr);


        //    % update QR factorization of H
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
  PHIST_OUT(PHIST_VERBOSE,"(Hj[j-1],Hj[j]) = (%8.4e, %8.4e)\n", Hj[j-1],Hj[j]);
  PHIST_OUT(PHIST_VERBOSE,"(c,s) = (%8.4e, %8.4e)\n", S[i]->cs_[j-1],S[i]->sn_[j-1]);
  PHIST_OUT(PHIST_VERBOSE,"r = %8.4e\n", tmp);
  _ST_ r_ = S[i]->cs_[j-1]*Hj[j-1] + S[i]->sn_[j-1]*Hj[j];
  _ST_ zero_ = -st::conj(S[i]->sn_[j-1])*Hj[j-1] + S[i]->cs_[j-1]*Hj[j];
  PHIST_OUT(PHIST_VERBOSE,"(r, 0) = (%8.4e, %8.4e)\n", r_, zero_);
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

    PHIST_SOUT(PHIST_VERBOSE,"GMRES iteration states\n");
    PHIST_SOUT(PHIST_VERBOSE,"======================\n");
#if PHIST_OUTLEV>=PHIST_VERBOSE
    for (int i=0;i<numSys;i++)
    {
    PHIST_SOUT(PHIST_VERBOSE,"[%d]: %d\t%8.4e\n",i,
          S[i]->curDimV_-1,S[i]->normR_);
    }
#endif
    PHIST_SOUT(PHIST_VERBOSE,"%d converged, %d failed.\n",anyConverged,anyFailed);
    PHIST_SOUT(PHIST_VERBOSE,"----------------------\n");

    (*nIter)++;
  }

  // delete work storage
  PHIST_CHK_IERR( SUBR(mvec_delete) (work_x, ierr), *ierr);
  PHIST_CHK_IERR( SUBR(mvec_delete) (work_y, ierr), *ierr);
  PHIST_CHK_IERR( SUBR(mvec_delete) (work_Ax, ierr), *ierr);

  if (anyConverged > 0)
    *ierr=0;
      
  if (anyFailed > 0)
    *ierr=1;
}

