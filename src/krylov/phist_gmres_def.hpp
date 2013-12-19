//! small helper function for GMRES to compute
// a Givens rotation such that
// | c s ||f| |r|
// | _   || |=| |
// |-s c ||g| |0|
// This function ses the convention for the BLAS,
// cf. for instance Bindel et al. 2001). The
// BLAS routine is called XROTG, however, it does
// not do what I want in the complex case (TODO is there
// another routine that gives complex coefficients??)
void SUBR(rotg)(_ST_ f, _ST_ g, _ST_& cs, _ST_& sn, _ST_& r)
{
#include "phist_std_typedefs.hpp"
  MT af=st::abs(f);
  MT ag=st::abs(g);
  MT d=mt::sqrt(af*af + ag*ag);
  cs=st::zero();
  sn=st::zero();
  r=st::zero();

  if (ag==mt::zero() ) // includes case f=0 and g=0
  {
    cs = st::one();
    sn = st::zero();
    r  = f;
  }
  else if (af == mt::zero())
  {
    cs=st::zero();
    sn=st::conj(g)/ag;
    r=ag;
  }
  else
  {
    cs=af/d;
    sn=(f/af)*(st::conj(g)/d);
    r=(f/af)*d;
  }
  return;
}

// create new state objects. We just get an array of (NULL-)pointers
void SUBR(gmresStates_create)(TYPE(gmresState_ptr) state[], int numSys,
        const_map_ptr_t map, int maxBas,int* ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  if (numSys==0) return;
  const_comm_ptr_t comm;
  PHIST_CHK_IERR(phist_map_get_comm(map,&comm,ierr),*ierr);

  // memory for the bases V is allocated in one big chunk
  int tot_nV = (maxBas+1)*numSys;
  TYPE(mvec_ptr) Vglob=NULL;
  PHIST_CHK_IERR(SUBR(mvec_create)(&Vglob, map,tot_nV, ierr),*ierr);
  
  for (int i=0;i<numSys;i++)
  {
    state[i] = new TYPE(gmresState);
    state[i]->id=i;
    // set some default options
    state[i]->tol=0.5; // typical starting tol for JaDa inner iterations...
    state[i]->ierr=-2;// not initialized
    state[i]->Vglob_=Vglob;
    state[i]->offsetVglob_=i*(maxBas+1);
    state[i]->maxBas_=maxBas;
    // we allow one additional vector to be stored in the basis so that
    // we can have V(:,i+1) = A*V(:,i) temporarily
    state[i]->V_=NULL;
    PHIST_CHK_IERR(SUBR(mvec_view_block)(Vglob,&state[i]->V_,
        state[i]->offsetVglob_, state[i]->offsetVglob_+maxBas, ierr),*ierr);
  
    PHIST_CHK_IERR(SUBR(sdMat_create)(&state[i]->H_, maxBas+1, maxBas, comm,ierr),*ierr);
    state[i]->cs_ = new ST[maxBas];
    state[i]->sn_ = new ST[maxBas];
    state[i]->rs_ = new ST[maxBas];
    state[i]->curDimV_=0;
    state[i]->normB_=-mt::one(); // not initialized
    // not initialized, set B pointer to NULL
    state[i]->B_=NULL;
  }
}


//! delete gmresState object
void SUBR(gmresStates_delete)(TYPE(gmresState_ptr) state[], int numSys, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  if (numSys==0) return;
  TYPE(mvec_ptr) Vglob=state[0]->Vglob_;
  for (int i=0;i<numSys;i++)
  {
    PHIST_CHK_IERR(SUBR(mvec_delete)(state[i]->V_,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_delete)(state[i]->H_,ierr),*ierr);
    delete [] state[i]->cs_;
    delete [] state[i]->sn_;
    delete [] state[i]->rs_;
    delete state[i];
  }
  PHIST_CHK_IERR(SUBR(mvec_delete)(Vglob,ierr),*ierr);
}

// reset gmres state.
void SUBR(gmresState_reset)(TYPE(gmresState_ptr) S, TYPE(const_mvec_ptr) b,
        TYPE(const_mvec_ptr) x0,int *ierr)
{
#include "phist_std_typedefs.hpp"  
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  
  if (b==NULL && S->B_==NULL)
  {
    PHIST_OUT(PHIST_ERROR,"on the first call to gmresState_reset you *must* provide the RHS vector");
    *ierr=-1;
    return;
  }
  else if (b!=S->B_ && b!=NULL)
  {
    // new rhs -> need to recompute ||B||
    S->B_=b;
    S->normB_=-mt::one(); // needs to be computed in next iterate call
  }
  S->curDimV_=0;
  S->normR0_=-mt::one(); // needs to be computed in nect iterate call
  // set V_0=X_0. iterate() will have to compute the residual and normalize it,
  // because the actual V_0 we want is r/||r||_2, but we can't easily apply the
  // operator to a single vector.
  PHIST_CHK_IERR(SUBR(mvec_set_block)(S->V_,x0,0,0,ierr),*ierr);
  for (int i=0;i<S->maxBas_;i++)
  {
    S->rs_[i]=st::zero();
  }
  PHIST_CHK_IERR(SUBR(sdMat_put_value)(S->H_,st::zero(),ierr),*ierr);
  return;
}

void SUBR(gmresStates_updateSol)(TYPE(gmresState_ptr) S_array[], int numSys, TYPE(mvec_ptr) x, int* ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);

  const_comm_ptr_t comm=NULL;
  PHIST_CHK_IERR(SUBR(mvec_get_comm)(x,&comm,ierr),*ierr);
  TYPE(mvec_ptr) x_i=NULL;

  // TODO - omp parallel for
  for (int i=0;i<numSys;i++)
  {
    PHIST_CHK_IERR(SUBR(mvec_view_block)(x,&x_i,i,i,ierr),*ierr);
    // compute y by solving the triangular system
    TYPE(const_gmresState_ptr) S = S_array[i];
    TYPE(sdMat_ptr) y=NULL;
    ST *H_raw=NULL, *y_raw=NULL;
    lidx_t ldH,ldy;

    int m=S->curDimV_-1;
    
    PHIST_CHK_IERR(SUBR(sdMat_create)(&y,m,1,comm,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_extract_view)(y,&y_raw,&ldy,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_extract_view)(S->H_,&H_raw,&ldH,ierr),*ierr);

#if PHIST_OUTLEV>=PHIST_DEBUG
    PHIST_SOUT(PHIST_DEBUG,"gmres_updateSol[%d], curDimV=%d, H=\n",i,S->curDimV_);
    PHIST_CHK_IERR(SUBR(sdMat_print)(S->H_,ierr),*ierr);
    PHIST_SOUT(PHIST_DEBUG,"rs=\n");
    for (int i=0;i<m;i++)
    {
      PHIST_SOUT(PHIST_DEBUG,"%16.8f+%16.8fi\n",st::real(S->rs_[i]),st::imag(S->rs_[i]));      
    }
#endif

    for (int i=0;i<m;i++)
    {
      y_raw[i]=S->rs_[i];
    }

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
    PHIST_CHK_IERR(SUBR(sdMat_delete)(y,ierr),*ierr);
  }
  PHIST_CHK_IERR(SUBR(mvec_delete)(x_i,ierr),*ierr);
  return;
}

// implementation of gmres on several systems simultaneously
void SUBR(gmresStates_iterate)(TYPE(const_op_ptr) Op,
        TYPE(gmresState_ptr) S_array[], int numSys,
        int* ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
  *ierr = 0;

  // multi-vectors to work with several of the systems simultaneously
  TYPE(mvec_ptr) V=NULL, W=NULL;
  // views into vector blocks of system i
  TYPE(mvec_ptr) Vprev=NULL, Vj=NULL;
  // for the orthog routine (again, a view into the H objects in the states
  TYPE(sdMat_ptr) R1=NULL,R2=NULL;

#if PHIST_OUTLEV>=PHIST_DEBUG
  PHIST_SOUT(PHIST_DEBUG,"starting function iterate() with %d systems\n curDimVs: ",numSys);
  for (int i=0;i<numSys;i++)
    {
    PHIST_SOUT(PHIST_DEBUG,"%d ",S_array[i]->curDimV_);
    }
  PHIST_SOUT(PHIST_DEBUG,"\n");
#endif
  
  // check if the given state objects have computed the norm of B,
  // and do so if not.
  for (int i=0;i<numSys;i++)
  {
    S_array[i]->ierr=1; // not converged yet
    if (S_array[i]->normB_<mt::zero())
    {
      if (S_array[i]->B_==NULL)
      {
        PHIST_OUT(PHIST_ERROR,"rhs vector not set in state %d, "
        "did you forget to call reset()? (file %s, line %d)",i,__FILE__,__LINE__);
        *ierr=-1;
        return;
      }
      SUBR(mvec_norm2)(S_array[i]->B_,&S_array[i]->normB_,&S_array[i]->ierr);
      PHIST_CHK_IERR(*ierr=S_array[i]->ierr,*ierr);
    }
  }
  
  // indices of next vectors V_j
  int v_idx[numSys];
  for (int i=0;i<numSys;i++)
  {
    v_idx[i] = S_array[i]->offsetVglob_ + S_array[i]->curDimV_;
  }
    
  // we return as soon as one system converges or reaches its
  // maximum permitted number of iterations. The decision about what to do
  // next is then left to the caller.
  bool anyConverged=false;
  bool anyFailed=false;

  while (!anyConverged && !anyFailed)
  {
    // create 'scattered view' V of all the vectors that we want to multiply our operator with
    PHIST_CHK_IERR(SUBR(mvec_view_scattered)(S_array[0]->Vglob_,&V,v_idx,numSys,ierr),*ierr);
    // get W as view of the locations V(:,j+1)
    for (int i=0;i<numSys;i++)
    {
      v_idx[i]++;
    }
    PHIST_CHK_IERR(SUBR(mvec_view_scattered)(S_array[0]->Vglob_,&W,v_idx,numSys,ierr),*ierr);
    // the indices are now already at j+1 for the next iteration.
    // We could avoid extracting V by swapping the V and W pointers,
    // but as the whole viewing business is about swapping around
    // pointers, it doesn't matter too much and may be more confusing.

    //W=A*(M\V(:,j));
    PHIST_CHK_IERR(Op->apply(st::one(),Op->A, V, st::zero(), W, ierr), *ierr);
    // TODO - maybe we could parallelize this loop by an OpenMP section
    //        so that the small stuff can be done in parallel too?
    for (int i=0;i<numSys;i++)
    {
      TYPE(gmresState_ptr) S = S_array[i];
      int j=S->curDimV_;
      
      PHIST_DEB("SYSTEM %d, ITERATE %d\n",i,j);

      // First check for failed and restarted systems.
      // In case a system was just (re-)started, we now only computed A*x0 (cf. reset 
      // function). Update this to the initial residual A*x0-b
      if (j>=S->maxBas_)
      {
        anyFailed=true;
        S->ierr=2;
        continue;
      }//TROET
      else if (j==0)
      {
        PHIST_OUT(PHIST_VERBOSE,"gmres state %d (re-)starts\n",i);
        TYPE(mvec_ptr) v0=NULL,Ax0=NULL;
        PHIST_CHK_IERR(SUBR(mvec_view_block)(S->V_,&v0,0,0,ierr),*ierr);
        PHIST_CHK_IERR(SUBR(mvec_view_block)(S->V_,&Ax0,1,1,ierr),*ierr);
        PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),S->B_,-st::one(),Ax0,ierr),*ierr);
        // set v0=b-A*x0
        PHIST_CHK_IERR(SUBR(mvec_set_block)(S->V_,Ax0,0,0,ierr),*ierr);
        // normalize, rs[0]=||r0||, v0=r0/||r0||
        PHIST_CHK_IERR(SUBR(mvec_normalize)(v0,&S->normR0_,ierr),*ierr);
        S->rs_[0]= S->normR0_;
        S->normR_= S->normR0_;
        // Note: this is a bit of a hack. As we have overwritten vector j=0, 
        // adjust the index array for the next matvec so that w1=A*v0, not A*v1
        v_idx[i]--;
      }
      else
      {
        // Arnoldi update. TODO - we could save some messages by
        // clustering the communication in orthog
              
        // view the existing basis vectors as Vprev
        PHIST_CHK_IERR(SUBR(mvec_view_block)(S->V_,&Vprev,0,j-1,ierr),*ierr);
        // view the next V vector as Vj (it currently contains A*v_{j-1}
        PHIST_CHK_IERR(SUBR(mvec_view_block)(S->V_,&Vj,j,j,ierr),*ierr);
        // view H(j,j) as R1
        PHIST_CHK_IERR(SUBR(sdMat_view_block)(S->H_,&R1,j,j,j-1,j-1,ierr),*ierr);
        // view H(1:j,j) as R2
        PHIST_CHK_IERR(SUBR(sdMat_view_block)(S->H_,&R2,0,j-1,j-1,j-1,ierr),*ierr);
#if 0 /*PHIST_OUTLEV>=PHIST_DEBUG*/
        PHIST_SOUT(PHIST_DEBUG,"SYSTEM %d, V(:,0:%d)=\n",i,j-1);
        PHIST_CHK_IERR(SUBR(mvec_print)(Vprev,ierr),*ierr);
        PHIST_SOUT(PHIST_DEBUG,"SYSTEM %d, V(:,%d)=\n",i,j);
        PHIST_CHK_IERR(SUBR(mvec_print)(Vj,ierr),*ierr);
#endif
        //orthogonalize
        PHIST_CHK_IERR(SUBR(orthog)(Vprev,Vj,R1,R2,2,ierr),*ierr);

#if 0 /*PHIST_OUTLEV>=PHIST_DEBUG*/
  PHIST_DEB("untransformed H");
  PHIST_CHK_IERR(SUBR(sdMat_print)(S->H_,ierr),*ierr);
#endif
        
        //    % update QR factorization of H

        // note: it is OK to work on raw data of serial dense matrices because 
        // they are typically not modified on an accelerator, however, we must 
        // somehow implement a check of this in the kernel interfaces (TODO)   
        // note also that we will access Hcol[j], which is strictly speaking   
        // R1, but by construction they should be aligned in memory. 
        ST *Hcol=NULL;
        lidx_t ldH; 
        PHIST_CHK_IERR(SUBR(sdMat_extract_view)(R2,&Hcol,&ldH,ierr),*ierr); 

        ST htmp;

        //    % apply previous (j-1) transformations to column j
        for (int jj=0;jj<j-1;jj++)
        {
          htmp = S->cs_[jj]*Hcol[jj] + 
                      S->sn_[jj]*Hcol[jj+1]; // H(j-1,j) in the last step, which is R1 above
          Hcol[jj+1] = -st::conj(S->sn_[jj])*Hcol[jj] + S->cs_[jj]*Hcol[jj+1];
          Hcol[jj]   = htmp;
        }

        // new Givens rotation for eliminating H(j+1,j)
        SUBR(rotg)(Hcol[j-1],Hcol[j],S->cs_[j-1],S->sn_[j-1],htmp);
        // eliminate H(j,j-1)
        Hcol[j-1] = htmp;
        Hcol[j]=st::zero();
        // apply to RHS
        htmp=S->cs_[j-1]*S->rs_[j-1];
        S->rs_[j] = -st::conj(S->sn_[j-1])*S->rs_[j-1];
        S->rs_[j-1]=htmp;

#if 0 /*PHIST_OUTLEV>=PHIST_DEBUG*/
  PHIST_DEB("transformed H");
  PHIST_CHK_IERR(SUBR(sdMat_print)(S->H_,ierr),*ierr);
#endif
        S->normR_=st::abs(S->rs_[j]);
      }// initial step or standard Arnoldi?

      MT relres=S->normR_/S->normB_;
      
      if (relres<S->tol)
      {
        S->ierr=0; // converged
        anyConverged=true;
      }
    S->curDimV_++;
    }// for all systems i
  PHIST_SOUT(PHIST_VERBOSE,"GMRES iteration states\n");
  PHIST_SOUT(PHIST_VERBOSE,"======================\n");
#if PHIST_OUTLEV>=PHIST_VERBOSE
  for (int i=0;i<numSys;i++)
  {
  PHIST_SOUT(PHIST_VERBOSE,"[%d]: %d\t%8.4e\n",i,
        S_array[i]->curDimV_,S_array[i]->normR_/S_array[i]->normB_);
  }
#endif
  PHIST_SOUT(PHIST_VERBOSE,"----------------------\n");
  }// while-loop (Arnoldi iterations)

  if (anyConverged)
  {
    *ierr=0;
  }
      
  if (anyFailed)
  {
    *ierr=1;
  }
return;
}
