// create new state objects. We just get an array of (NULL-)pointers
void SUBR(cgStates_create)(TYPE(cgState_ptr) state[], int numSys,
        const_map_ptr_t map, int maxIters,int* ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  if (numSys==0) return;

  for (int i=0;i<numSys;i++)
  {
    state[i] = new TYPE(cgState);
    state[i]->id=i;
    // set some default options
    state[i]->tol=0.5; // typical starting tol for JaDa inner iterations...
    state[i]->ierr=-2;// not initialized
    state[i]->maxIters_=maxIters; // max num iters, will just be increased if necessary
  
    PHIST_CHK_IERR(SUBR(mvec_create)(&state[i]->x0_,map,1,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_create)(&state[i]->b_,map,1,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_create)(&state[i]->q_,map,1,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_create)(&state[i]->r_,map,1,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_create)(&state[i]->p_,map,1,ierr),*ierr);
    
    state[i]->alpha_ = new ST[maxIters];
    state[i]->beta_ = new MT[maxIters];
    state[i]->curDimV_=0;
    state[i]->normR0_=-mt::one(); // not initialized
    state[i]->normR_=-mt::one();
  }
}


//! delete cgState object
void SUBR(cgStates_delete)(TYPE(cgState_ptr) state[], int numSys, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  for (int i=0;i<numSys;i++)
  {
    PHIST_CHK_IERR(SUBR(mvec_delete)(state[i]->x0_,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_delete)(state[i]->b_,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_delete)(state[i]->q_,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_delete)(state[i]->r_,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_delete)(state[i]->p_,ierr),*ierr);
    delete [] state[i]->alpha_;
    delete [] state[i]->beta_;
    delete state[i];
  }
}

// reset pcg state.
void SUBR(cgState_reset)(TYPE(cgState_ptr) S, 
        TYPE(const_mvec_ptr) b,
        TYPE(const_mvec_ptr) x0,int *ierr)
{
#include "phist_std_typedefs.hpp"  
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  
  if (b==NULL && S->normR0_ == -mt::one())
  {
    PHIST_OUT(PHIST_ERROR,"on the first call to cgState_reset you *must* provide the\n" 
                          "RHS vector b and sparse matrix A.");
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
}


// implementation of pcg on several systems simultaneously
void SUBR(cgStates_iterate)(TYPE(const_op_ptr) op,
        TYPE(cgState_ptr) S[], TYPE(mvec_ptr) X,
        int* nIter, int* ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
  *ierr = 0;

  int numSys;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(X,&numSys,ierr),*ierr);

#if PHIST_OUTLEV>=PHIST_DEBUG
  PHIST_SOUT(PHIST_DEBUG,"starting function iterate() with %d systems\n curDimVs: ",numSys);
  for (int i=0;i<numSys;i++)
  {
    PHIST_SOUT(PHIST_DEBUG,"%d ",S[i]->curDimV_);
  }
  PHIST_SOUT(PHIST_DEBUG,"\n");
#endif

  int anyConverged=0;
  int anyFailed=0;
  
  *ierr=-99; //TODO: our own CG implementation

  if (anyConverged > 0)
    *ierr=0;
      
  if (anyFailed > 0)
    *ierr=1;
  return;
}

