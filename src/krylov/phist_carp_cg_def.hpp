// create new state objects. We just get an array of (NULL-)pointers
void SUBR(carp_cgStates_create)(TYPE(carp_cgState_ptr) state[], int numSys,
        const_map_ptr_t map, int nvec,int* ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  if (numSys==0) return;

  for (int i=0;i<numSys;i++)
  {
    state[i] = new TYPE(carp_cgState);
    state[i]->id=i;
    // set some default options
    state[i]->tol=(MT)100*mt::eps(); // default value
    state[i]->ierr=-2;// not initialized
    state[i]->maxIters_=1000; // default value
    state[i]->nvec_=nvec;
  
    PHIST_CHK_IERR(SUBR(mvec_create)(&state[i]->x0_,map,nvec,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_create)(&state[i]->q_,map,nvec,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_create)(&state[i]->r_,map,nvec,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_create)(&state[i]->p_,map,nvec,ierr),*ierr);

#ifndef IS_COMPLEX
    // separate imaginary parts of the vectors
    PHIST_CHK_IERR(SUBR(mvec_create)(&state[i]->x0i_,map,nvec,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_create)(&state[i]->qi_,map,nvec,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_create)(&state[i]->ri_,map,nvec,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_create)(&state[i]->pi_,map,nvec,ierr),*ierr);
#endif
    
    state[i]->alpha_r_ = new MT[nvec];
    state[i]->alpha_i_ = new MT[nvec];
    state[i]->beta_ =  new MT[nvec];
    state[i]->normR0_= new MT[nvec];
    state[i]->normR_= new MT[nvec];
    for (int j=0;j<nvec;j++)
    {
      state[i]->normR0_[j]=-mt::one(); // not initialized
      state[i]->normR_[j]=-mt::one(); // not initialized
    }
  }
}


//! delete cgState object
void SUBR(carp_cgStates_delete)(TYPE(carp_cgState_ptr) state[], int numSys, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  for (int i=0;i<numSys;i++)
  {
    PHIST_CHK_IERR(SUBR(mvec_delete)(state[i]->x0_,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_delete)(state[i]->q_,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_delete)(state[i]->r_,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_delete)(state[i]->p_,ierr),*ierr);
#ifndef IS_COMPLEX
    PHIST_CHK_IERR(SUBR(mvec_delete)(state[i]->x0i_,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_delete)(state[i]->qi_,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_delete)(state[i]->ri_,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_delete)(state[i]->pi_,ierr),*ierr);
#endif
    delete [] state[i]->alpha_r_;
    delete [] state[i]->alpha_i_;
    delete [] state[i]->beta_;
    delete state[i];
  }
}

// reset pcg state.
void SUBR(carp_cgState_reset)(TYPE(carp_cgState_ptr) S, 
        TYPE(const_crsMat_ptr) A,
        TYPE(const_mvec_ptr) b,
        TYPE(const_mvec_ptr) x0_r,
        TYPE(const_mvec_ptr) x0_i, 
        int *ierr)
{
#include "phist_std_typedefs.hpp"  
  ENTER_FCN(__FUNCTION__);
#if 1
  *ierr=-99;
#else
  *ierr=0;
  
  if (b==NULL && S->normR0_ == -mt::one())
  {
    PHIST_OUT(PHIST_ERROR,"on the first call to carp_cgState_reset you *must* provide the\n" 
                          "RHS vector b and sparse matrix A.");
    *ierr=-1;
    return;
  }
  else if (b!=NULL)
  {
    // new rhs -> need to recompute ||b-A*x0||
    S->b_=b;
    S->ierr = -1;
    S->totalIter = 0;
  }
  
  S->normR0_=-mt::one(); // needs to be computed in next iterate call
  // set V_0=X_0. iterate() will have to compute the residual and normalize it,
  // because the actual V_0 we want is r/||r||_2, but we can't easily apply the
  // operator to a single vector.
  PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),x0_r,st::zero(), S->x0_r_, ierr),*ierr);
#ifndef IS_COMPLEX
  PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),x0_r,st::zero(), S->x0_i_, ierr),*ierr);
#endif
  //TODO - recompute if shift or A changes
  S->nrm_ai2=nrms_ai2(A,sigma);

#endif
  return;
}


// implementation of pcg on several systems simultaneously
void SUBR(carp_cgStates_iterate)(TYPE(const_crsMat_ptr) A,
        TYPE(const_mvec_ptr) rhs,
        TYPE(carp_cgState_ptr) S_array[], int numSys,
        TYPE(mvec_ptr) X_r[], TYPE(mvec_ptr) X_i[],
        int* nIter, int* ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
#if 1
  *ierr=-99;
#else
  *ierr = 0;

  int numSys;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(X,&numSys,ierr),*ierr);

#if PHIST_OUTLEV>=PHIST_DEBUG
  PHIST_SOUT(PHIST_DEBUG,"starting function iterate() with %d systems\n curDimVs: ",numSys);
  for (int i=0;i<numSys;i++)
  {
    PHIST_SOUT(PHIST_DEBUG,"%d ",S[i]->totalIter_);
  }
  PHIST_SOUT(PHIST_DEBUG,"\n");
#endif

  int anyConverged=0;
  int anyFailed=0;

  if (numSys>1)
  { 
    *ierr=-99; //TODO: implement pipelining
  }
  
int itprint=1;
int itcheck=1;

*ierr=0;

TYPE(carp_cgState_ptr) S = s_array[0];

MT reltol2=S->tol*S->tol*S->normR0_*S->normR0_;

r=dkswp(A,sigma,B,b,x,omega,nrm_ai2)-x;
if (deflMethod==2)
  Vtr=V'*r;
  vr=E\Vtr;
  rtil=r-AV*(vr);
  z=apply_op(rtil,M);
  z=z-AV*(E\(V'*z));
  z=z+V*vr;
else
  z=apply_op(r,M);
end
p=z;

if (deflMethod==1)
  p = p - V*(AV'*p);
end

r2_new = r'*z;

alpha_old=0;
beta_old=0;

disp(sprintf('%d\t%e\t%e',0,sqrt(r2_new),sqrt(r2_new)/nrm_b));
for k=1:maxIter
  q=p-dkswp(A,sigma,B,bnul,p,omega,nrm_ai2);
  alpha = (r'*z)/(p'*q);
  x=x+alpha*p;
  if (mod(k-1,itcheck)==0)
    nrm_r = norm(A*x-sigma*x-b);
    if (mod(k-1,itprint)==0)
      disp(sprintf('%d\t%e\t%e',k,nrm_r,nrm_r/nrm_b));
    end
    relres=nrm_r/nrm_b;
    resvec=[resvec,nrm_r];
    if (nrm_r<tol*nrm_r0)
      break;
    end
  end
  r=r-alpha*q;
  if (deflMethod==2)
    Vtr=V'*r;
    vr=E\Vtr;
    rtil=r-AV*(vr);
    z=apply_op(rtil,M);
    z=z-AV*(E\(V'*z));
    z=z+V*vr;
  else
    z=apply_op(r,M);
  end
  r2_old=r2_new;
  r2_new=r'*z;
  beta=r2_new/r2_old;
  p=z+beta*p;

  relres=sqrt(r2_old)/nrm_b;
  %fprintf('\t%d\t%e\n',k,relres);
  %if (r2_old<reltol2) 
  %  break;
  %end
  } // while


  if (anyConverged > 0)
    *ierr=0;
      
  if (anyFailed > 0)
    *ierr=1;
#endif
  return;
}

