// in this implementation, there is no mixed real/complex arithmetic
// simply since we don't have a kernel lib interface and we also want
// to support kernel libs like epetra or our own fortran variant, which
// do not have complex arithmetic let alone mixed functionality.

//TODO - our pseudo complex dot product and MVM are probably very inefficient 
//      compared to using proper complex vectors, we should check if the kernel 
//      lib supports mixed real/complex arithmetic and use it if possible.

//
void SUBR(private_dotProd)(TYPE(const_mvec_ptr) v, TYPE(const_mvec_ptr) vi,
                           TYPE(const_mvec_ptr) w, TYPE(const_mvec_ptr) wi,
                           int nvec, _ST_   *dots, _MT_* dotsi, int *ierr);
//
void SUBR(private_compResid)(TYPE(const_crsMat_ptr) A, int nvec, _ST_ sigma, _MT_ sigma_i,
                       TYPE(const_mvec_ptr) B,
                       TYPE(const_mvec_ptr) x, TYPE(const_mvec_ptr) xi,
                       TYPE(mvec_ptr) r,        TYPE(mvec_ptr) ri,
                       _MT_  *nrms, int *ierr);

void SUBR(private_printResid)(int it, int nvec, _ST_ const* normR, 
        _MT_ const* normR0, _MT_ const* normB);

// create new state objects. We just get an array of (NULL-)pointers
void SUBR(carp_cgStates_create)(TYPE(carp_cgState_ptr) state[], int numSys,
        _MT_ sigma_r[], _MT_ sigma_i[],
        TYPE(const_crsMat_ptr) A, int nvec,int* ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  if (numSys==0) return;
  
    // setup the CARP kernel and get the required data structures:
    MT* nrms_ai2i=NULL;
    void* aux=NULL;
    PHIST_CHK_IERR(SUBR(carp_setup)(A,numSys,sigma_r,sigma_i,
        &nrms_ai2i, &aux, ierr),*ierr);    

  const_map_ptr_t map=NULL;
  PHIST_CHK_IERR(SUBR(crsMat_get_row_map)(A,&map,ierr),*ierr);
  lidx_t nloc;
  PHIST_CHK_IERR(phist_map_get_local_length(map,&nloc,ierr),*ierr);
  
  for (int i=0;i<numSys;i++)
  {
    state[i] = new TYPE(carp_cgState);
    state[i]->id=i;
    state[i]->ierr=-2;// not initialized (need to call reset())
    
    state[i]->q_=NULL; state[i]->qi_=NULL;
    state[i]->r_=NULL; state[i]->ri_=NULL;
    state[i]->p_=NULL; state[i]->pi_=NULL;
    state[i]->z_=NULL; state[i]->zi_=NULL;
    
    state[i]->A_=A;
    state[i]->sigma_r_=sigma_r[i];
    state[i]->sigma_i_=sigma_i[i];
    state[i]->nvec_=nvec;
    
    state[i]->nrms_ai2i_=nrms_ai2i+i*nloc;
    state[i]->aux_=aux;
    state[i]->omega_=mt::one(); // relaxation parameter, for the
                                // moment just set it to 1, which
                                // gave good results for Graphene
                                // in the matlab tests.

    PHIST_CHK_IERR(SUBR(mvec_create)(&state[i]->q_,map,nvec,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_create)(&state[i]->r_,map,nvec,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_create)(&state[i]->p_,map,nvec,ierr),*ierr);
    // z is the preconditioned residual in CG, but as we don't have
    // additional preconditioning, we set z=r.
    state[i]->z_=state[i]->r_;
    
    state[i]->beta_ =  new MT[nvec];
    state[i]->alpha_ = new ST[nvec];

#ifndef IS_COMPLEX
    // separate imaginary parts of the vectors
    state[i]->alpha_i_ = new MT[nvec];
    PHIST_CHK_IERR(SUBR(mvec_create)(&state[i]->qi_,map,nvec,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_create)(&state[i]->ri_,map,nvec,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_create)(&state[i]->pi_,map,nvec,ierr),*ierr);
    state[i]->zi_=state[i]->ri_;
#else
    state[i]->qi_=NULL;
    state[i]->ri_=NULL;
    state[i]->pi_=NULL;
    state[i]->zi_=NULL;
#endif
    
    state[i]->normR0_= new MT[nvec];
    state[i]->normB_= new MT[nvec];
    state[i]->normR= new MT[nvec];
    for (int j=0;j<nvec;j++)
    {
      state[i]->normR0_[j]=-mt::one(); // not initialized
      state[i]->normB_[j]=-mt::one(); // not initialized
      state[i]->normR[j]=-mt::one(); // not initialized
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
    PHIST_CHK_IERR(SUBR(mvec_delete)(state[i]->q_,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_delete)(state[i]->r_,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_delete)(state[i]->p_,ierr),*ierr);
#ifndef IS_COMPLEX
    PHIST_CHK_IERR(SUBR(mvec_delete)(state[i]->qi_,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_delete)(state[i]->ri_,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_delete)(state[i]->pi_,ierr),*ierr);
#endif
    delete [] state[i]->alpha_;
    delete [] state[i]->alpha_i_;
    delete [] state[i]->beta_;
    delete [] state[i]->normR;
    delete [] state[i]->normR0_;
    delete [] state[i]->normB_;
    delete state[i];
  }
}

// reset pcg state. If normsB==NULL, the two-norm of 
// B is computed, otherwise it is copied from the given
// pointer (length S->nvec_). B must have the same number
// of vectors (S->nvec_).
void SUBR(carp_cgState_reset)(TYPE(carp_cgState_ptr) S,
        TYPE(const_mvec_ptr) B,
        _MT_* normsB,
        int *ierr)
{
#include "phist_std_typedefs.hpp"  
  ENTER_FCN(__FUNCTION__);
  *ierr=0;

  // new rhs -> need to recompute ||b-A*x0||
  S->b_=B;
  int nvec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(B,&nvec,ierr),*ierr);
  if (nvec!=S->nvec_)
  {
    PHIST_SOUT(PHIST_ERROR,"number of vectors in b must not change in %s",__FUNCTION__);
    *ierr=PHIST_INVALID_INPUT;
    return;
  }
  S->ierr = -1;
  S->numIter = 0;

  for (int i=0; i<nvec; i++)
  {
    S->normR0_[i]=-mt::one(); // needs to be computed in next iterate call
  }
  if (normsB==NULL)
  {
    PHIST_CHK_IERR(SUBR(mvec_norm2)(S->b_,S->normB_,ierr),*ierr);
  }
  else
  {
    for (int j=0;j<nvec;j++)
    {
      S->normB_[j]=normsB[j];
    }
  }
  return;
}

// implementation of pcg on several systems with multiple RHS each
void SUBR(carp_cgStates_iterate)(
        TYPE(carp_cgState_ptr) S_array[], int numSys,
        TYPE(mvec_ptr) X_r[], TYPE(mvec_ptr) X_i[],
        _MT_ tol, int maxIter, 
        int* ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
  *ierr = 0;

  // some internal settings  
  int itprint=1; // how often to print the current impl. residual norms
  int itcheck=1; // how often to check the actual expl. res norms. We
                  // only stop iterating if *all* expl. res norms for
                  // a given shift are below the tolerance. The impl.
                  // res norm is based on the carp operator, and I'm not
                  // sure how much sense it would make to use it as an
                  // indication of convergence.

  int numSolved=0;
  TYPE(mvec_ptr) bnul=NULL; // we can just pass in b=NULL if all entries for a carp_sweep
                            // are 0 (during CG iteration), for clarity we give it a name
/*
  const void* map;
  PHIST_CHK_IERR(SUBR(crsMat_get_range_map)(S_array[0]->A_,&map,ierr),*ierr);           
  PHIST_CHK_IERR(SUBR(mvec_create)(&bnul,map,nvec,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_put_value)(bnul,st::zero(),ierr),*ierr);
*/

////////////////////////////////////////////////////
// solve systems for one shift at a time, here we //
// could have some queuing, load balancing, extra //
// level of parallelism etc.                      //
// The multiple RHS per shift are treated with a  //
// separate Krylov space per shift, like in       //
// pgmres, but here there is no memory overhead   //
// because we're doing CG.                        //
////////////////////////////////////////////////////
  for (int ishift=0;ishift<numSys; ishift++)
  {
    // get some pointers to avoid the 'S_array[ishift]->' all the time
    TYPE(carp_cgState_ptr) S = S_array[ishift];
    TYPE(const_crsMat_ptr) A=S->A_;
    const MT sigma_r = S->sigma_r_;
    const MT sigma_i = S->sigma_i_;
    const ST sigma=sigma_r+st::cmplx_I()*sigma_i;
    int nvec=S->nvec_;
    TYPE(const_mvec_ptr) b=S->b_;
    TYPE(mvec_ptr) x=X_r[ishift];
    TYPE(mvec_ptr) xi=X_i[ishift];
    
    TYPE(mvec_ptr) r =S->r_;
    TYPE(mvec_ptr) ri =S->ri_;

    TYPE(mvec_ptr) q =S->q_;
    TYPE(mvec_ptr) qi =S->qi_;

    TYPE(mvec_ptr) p =S->p_;
    TYPE(mvec_ptr) pi =S->pi_;

    TYPE(mvec_ptr) z =S->z_;
    TYPE(mvec_ptr) zi =S->zi_;
    
    // CG (Lanczos) coefficients
    ST* alpha=S->alpha_;
    MT* alpha_i=S->alpha_i_;
    MT* beta=S->beta_;

    int numConverged=0;

#ifndef IS_COMPLEX
    // this variable indicates that while we're in
    // real arithmetic, we still have a complex shift
    // and thus xi!=NULL
    const bool rc_variant=(sigma_i!=mt::zero());
#endif
    
    if (S->ierr==-2)
    {
      *ierr=-2;
      PHIST_SOUT(PHIST_ERROR, "for carp_cgState[%d], reset has not been called!\n"
                              "(file %s, line %d)\n",S->id, __FILE__,__LINE__);
      return;
    }
    else if (S->ierr==-1)
    {
      // compute initial residual normR0 after reset
      PHIST_CHK_IERR(SUBR(private_compResid)(A, nvec, sigma, sigma_i,
                       b, x, xi, NULL, NULL, S->normR, ierr),*ierr);
      for (int j=0;j<nvec;j++)
      {
        S->normR0_[j]=std::sqrt(S->normR[j]);
      }
    }
    S->ierr=1; // unumConverged
    MT reltol2[nvec];
    for (int j=0;j<nvec;j++)
    {
      //reltol2[j]=tol*tol*S->normR[j];
      reltol2[j]=tol*tol*S->normB_[j]*S->normB_[j];
    }

    // initial Kaczmarz/CARP sweep. Note that our function carp_sweep operates
    // in place, but we do not want to update x right now, we just want to
    // get a CG direction.
    //r=carp_sweep(A,sigma,B,b,x,omega,nrm_ai2)-x;
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),x,st::zero(),r,ierr),*ierr);
#ifndef IS_COMPLEX
    if (rc_variant)
    {
      PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),xi,st::zero(),ri,ierr),*ierr);
    }
#endif
    // double carp sweep in place, updates r=dkswp(sI-A,omega,r)
    PHIST_CHK_IERR(SUBR(carp_sweep)(A, 1, &sigma_r, &sigma_i,b,&r,&ri,
          S->nrms_ai2i_,S->aux_,&S->omega_,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(-st::one(),x,st::one(),r,ierr),*ierr);
#ifndef IS_COMPLEX
    if (rc_variant)
    {
      PHIST_CHK_IERR(SUBR(mvec_add_mvec)(-st::one(),xi,st::one(),ri,ierr),*ierr);
    }
#endif
    // z=precond(r)
    // ... we currently have z pointing to r, as there is no additional preconditioning.
    //p=z
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),z,st::zero(),p,ierr),*ierr);
#ifndef IS_COMPLEX
    if (rc_variant)
    {
      PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),zi,st::zero(),pi,ierr),*ierr);
    }
#endif
    //r2_new = r'*z. For technical reasons we store (r'z)^2 as complex type if IS_COMPLEX
    ST r2_new[nvec];
    MT r2_old[nvec];
    PHIST_CHK_IERR(SUBR(private_dotProd)(r,ri,z,zi,nvec,r2_new,NULL,ierr),*ierr);

    if (itprint>0)
    {
      // print header
      SUBR(private_printResid)(0,0,NULL,NULL,NULL);
      SUBR(private_printResid)(0,nvec,r2_new,S->normR0_,S->normB_);
    }
    for (int it=1;it<maxIter; it++)
    {
      // apply operator, I-carp_sweep(...) to p. Note that our function carp_sweep operates
      // in place, so we first copy p to q again. The rhs vector is 0, which the
      // kernel lib should understand if we pass in NULL.

      //q=p-carp_sweep(A,sigma,B,bnul,p,omega,nrm_ai2);
      PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),p,st::zero(),q,ierr),*ierr);
#ifndef IS_COMPLEX
      if (rc_variant)
      {
        PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),pi,st::zero(),qi,ierr),*ierr);
      }
#endif
      // double carp sweep in place, updates q to carp_sweep(p)
      PHIST_CHK_IERR(SUBR(carp_sweep)(A, 1, &sigma_r, &sigma_i,bnul,&q,&qi,
          S->nrms_ai2i_,S->aux_,&S->omega_,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),p,-st::one(),q,ierr),*ierr);
#ifndef IS_COMPLEX
      if (rc_variant)
      {
        PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),pi,-st::one(),qi,ierr),*ierr);
      }
#endif
      ////////////////////////////
      // update solution x      //
      ////////////////////////////
      
      //alpha = (r'*z)/(p'*q);
      ST denom  [nvec];
      MT denom_i[nvec];
      PHIST_CHK_IERR(SUBR(private_dotProd)(r,ri,z,zi,nvec,alpha,alpha_i,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(private_dotProd)(p,pi,q,qi,nvec,denom,denom_i,ierr),*ierr);
fprintf(stdout,"TROET RtZ=%e %+ei\n",alpha[0],alpha_i[0]);
fprintf(stdout,"TROET PtQ=%e %+ei\n",denom[0],denom_i[0]);
      MT minus_alpha_i[nvec];
#ifdef IS_COMPLEX
        for (int j=0;j<nvec;j++)
        {
          alpha[j]=alpha[j]/denom[j];
        }
        // update x <- x + alpha*p
        PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(alpha,p,st::one(),x,ierr),*ierr);
#else
      if (!rc_variant)
      {
        for (int j=0;j<nvec;j++)
        {
          alpha[j]=alpha[j]/denom[j];
        }
        // update x <- x + alpha*p
        PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(alpha,p,st::one(),x,ierr),*ierr);
      }
      else
      {
        for (int j=0;j<nvec;j++)
        {
          CT tmp1(alpha[j],alpha_i[j]);
          CT tmp2(denom[j],denom_i[j]);
          CT tmp3=tmp1/tmp2;
          alpha[j]=ct::real(tmp3);
          alpha_i[j]    = ct::imag(tmp3);
          minus_alpha_i[j]=-ct::imag(tmp3);
        }
        // update x <- x + alpha*p
        PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(alpha,p,st::one(),x,ierr),*ierr);
        PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(minus_alpha_i,pi,st::one(),x,ierr),*ierr);
        PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(alpha,pi,st::one(),xi,ierr),*ierr);
        PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(alpha_i,p,st::one(),xi,ierr),*ierr);
      }
#endif
fprintf(stdout,"TROET alpha=%e %+ei\n",alpha[0],alpha_i[0]);
      if ( it%itcheck == 0)
      {
        PHIST_CHK_IERR(SUBR(private_compResid)(A, nvec, sigma, sigma_i,
                         b, x, xi, NULL, NULL, S->normR, ierr),*ierr);
        if ( it%itprint==0)
        {
#ifdef IS_COMPLEX
          ST tmp[nvec];
          for (int j=0;j<nvec;j++)
          {
            tmp[j]=(ST)S->normR[j];
          }
#else
          ST* tmp=S->normR;
#endif          
          SUBR(private_printResid)(it, nvec, tmp, S->normR0_, S->normB_);
        }
        
        // check for convergence. 
        // TODO - which convergence criterion should we use?
        // For the moment, we use ||r||_2/||b||_2 < tol
        numConverged=0;
        for (int j=0;j<nvec;j++)
        {
          if (S->normR[j]<reltol2[j])
          {
            numConverged++;
          }
        }
        if (numConverged==nvec)
        {
          numSolved++; 
          break;
        }
      }

      //r=r-alpha*q;
      ST minus_alpha[nvec];
      for (int j=0;j<nvec;j++)
      {
        minus_alpha[j]=-alpha[j];
      }
      PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(minus_alpha,q,st::one(),r,ierr),*ierr);
#ifndef IS_COMPLEX
      if (rc_variant&&false)
      {
        PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(      alpha_i,qi,st::one(),r,ierr),*ierr);
        PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(minus_alpha,  qi,st::one(),ri,ierr),*ierr);
        PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(minus_alpha_i,q, st::one(),ri,ierr),*ierr);
      }
#endif
//  z=apply_op(r,M):
// .. do nothing ...

      for (int j=0;j<nvec;j++)
      {
        r2_old[j]=std::abs(r2_new[j]);
      }
      PHIST_CHK_IERR(SUBR(private_dotProd)(r,ri,z,zi,nvec,r2_new,NULL,ierr),*ierr);
      //note: for a precond M!=I I think we need complex
      // arithmetic here because r!=z => above dotProd has imag!=0.
      // Only for symmetric preconditioning would we get a Hermitian
      // Lanczos matrix and thus a real beta on the diagonal.
      for (int j=0;j<nvec;j++)
      {
        beta[j]=std::abs(r2_new[j])/r2_old[j];
        fprintf(stdout,"TROET beta nom  =%e\n",r2_new[j]);
        fprintf(stdout,"TROET beta denom=%e\n",r2_old[j]);
        fprintf(stdout,"TROET beta=%e\n",beta[j]);
      }
      //p=z+beta*p;
#ifdef IS_COMPLEX
      ST tmp[nvec];
      for (int j=0;j<nvec;j++)
      {
        tmp[j]=(ST)beta[j];
      }
      PHIST_CHK_IERR(SUBR(mvec_vscale)(p,tmp,ierr),*ierr);
#else
      PHIST_CHK_IERR(SUBR(mvec_vscale)(p,beta,ierr),*ierr);
      if (rc_variant)
      {
        PHIST_CHK_IERR(SUBR(mvec_vscale)(pi,beta,ierr),*ierr);
      }
#endif
      PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),z,st::one(),p,ierr),*ierr);
#ifndef IS_COMPLEX
      if (rc_variant)
      {
        PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),zi,st::one(),pi,ierr),*ierr);
      }
#endif
    }// CG iterations (it)
  
    //...
    
    // take square root so that normR is indeed the norm, not the
    // squared 2-norm:
    for (int j=0;j<nvec;j++)
    {
      S->normR[j]=std::sqrt(S->normR[j]);
    }
    
  } // for all shifts, solve (s[j]I-A)X[j]=B


  *ierr=numSys-numSolved;
      
  return;
}

// compute residual r=b-(sI-A)x and ||r||_2^2 in nrms2. If r and ri are NULL, a temporary
// vector is used and discarded.
void SUBR(private_compResid)(TYPE(const_crsMat_ptr) A, int nvec, _ST_ sigma, _MT_ sigma_i,
                       TYPE(const_mvec_ptr) B,
                       TYPE(const_mvec_ptr) x, TYPE(const_mvec_ptr) xi,
                       TYPE(mvec_ptr) r,        TYPE(mvec_ptr) ri,
                       _MT_  *nrms2, int *ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);

#ifndef IS_COMPLEX
  bool rc_variant= (sigma_i!=mt::zero()) ||
                 (xi!=NULL);
#else
  bool rc_variant=false;
#endif

  ST shifts[nvec];
  for (int i=0;i<nvec;i++)
  {
    shifts[i]=-sigma;
  }

  TYPE(mvec_ptr) R=r, RI=ri;
  if (R==NULL)
  {
    const_map_ptr_t map;
    PHIST_CHK_IERR(SUBR(mvec_get_map)(x,&map,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_create)(&R,map,nvec,ierr),*ierr);
  }
  if (RI==NULL)
  {
    const_map_ptr_t map;
    PHIST_CHK_IERR(SUBR(mvec_get_map)(xi,&map,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_create)(&RI,map,nvec,ierr),*ierr);
  }

  // r = b-(sI-A)x=b+(A-sI)x
  // r = b+(A-sr)xr + si*xi
  //    +i[(A-sr)xi - si*xr]
  
  // r=b
  PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),B,st::zero(),R,ierr),*ierr);
  // r=b-(sI-A)*x
  PHIST_CHK_IERR(SUBR(crsMat_times_mvec_vadd_mvec)
    (st::one(),A,shifts,x,st::one(),R,ierr),*ierr);
  if (rc_variant)
  {
    //r=r+si*xi
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(sigma_i, xi,st::one(),R,ierr),*ierr);
    //ri=(A-srI)xi
    PHIST_CHK_IERR(SUBR(crsMat_times_mvec_vadd_mvec)
      (st::one(),A,shifts,xi,st::zero(),RI,ierr),*ierr);
    // ri-=si*xr
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(-sigma_i,x,st::one(),RI,ierr),*ierr);
  }
  // now compute the 2-norm of each column
#ifdef IS_COMPLEX
          ST tmp[nvec];
#else
          ST* tmp=nrms2;
#endif          
  PHIST_CHK_IERR(SUBR(private_dotProd)(R, RI, R, RI, nvec, tmp, NULL, ierr),*ierr);
#ifdef IS_COMPLEX
  for (int j=0;j<nvec;j++)
  {
    nrms2[j]=st::real(tmp[j]);
  }
#endif
  if (R!=r && R!=NULL)
  {
    PHIST_CHK_IERR(SUBR(mvec_delete)(R,ierr),*ierr);
  }
  if (RI!=ri && RI!=NULL)
  {
    PHIST_CHK_IERR(SUBR(mvec_delete)(RI,ierr),*ierr);
  }
  return;
}

// compute v'*w, where v and w may have an imaginary part even in real (S/D) case.
// in the complex case (C/Z), vi and wi are ignored, the result dots is complex and
// dotsi is ignored. In the real case, if vi and wi are NULL, they are assumed
// to be 0 and dotsi is ignored.
void SUBR(private_dotProd)(TYPE(const_mvec_ptr) v, TYPE(const_mvec_ptr) vi,
                           TYPE(const_mvec_ptr) w, TYPE(const_mvec_ptr) wi,
                           int nvec, _ST_ *dots, _MT_ *dotsi, int *ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
#ifndef IS_COMPLEX
  // a NULL pointer for vi or wi is interpreted as the
  // imaginary part of v (w) being 0, which means there
  // is no contribution from that part.
  const bool rc_variant = (vi!=NULL || wi!=NULL);
#endif
  *ierr=0;
  PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(v,w,dots,ierr),*ierr);
  #ifndef IS_COMPLEX
  if (rc_variant)
  {
    ST tmp[nvec];
    PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(vi,wi,tmp,ierr),*ierr);
    for (int j=0;j<nvec;j++)
    {
      dots[j]+=tmp[j];
    }
    if ( (v!=w || vi!=wi) && dotsi!=NULL)
    {
      PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(v,wi,dotsi,ierr),*ierr);      
      PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(vi,w,tmp,ierr),*ierr);      
      for (int j=0;j<nvec;j++)
      {
        dotsi[j]-=tmp[j];
      }
    }
    else if (dotsi!=NULL)
    {
      for (int j=0;j<nvec;j++)
      {
        dotsi[j]=mt::zero();
      }
    }
  }
#endif
  return;
}

// print residual info. If we get nvec=0, the input arrays are ignored
// and a header is printed.
void SUBR(private_printResid)(int it, int nvec, _ST_ const* normR, 
        _MT_ const* normR0, _MT_ const* normB)
{
#include "phist_std_typedefs.hpp"
  const char* carp_label = "CARP_CG";
  if (nvec==0)
  {
    PHIST_SOUT(PHIST_INFO,"%s\tit\t||r||/||b||\t||r||/||r0||\n",carp_label);
  }
  else
  {
    ST tmp=mt::sqrt(st::real(normR[0]));
    PHIST_SOUT(PHIST_INFO,"%s %d\t%e\t%e\n",carp_label,it,
          tmp/normB[0],tmp/normR0[0])
    for (int j=1;j<nvec;j++)
    {
      MT tmp=mt::sqrt(st::real(normR[j]));
      PHIST_SOUT(PHIST_INFO,"%s\t\t%e\t%e\n",carp_label,
             tmp/normB[j],tmp/normR0[j]);
    }
  }
}
