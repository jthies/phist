// in this implementation, there is no mixed real/complex arithmetic
// simply since we don't have a kernel lib interface and we also want
// to support kernel libs like epetra or our own fortran variant, which
// do not have complex arithmetic let alone mixed functionality.

//TODO - our pseudo complex dot product and MVM are probably very inefficient 
//      compared to using proper complex vectors, we should check if the kernel 
//      lib supports mixed real/complex arithmetic and use it if possible.

///////////////////////////////////////////////////////////////////////////////////////////////////
// private helper functions                                                                      //
///////////////////////////////////////////////////////////////////////////////////////////////////

//
void SUBR(private_dotProd)(TYPE(const_mvec_ptr) v, TYPE(const_mvec_ptr) vi, TYPE(const_sdMat_ptr) vp, TYPE(const_sdMat_ptr) vpi,
                           TYPE(const_mvec_ptr) w, TYPE(const_mvec_ptr) wi, TYPE(const_sdMat_ptr) wp, TYPE(const_sdMat_ptr) wpi,
                           int nvec, _ST_   *dots, _MT_* dotsi, bool rc_variant, int *iflag);
//
void SUBR(private_compResid)(TYPE(const_sparseMat_ptr) A, int nvec, _MT_ const sigma_r[], _MT_ const sigma_i[],
                       TYPE(const_mvec_ptr) Vproj,
                       TYPE(const_mvec_ptr) Rhs,
                       TYPE(const_mvec_ptr) x, TYPE(const_mvec_ptr) xi, TYPE(const_sdMat_ptr) xp, TYPE(const_sdMat_ptr) xpi,
                       TYPE(mvec_ptr) r,        TYPE(mvec_ptr) ri, TYPE(sdMat_ptr) rp, TYPE(sdMat_ptr) rpi,
                       _MT_  *nrms, bool rc_variant, int *iflag);

// pretty-print convergence history. nvec=0: print header. nvec=-1: print footer
void SUBR(private_printResid)(int it, int nvec, _ST_ const* normR, 
        _MT_ const* normR0, _MT_ const* normB, int const* locked);

// allocate CG vector blocks
void SUBR(private_carp_cgState_alloc)(TYPE(carp_cgState_ptr) S, int* iflag);
// allocate CG vector blocks
void SUBR(private_carp_cgState_dealloc)(TYPE(carp_cgState_ptr) S, int* iflag);

// destroy vector elements to test automatic fault-detection
void SUBR(destroy_vector_elements_at_random)(TYPE(mvec_ptr) V, double probability, int nloc, int *iflag); 
void SUBR(show_vector_elements)(TYPE(mvec_ptr) V, int *iflag) ;

///////////////////////////////////////////////////////////////////////////////////////////////////
// public interface                                                                              //
///////////////////////////////////////////////////////////////////////////////////////////////////

// create new state object
void SUBR(carp_cgState_create)(TYPE(carp_cgState_ptr) *state,
        TYPE(const_sparseMat_ptr) A, TYPE(const_mvec_ptr) Vproj,
        int nvec, _MT_ sigma_r[], _MT_ sigma_i[],
        int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
  if (nvec==0) return;
  
    // setup the CARP kernel and get the required data structures:
    void* aux=NULL;
    PHIST_CHK_IERR(SUBR(carp_setup)(A,nvec,sigma_r,sigma_i,
        &aux, iflag),*iflag);    

  const_map_ptr_t map=NULL;
  PHIST_CHK_IERR(SUBR(sparseMat_get_row_map)(A,&map,iflag),*iflag);
  lidx_t nloc;
  PHIST_CHK_IERR(phist_map_get_local_length(map,&nloc,iflag),*iflag);
  
  *state = new TYPE(carp_cgState);
  (*state)->iflag=-2;// not initialized (need to call reset())
  
  (*state)->rc_variant_=false;
#ifndef IS_COMPLEX
  for (int i=0; i<nvec;i++)
  {
    if (sigma_i[i]!=mt::zero()) (*state)->rc_variant_=true;
  }
#endif
    
  (*state)->q_=NULL; (*state)->qi_=NULL;
    (*state)->r_=NULL; (*state)->ri_=NULL;
    (*state)->p_=NULL; (*state)->pi_=NULL;
    (*state)->z_=NULL; (*state)->zi_=NULL;
    
    (*state)->A_=A;
    (*state)->Vproj_=Vproj;

    (*state)->conv = new int[nvec];
    (*state)->sigma_r_ = new MT[nvec];
    (*state)->sigma_i_ = new MT[nvec];
    
    (*state)->normR0_= new MT[nvec];
    (*state)->normB_= new MT[nvec];
    (*state)->normR= new MT[nvec];
    (*state)->normR_old= new MT[nvec];

    (*state)->beta_ =  new MT[nvec];
    (*state)->alpha_ = new ST[nvec];
    (*state)->alpha_i_ = new MT[nvec];
    (*state)->omega_=new MT[nvec]; // relaxation parameter for col j

    for (int i=0; i<nvec; i++)
    {
      (*state)->sigma_r_[i]=sigma_r[i];
      (*state)->sigma_i_[i]=sigma_i[i];
      (*state)->omega_[i]=mt::one();
    }

    (*state)->nvec_=nvec;

    (*state)->aux_=aux;
                                // moment just set it to 1, which
                                // gave good results for Graphene
                                // in the matlab tests.

    for (int j=0;j<nvec;j++)
    {
      (*state)->conv[j]=0;
      (*state)->normR0_[j]=-mt::one(); // not initialized
      (*state)->normB_[j]=-mt::one(); // not initialized
      (*state)->normR[j]=-mt::one(); // not initialized
    }
}

void SUBR(private_carp_cgState_alloc)(TYPE(carp_cgState_ptr) S, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  
    int nvec=S->nvec_;
    const void* map=NULL;
    PHIST_CHK_IERR(SUBR(mvec_get_map)(S->b_,&map,iflag),*iflag);
  
    PHIST_CHK_IERR(SUBR(mvec_create)(&S->q_,map,nvec,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_create)(&S->r_,map,nvec,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_create)(&S->p_,map,nvec,iflag),*iflag);
    // z is the preconditioned residual in CG, but as we don't have
    // additional preconditioning, we set z=r.
    S->z_=S->r_;    

    if (S->rc_variant_)
    {
      // separate imaginary parts of the vectors
      PHIST_CHK_IERR(SUBR(mvec_create)(&S->qi_,map,nvec,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(mvec_create)(&S->ri_,map,nvec,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(mvec_create)(&S->pi_,map,nvec,iflag),*iflag);
      S->zi_=S->ri_;
    }
    else
    {
      S->qi_=NULL;
      S->ri_=NULL;
      S->pi_=NULL;
      S->zi_=NULL;
    }
}

void SUBR(private_carp_cgState_dealloc)(TYPE(carp_cgState_ptr) S, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  
    PHIST_CHK_IERR(SUBR(mvec_delete)(S->q_,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_delete)(S->r_,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_delete)(S->p_,iflag),*iflag);
    // z is the preconditioned residual in CG, but as we don't have
    // additional preconditioning, we set z=r.
    
#ifndef IS_COMPLEX
    // separate imaginary parts of the vectors
    PHIST_CHK_IERR(SUBR(mvec_delete)(S->qi_,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_delete)(S->ri_,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_delete)(S->pi_,iflag),*iflag);
#endif
    
}

//! delete cgState object
void SUBR(carp_cgState_delete)(TYPE(carp_cgState_ptr) state, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
  
  PHIST_CHK_IERR(SUBR(carp_destroy)(state->A_,
        state->aux_, iflag),*iflag);
  
  
    delete []  state->conv;

    delete [] state->normR0_;
    delete [] state->normB_;
    delete [] state->normR;
    delete [] state->normR_old;

    delete [] state->sigma_r_;
    delete [] state->sigma_i_;
    delete [] state->omega_;

    delete [] state->beta_;
    delete [] state->alpha_;
    delete [] state->alpha_i_;

    delete state;
}

// reset pcg state. If normsB==NULL, the two-norm of 
// B is computed, otherwise it is copied from the given
// pointer (length num_vectors of B).
void SUBR(carp_cgState_reset)(TYPE(carp_cgState_ptr) S,
        TYPE(const_mvec_ptr) B,
        _MT_* normsB,
        int *iflag)
{
#include "phist_std_typedefs.hpp"  
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;

  // new rhs -> need to recompute ||b-A*x0||
  S->b_=B;
  int nvec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(B,&nvec,iflag),*iflag);
  S->iflag = -1;
  S->numIter = 0;

  if (nvec!=S->nvec_)
  {
    *iflag=PHIST_NOT_IMPLEMENTED;
    PHIST_SOUT(PHIST_ERROR,"cannot reset with different #vectors in the current implementation,\n"
                           "you will have to destroy the object and create it anew.\n"
                           "(file %s, line %d)\n",__FILE__,__LINE__);
    return;
/*  
    delete [] S->conv;
    S->conv= new int[nvec];

    delete [] S->normR0_;
    S->normR0_= new MT[nvec];
    
    delete [] S->normB_;
    S->normB_= new MT[nvec];

    delete [] S->normR;
    S->normR= new MT[nvec];

    delete [] S->normR_old;
    S->normR_old= new MT[nvec];

    delete [] S->beta_;
    S->beta_ =  new MT[nvec];
    delete [] S->alpha_;
    S->alpha_ = new ST[nvec];
    delete [] S->alpha_i_;
    S->alpha_i_ = new MT[nvec];
  */
  }// reallocate scalars if nvec changes

  for (int i=0; i<nvec; i++)
  {
    S->conv[i]=0;
    S->normR0_[i]=-mt::one(); // needs to be computed in next iterate call
  }
  if (normsB==NULL)
  {
    PHIST_CHK_IERR(SUBR(mvec_norm2)(S->b_,S->normB_,iflag),*iflag);
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

// implementation of pcg on several systems with multiple RHS each.
//
void SUBR(carp_cgState_iterate)(
        TYPE(carp_cgState_ptr) S,
        TYPE(mvec_ptr) X_r, TYPE(mvec_ptr) X_i,
        _MT_ tol, int maxIter, bool abortIfOneConverges,
        int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag = 0;
  //getchar();
  // some internal settings 
  bool correction_needed = false; 
  int cor_count = 0;
  double cor_tol = 0.99;
//  double cor_tol = 0.999999999999999999999e99;
  int debugstate = 999999999;
  //std::cout << tol << std:: endl;

  // how often to print the current impl. residual norms
#if PHIST_OUTLEV>=PHIST_VERBOSE
  int itprint=1; 
  int itcheck=1;
#else
  int itprint=-1; 
  int itcheck=10; // how often to check the actual expl. res norms. We
                  // only stop iterating if *all* expl. res norms for
                  // a given shift are below the tolerance. The impl.
                  // res norm is based on the carp operator, and I'm not
                  // sure how much sense it would make to use it as an
                  // indication of convergence.
#endif
  itcheck=std::min(itcheck,maxIter);
  int numSolved=0;
  TYPE(mvec_ptr) bnul=NULL; // we can just pass in b=NULL if all entries for a carp_sweep
                            // are 0 (during CG iteration), for clarity we give it a name

    // get some pointers to ease the notation
    TYPE(const_sparseMat_ptr) A=S->A_;
    const MT* sigma_r = S->sigma_r_;
    const MT* sigma_i = S->sigma_i_;
    int nvec=S->nvec_;
    TYPE(const_mvec_ptr) b=S->b_;
    TYPE(mvec_ptr) x=X_r;
    TYPE(mvec_ptr) xi=X_i;

  // used for premature termination of the loop if requested
  int minConv=abortIfOneConverges? 1: nvec;
   
    int nvecX,nvecB;
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(b,&nvecB,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(x,&nvecX,iflag),*iflag);
    
    if (nvec!=nvecB || nvec!=nvecX)
    {
      PHIST_SOUT(PHIST_ERROR,"input vectors to %s must have same num vectors as block size \n"
                             "passed to constructor (expected %d, found nvec(X)=%d, nvec(B)=%d instead)\n"
                             "(in %s, line %d)\n",
        __FUNCTION__,nvec,nvecX,nvecB,__FILE__,__LINE__);
      *iflag=-1;
      return;      
    }
#ifndef IS_COMPLEX
    if (xi!=NULL)
    {
      PHIST_CHK_IERR(SUBR(mvec_num_vectors)(xi,&nvecX,iflag),*iflag);
      if (nvec!=nvecX)
      {
        PHIST_SOUT(PHIST_ERROR,"input vectors X_i to %s must have same num vectors as X_r and B\n",__FUNCTION__);
        *iflag=-1;
        return;
      }
    }
#endif    

    // allocate CG vectors: one per rhs
    PHIST_CHK_IERR(SUBR(private_carp_cgState_alloc)(S,iflag),*iflag);

    TYPE(mvec_ptr) r =S->r_;
    TYPE(mvec_ptr) ri =S->ri_;

    TYPE(mvec_ptr) q =S->q_;
    TYPE(mvec_ptr) qi =S->qi_;

    TYPE(mvec_ptr) p =S->p_;
    TYPE(mvec_ptr) pi =S->pi_;

    TYPE(mvec_ptr) z =S->z_;
    TYPE(mvec_ptr) zi =S->zi_;
    
    // for additionally projected CARP-CG
    TYPE(const_mvec_ptr) Vproj=S->Vproj_;
    TYPE(sdMat_ptr) rp=S->rp_;
    TYPE(sdMat_ptr) rpi=S->rpi_;

    TYPE(sdMat_ptr) qp=S->qp_;
    TYPE(sdMat_ptr) qpi=S->qpi_;

    TYPE(sdMat_ptr) zp=S->zp_;
    TYPE(sdMat_ptr) zpi=S->zpi_;
    
    // CG (Lanczos) coefficients
    ST* alpha=S->alpha_;
    MT* alpha_i=S->alpha_i_;
    MT* beta=S->beta_;
    
    int *conv=S->conv;

    int numConverged=0;

    if (S->iflag==-2)
    {
      *iflag=-2;
      PHIST_SOUT(PHIST_ERROR, "for carp_cgState[%d], reset has not been called!\n"
                              "(file %s, line %d)\n",S->id, __FILE__,__LINE__);
      return;
    }
    else if (S->iflag==-1)
    {
      // compute initial residual normR0 after reset
      PHIST_CHK_IERR(SUBR(private_compResid)(A, nvec, sigma_r, sigma_i,
                       b, x, xi, NULL, NULL, S->normR, S->rc_variant_,iflag),*iflag);
      for (int j=0;j<nvec;j++)
      {
        S->normR0_[j]=std::sqrt(S->normR[j]);
        S->normR_old[j] = S->normR0_[j];
      }
    }
    S->iflag=1; // unumConverged
    MT reltol2[nvec];

    for (int j=0;j<nvec;j++)
    {
      reltol2[j]=tol*tol*S->normR[j];
      reltol2[j]=std::max(reltol2[j],tol*tol*S->normB_[j]*S->normB_[j]);
    }

    // initial Kaczmarz/CARP sweep. Note that our function carp_sweep operates
    // in place, but we do not want to update x right now, we just want to
    // get a CG direction.
    //r=carp_sweep(A,sigma,B,b,x,omega)-x;
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),x,st::zero(),r,iflag),*iflag);
#ifndef IS_COMPLEX
    if (S->rc_variant_)
    {
      PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),xi,st::zero(),ri,iflag),*iflag);
    }
#endif
    // double carp sweep in place, updates r=dkswp(sI-A,omega,r)
    if (Vproj!=NULL)
    {
      PHIST_CHK_IERR(SUBR(carp_sweep_aug)(A, sigma_r, sigma_i,Vproj,b,r,ri,
            rp,rpi,S->aux_,S->omega_,iflag),*iflag);
    }
    else
    {
      PHIST_CHK_IERR(SUBR(carp_sweep)(A, sigma_r, sigma_i,b,r,ri,
            S->aux_,S->omega_,iflag),*iflag);
    }
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(-st::one(),x,st::one(),r,iflag),*iflag);
#ifndef IS_COMPLEX
    if (S->rc_variant_)
    {
      PHIST_CHK_IERR(SUBR(mvec_add_mvec)(-st::one(),xi,st::one(),ri,iflag),*iflag);
    }
#endif
    // z=precond(r)
    // ... we currently have z pointing to r, as there is no additional preconditioning.
    //p=z
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),z,st::zero(),p,iflag),*iflag);
#ifndef IS_COMPLEX
    if (S->rc_variant_)
    {
      PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),zi,st::zero(),pi,iflag),*iflag);
    }
#endif
    //r2_new = r'*z. For technical reasons we store (r'z)^2 as complex type if IS_COMPLEX
    ST r2_new[nvec];
    MT r2_old[nvec];
    PHIST_CHK_IERR(SUBR(private_dotProd)(r,ri,z,zi,nvec,r2_new,NULL,S->rc_variant_,iflag),*iflag);

    if (itprint>0)
    {
      // print header
      SUBR(private_printResid)(0,0,NULL,NULL,NULL,NULL);
      SUBR(private_printResid)(0,nvec,r2_new,S->normR0_,S->normB_,S->conv);
    }

    //getchar();
    //PHIST_CHK_IERR(SUBR(show_vector_elements)(r, iflag),*iflag);
    //PHIST_CHK_IERR(SUBR(show_vector_elements)(p, iflag),*iflag);
    //PHIST_CHK_IERR(SUBR(show_vector_elements)(x, iflag),*iflag);

    //PHIST_SOUT(PHIST_INFO,"HIER FAENGT DER EIGENTLICHE CARP-CG TEIL AN UND ICH BRAUCHE EINE AUSGABE DAMIT ICH DIESEN TEIL WIEDERFINDEN KANN.\n");                                                            ///////////////////////////////////////////////////////////////////
    //maxIter =10000;

    for (int it=1;it<maxIter; it++)
    {    
      //if(it % 100 == 0){
      // getchar();  
      //}
      
      if(it>debugstate){
      std::cout << "Initial vectors. " << std::endl;
      getchar();
      PHIST_CHK_IERR(SUBR(show_vector_elements)(r, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(show_vector_elements)(p, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(show_vector_elements)(q, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(show_vector_elements)(x, iflag),*iflag);
      }

      /* This code implements the self-stabilizing CG described in Sao & Vuduc, ScalA'14
      proceedings. The implementation was done by Florian Fritzen in an internship and seems
      way too invasive to me, this if statement should be reduced to the few lines it is
      actually supposed to touch (TODO).
      */
      if(correction_needed == true ){
        PHIST_SOUT(PHIST_INFO,"CARP_CG: correction step\n");
        correction_needed = false;
        cor_count++;
      //q=p-carp_sweep(A,sigma,B,bnul,p,omega);
      PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),p,st::zero(),q,iflag),*iflag);
#ifndef IS_COMPLEX
      if (S->rc_variant_)
      {
        PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),pi,st::zero(),qi,iflag),*iflag);
      }
#endif
      // double carp sweep in place, updates q to carp_sweep(p)
      PHIST_CHK_IERR(SUBR(carp_sweep)(A, sigma_r, sigma_i,bnul,q,qi,
          S->aux_,S->omega_,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),p,-st::one(),q,iflag),*iflag);
#ifndef IS_COMPLEX
      if (S->rc_variant_)
      {
        PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),pi,-st::one(),qi,iflag),*iflag);
      }
#endif
      
      if(it>debugstate){
      std::cout << "Nach Q-Sweep. " << std::endl;
      getchar();
      PHIST_CHK_IERR(SUBR(show_vector_elements)(r, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(show_vector_elements)(p, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(show_vector_elements)(q, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(show_vector_elements)(x, iflag),*iflag);
      }
      

//r=carp_sweep(A,sigma,B,b,x,omega)-x;
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),x,st::zero(),r,iflag),*iflag);
#ifndef IS_COMPLEX
    if (S->rc_variant_)
    {
      PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),xi,st::zero(),ri,iflag),*iflag);
    }
#endif
    // double carp sweep in place, updates r=dkswp(sI-A,omega,r)
    PHIST_CHK_IERR(SUBR(carp_sweep)(A, sigma_r, sigma_i,b,r,ri,
          S->aux_,S->omega_,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(-st::one(),x,st::one(),r,iflag),*iflag);
#ifndef IS_COMPLEX
    if (S->rc_variant_)
    {
      PHIST_CHK_IERR(SUBR(mvec_add_mvec)(-st::one(),xi,st::one(),ri,iflag),*iflag);
    }
#endif
      
      if(it>debugstate){
      std::cout << "Nach RC-Sweep. " << std::endl;
      getchar();
      PHIST_CHK_IERR(SUBR(show_vector_elements)(r, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(show_vector_elements)(p, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(show_vector_elements)(q, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(show_vector_elements)(x, iflag),*iflag);
      }

      ////////////////////////////
      // update solution x      //
      ////////////////////////////
      
      //alpha = (r'*z)/(p'*q);
      ST denom  [nvec];
      MT denom_i[nvec];
      PHIST_CHK_IERR(SUBR(private_dotProd)(r,ri,p,pi,nvec,alpha,alpha_i,S->rc_variant_,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(private_dotProd)(p,pi,q,qi,nvec,denom,denom_i,S->rc_variant_,iflag),*iflag);
      MT minus_alpha_i[nvec];

      // stop updating x if the system is already converged (1-conv[j])
#ifdef IS_COMPLEX
        for (int j=0;j<nvec;j++)
        {
          alpha[j]=(ST)(1-conv[j])*(alpha[j]/denom[j]);
        }
        // update x <- x + alpha*p
        PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(alpha,p,st::one(),x,iflag),*iflag);
#else
      if (!S->rc_variant_)
      {
        for (int j=0;j<nvec;j++)
        {
          alpha[j]=(ST)(1-conv[j])*(alpha[j]/denom[j]);
        }
        // update x <- x + alpha*p
        PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(alpha,p,st::one(),x,iflag),*iflag);
      }
      else
      {
        for (int j=0;j<nvec;j++)
        {
          CT tmp1(alpha[j],alpha_i[j]);
          CT tmp2(denom[j],denom_i[j]);
          CT tmp3=tmp1/tmp2;
          tmp1*=(MT)(1-conv[j]);
          alpha[j]=ct::real(tmp3);
          alpha_i[j]    = ct::imag(tmp3);
          minus_alpha_i[j]=-ct::imag(tmp3);
        }
        // update x <- x + alpha*p
        PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(alpha,p,st::one(),x,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(minus_alpha_i,pi,st::one(),x,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(alpha,pi,st::one(),xi,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(alpha_i,p,st::one(),xi,iflag),*iflag);
      }
#endif
      
      if(it>debugstate){
      std::cout << "Nach X-Sweep. " << std::endl;
      getchar();
      PHIST_CHK_IERR(SUBR(show_vector_elements)(r, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(show_vector_elements)(p, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(show_vector_elements)(q, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(show_vector_elements)(x, iflag),*iflag);
      }
      

      if ( it%itcheck == 0)
      {  
        for (int j=0;j<nvec;j++)
        {
          S->normR_old[j] = (S->normR[j]);
        }
        PHIST_CHK_IERR(SUBR(private_compResid)(A, nvec, sigma_r, sigma_i,
                         b, x, xi, NULL, NULL, S->normR, S->rc_variant_,iflag),*iflag);
        //std::cout << *(S->normR) <<" " << *(S->normR_old)<< " "<< (*(S->normR)/(*(S->normR_old))) << std::endl;
        // check for convergence. 
        // TODO - which convergence criterion should we use?
        // For the moment, we use ||r||_2/||b||_2 < tol
        numConverged=0;
        for (int j=0;j<nvec;j++)
        {
          if (S->normR[j]<reltol2[j])
          {
            conv[j]=1;
            numConverged++;
          }
        }

        if ( (it%itprint==0 && itprint>0) || (numConverged>=minConv))
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
          SUBR(private_printResid)(it, nvec, tmp, S->normR0_, S->normB_,S->conv);  
                                                                  ///////////////////////////////////////////////////////////////////
          if(std::sqrt(*(S->normR))/std::sqrt(*(S->normR_old)) > cor_tol){
            correction_needed = true;
            //std::cout <<  std::sqrt(*(S->normR)) << " " << std::sqrt((*(S->normR_old))) << " "<< std::sqrt(*(S->normR))/std::sqrt((*(S->normR_old))) <<std::endl;
            //std::cout << "Cstep needed, because " << std::sqrt(*(S->normR))/std::sqrt((*(S->normR_old))) <<" > " << cor_tol <<std::endl;
          }
        }
        
        if (numConverged>=minConv)
        {
          // print footer
          SUBR(private_printResid)(it,-1,NULL,NULL,NULL,NULL);
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
      PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(minus_alpha,q,st::one(),r,iflag),*iflag);

#ifndef IS_COMPLEX
      if (S->rc_variant_)
      {
        PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(      alpha_i,qi,st::one(),r,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(minus_alpha,  qi,st::one(),ri,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(minus_alpha_i,q, st::one(),ri,iflag),*iflag);
      }
#endif
      
      if(it>debugstate){
      std::cout << "Nach R-Sweep. " << *alpha << std::endl;
      getchar();
      PHIST_CHK_IERR(SUBR(show_vector_elements)(r, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(show_vector_elements)(p, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(show_vector_elements)(q, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(show_vector_elements)(x, iflag),*iflag);
    }
      
//  z=apply_op(r,M):
// .. do nothing ...

      for (int j=0;j<nvec;j++)
      {
        r2_old[j]=std::abs(r2_new[j]);
      }
      PHIST_CHK_IERR(SUBR(private_dotProd)(r,ri,z,zi,nvec,r2_new,NULL,S->rc_variant_,iflag),*iflag);
      //note: for a precond M!=I I think we need complex
      // arithmetic here because r!=z => above dotProd has imag!=0.
      // Only for symmetric preconditioning would we get a Hermitian
      // Lanczos matrix and thus a real beta on the diagonal.

      ST upper  [nvec];
      MT upper_i[nvec];
      ST lower  [nvec];
      MT lower_i[nvec];

      PHIST_CHK_IERR(SUBR(private_dotProd)(r,ri,q,qi,nvec,upper,upper_i,S->rc_variant_,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(private_dotProd)(p,pi,q,qi,nvec,lower,lower_i,S->rc_variant_,iflag),*iflag);

      for (int j=0;j<nvec;j++)
        {
#ifndef IS_COMPLEX
          CT tmp1(upper[j],upper_i[j]);
          CT tmp2(lower[j],lower_i[j]);
          CT tmp3=-tmp1/tmp2;
          //tmp1*=(MT)(1-conv[j]);
          beta[j]=ct::real(tmp3);
#else
          beta[j]=ct::real(upper[j]/lower[j]);
#endif       
        }

      //p=z+beta*p;
#ifdef IS_COMPLEX
      ST tmp[nvec];
      for (int j=0;j<nvec;j++)
      {
        tmp[j]=(ST)beta[j];
      }
      PHIST_CHK_IERR(SUBR(mvec_vscale)(p,tmp,iflag),*iflag);
#else
      PHIST_CHK_IERR(SUBR(mvec_vscale)(p,beta,iflag),*iflag);
      if (S->rc_variant_)
      {
        PHIST_CHK_IERR(SUBR(mvec_vscale)(pi,beta,iflag),*iflag);
      }
#endif
      PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),z,st::one(),p,iflag),*iflag);
#ifndef IS_COMPLEX
      if (S->rc_variant_)
      {
        PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),zi,st::one(),pi,iflag),*iflag);
      }
#endif
      
      if(it>debugstate){
      std::cout << "Nach P-Sweep. " << *beta << std::endl;
      getchar();
      PHIST_CHK_IERR(SUBR(show_vector_elements)(r, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(show_vector_elements)(p, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(show_vector_elements)(q, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(show_vector_elements)(x, iflag),*iflag);
      }





      }else{
      correction_needed = false;  
      //PHIST_SOUT(PHIST_INFO,"DIESE ZEILE MUESSTE JEDE ITERATION ERSCHEINEN, WAS SIE AUCH TUT.\n");                                                            ///////////////////////////////////////////////////////////////////
      // apply operator, I-carp_sweep(...) to p. Note that our function carp_sweep operates
      // in place, so we first copy p to q again. The rhs vector is 0, which the
      // kernel lib should understand if we pass in NULL.

      //q=p-carp_sweep(A,sigma,B,bnul,p,omega);
      PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),p,st::zero(),q,iflag),*iflag);
#ifndef IS_COMPLEX
      if (S->rc_variant_)
      {
        PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),pi,st::zero(),qi,iflag),*iflag);
      }
#endif
      // double carp sweep in place, updates q to carp_sweep(p)
      PHIST_CHK_IERR(SUBR(carp_sweep)(A, sigma_r, sigma_i,bnul,q,qi,
          S->aux_,S->omega_,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),p,-st::one(),q,iflag),*iflag);
#ifndef IS_COMPLEX
      if (S->rc_variant_)
      {
        PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),pi,-st::one(),qi,iflag),*iflag);
      }
#endif
      
      if(it>debugstate){
      std::cout << "Nach Q-Sweep. " << std::endl;
      getchar();
      PHIST_CHK_IERR(SUBR(show_vector_elements)(r, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(show_vector_elements)(p, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(show_vector_elements)(q, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(show_vector_elements)(x, iflag),*iflag);
      }

      ////////////////////////////
      // update solution x      //
      ////////////////////////////
      
      //alpha = (r'*z)/(p'*q);
      ST denom  [nvec];
      MT denom_i[nvec];
      PHIST_CHK_IERR(SUBR(private_dotProd)(r,ri,z,zi,nvec,alpha,alpha_i,S->rc_variant_,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(private_dotProd)(p,pi,q,qi,nvec,denom,denom_i,S->rc_variant_,iflag),*iflag);
      MT minus_alpha_i[nvec];

      // stop updating x if the system is already converged (1-conv[j])
#ifdef IS_COMPLEX
        for (int j=0;j<nvec;j++)
        {
          alpha[j]=(ST)(1-conv[j])*(alpha[j]/denom[j]);
        }
        // update x <- x + alpha*p
        PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(alpha,p,st::one(),x,iflag),*iflag);
#else
      if (!S->rc_variant_)
      {
        for (int j=0;j<nvec;j++)
        {
          alpha[j]=(ST)(1-conv[j])*(alpha[j]/denom[j]);
        }
        // update x <- x + alpha*p
        PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(alpha,p,st::one(),x,iflag),*iflag);
      }
      else
      {
        for (int j=0;j<nvec;j++)
        {
          CT tmp1(alpha[j],alpha_i[j]);
          CT tmp2(denom[j],denom_i[j]);
          CT tmp3=tmp1/tmp2;
          tmp1*=(MT)(1-conv[j]);
          alpha[j]=ct::real(tmp3);
          alpha_i[j]    = ct::imag(tmp3);
          minus_alpha_i[j]=-ct::imag(tmp3);
        }
        // update x <- x + alpha*p
        PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(alpha,p,st::one(),x,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(minus_alpha_i,pi,st::one(),x,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(alpha,pi,st::one(),xi,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(alpha_i,p,st::one(),xi,iflag),*iflag);
      }
#endif
      
      if(it>debugstate){
      std::cout << "Nach X-Sweep. " << std::endl;
      getchar();
      PHIST_CHK_IERR(SUBR(show_vector_elements)(r, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(show_vector_elements)(p, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(show_vector_elements)(q, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(show_vector_elements)(x, iflag),*iflag);
    }
    
      

      if ( it%itcheck == 0)
      {
        for (int j=0;j<nvec;j++)
        {
          S->normR_old[j] = S->normR[j];
        }
        PHIST_CHK_IERR(SUBR(private_compResid)(A, nvec, sigma_r, sigma_i,
                         b, x, xi, NULL, NULL, S->normR,S->rc_variant_, iflag),*iflag);
        //std::cout << *(S->normR) <<" " << *(S->normR_old)<< " "<< (*(S->normR)/(*(S->normR_old))) << std::endl;
        


        // check for convergence. 
        // TODO - which convergence criterion should we use?
        // For the moment, we use ||r||_2/||b||_2 < tol
        numConverged=0;
        for (int j=0;j<nvec;j++)
        {
          if (S->normR[j]<reltol2[j])
          {
            conv[j]=1;
            numConverged++;
          }
        }

        
        if ( (it%itprint==0 && itprint>0) || (numConverged>=minConv))
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
          SUBR(private_printResid)(it, nvec, tmp, S->normR0_, S->normB_,S->conv);
          //std::cout <<  std::sqrt(*(S->normR)) << " " << std::sqrt((*(S->normR_old))) << " "<< std::sqrt(*(S->normR))/std::sqrt((*(S->normR_old))) <<std::endl;
          if(std::sqrt(*(S->normR))/std::sqrt(*(S->normR_old)) > cor_tol){
            correction_needed = true;
            std::cout <<  std::sqrt(*(S->normR)) << " " << std::sqrt((*(S->normR_old))) << " "<< std::sqrt(*(S->normR))/std::sqrt((*(S->normR_old))) <<std::endl;
            std::cout << "Cstep needed, because " << std::sqrt(*(S->normR))/std::sqrt((*(S->normR_old))) <<" > " << cor_tol <<std::endl;
          }
        }

        if (numConverged>=minConv)
        {
          // print footer
          SUBR(private_printResid)(it,-1,NULL,NULL,NULL,NULL);
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
      PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(minus_alpha,q,st::one(),r,iflag),*iflag);


#ifndef IS_COMPLEX
      if (S->rc_variant_)
      {
        PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(      alpha_i,qi,st::one(),r,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(minus_alpha,  qi,st::one(),ri,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(minus_alpha_i,q, st::one(),ri,iflag),*iflag);
      }
#endif
      
      if(it>debugstate){
      std::cout << "Nach R-Sweep. " << *alpha << std::endl;
      getchar();
      PHIST_CHK_IERR(SUBR(show_vector_elements)(r, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(show_vector_elements)(p, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(show_vector_elements)(q, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(show_vector_elements)(x, iflag),*iflag);
      }
      

//  z=apply_op(r,M):
// .. do nothing ...

      for (int j=0;j<nvec;j++)
      {
        //std::cout << it << " "<<r2_old[j]<< " "<< std::abs(r2_new[j])<< " "<< std::abs(r2_new[j]/r2_old[j]) << " "<<std::sqrt(std::abs(r2_old[j]))<< " "<<std::sqrt(std::abs(r2_new[j])) << " "<<std::sqrt(std::abs(r2_new[j]))/std::sqrt(std::abs(r2_old[j]))<< std::endl;;
        r2_old[j]=std::abs(r2_new[j]);
      }
      PHIST_CHK_IERR(SUBR(private_dotProd)(r,ri,z,zi,nvec,r2_new,NULL,S->rc_variant_,iflag),*iflag);
      //note: for a precond M!=I I think we need complex
      // arithmetic here because r!=z => above dotProd has imag!=0.
      // Only for symmetric preconditioning would we get a Hermitian
      // Lanczos matrix and thus a real beta on the diagonal.
      for (int j=0;j<nvec;j++)
      {
        beta[j]=std::abs(r2_new[j])/r2_old[j];
      }

      //p=z+beta*p;
#ifdef IS_COMPLEX
      ST tmp[nvec];
      for (int j=0;j<nvec;j++)
      {
        tmp[j]=(ST)beta[j];
      }
      PHIST_CHK_IERR(SUBR(mvec_vscale)(p,tmp,iflag),*iflag);
#else
      PHIST_CHK_IERR(SUBR(mvec_vscale)(p,beta,iflag),*iflag);
      if (S->rc_variant_)
      {
        PHIST_CHK_IERR(SUBR(mvec_vscale)(pi,beta,iflag),*iflag);
      }
#endif
      PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),z,st::one(),p,iflag),*iflag);
#ifndef IS_COMPLEX
      if (S->rc_variant_)
      {
        PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),zi,st::one(),pi,iflag),*iflag);
      }
#endif
      
      if(it>debugstate){
      std::cout << "Nach P-Sweep. " << *beta<< std::endl;
      getchar();
      PHIST_CHK_IERR(SUBR(show_vector_elements)(r, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(show_vector_elements)(p, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(show_vector_elements)(q, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(show_vector_elements)(x, iflag),*iflag);
      }
      
      //std::cout << std::sqrt(*(S->normR)) << std::endl;

      // the code below injects some faults in the computation.
      // (see Praktikumsbericht by Florian Fritzen @DLR 2014)
#if 0
    if(it % 70 == 0){
      //std::cout << "Destruction completed " <<std::endl;
      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(r, 0.0, 2000, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(ri, 0.0, 2000, iflag),*iflag);

      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(p, 0.0, 2000,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(pi, 0.0, 2000,iflag),*iflag);
      
      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(q, 0.0, 2000, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(qi, 0.0, 2000, iflag),*iflag);

      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(z, 0.0, 2000, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(zi, 0.0, 2000, iflag),*iflag);
      
      //TODO
    }
    
    if(true){
      if(it==3 || it==9|| it==27|| it==31|| it==56|| it==64 || it==81 || it==100){
        PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(x, 0.0, 2000, iflag),*iflag);
        PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(xi, 0.0, 2000, iflag),*iflag);
      }
    }
    
    
    //TEST CASES IF PROCESS IS KILLED BY WHATEVER
    int mynode;
    MPI_Comm_rank(MPI_COMM_WORLD,&mynode);
    //std::cout << "Rank: " << mynode << std::endl;
    if(mynode == 2 && it==77){
      //2.001
      //8.001
      //17.001
      std::cout << std::endl <<" DONE " << std::endl << std::endl;
      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(r, 0.0, 17000, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(ri, 0.0, 17000, iflag),*iflag);

      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(p, 0.0, 17000,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(pi, 0.0, 17000,iflag),*iflag);
      
      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(q, 0.0, 17000, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(qi, 0.0, 17000, iflag),*iflag);

      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(z, 0.0, 17000, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(zi, 0.0, 17000, iflag),*iflag);

      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(x, 0.0, 17000, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(xi, 0.0, 17000, iflag),*iflag);
    }
    
    if(mynode == 3 && it==277){
      std::cout << std::endl <<" DONE " << std::endl << std::endl;
      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(r, 0.0, 17000, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(ri, 0.0, 17000, iflag),*iflag);

      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(p, 0.0, 17000,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(pi, 0.0, 17000,iflag),*iflag);
      
      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(q, 0.0, 17000, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(qi, 0.0, 17000, iflag),*iflag);

      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(z, 0.0, 17000, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(zi, 0.0, 17000, iflag),*iflag);

      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(x, 0.0, 17000, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(xi, 0.0, 17000, iflag),*iflag);
    }

    if(mynode == 1 && it==15){
      std::cout << std::endl <<" DONE " << std::endl << std::endl;
      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(r, 0.0, 17000, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(ri, 0.0, 17000, iflag),*iflag);

      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(p, 0.0, 17000,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(pi, 0.0, 17000,iflag),*iflag);
      
      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(q, 0.0, 17000, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(qi, 0.0, 17000, iflag),*iflag);

      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(z, 0.0, 17000, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(zi, 0.0, 17000, iflag),*iflag);

      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(x, 0.0, 17000, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(destroy_vector_elements_at_random)(xi, 0.0, 17000, iflag),*iflag);
    }
#endif //fault injection

    }// endof else
    }// CG iterations (it)
    //...
    
    // take square root so that normR is indeed the norm, not the
    // squared 2-norm:
    for (int j=0;j<nvec;j++)
    {
      S->normR[j]=std::sqrt(S->normR[j]);
    }
  
    // to save memory, deallocate CG vector data
    PHIST_CHK_IERR(SUBR(private_carp_cgState_dealloc)(S,iflag),*iflag);

//  std::cout << "CSteps: " << cor_count << std::endl;
  *iflag=nvec-numSolved;
  return;
}

// compute residual r=b-(A-sI)x and ||r||_2^2 in nrms2. If r and ri are NULL, a temporary
// vector is used and discarded.
void SUBR(private_compResid)(TYPE(const_sparseMat_ptr) A, int nvec, _MT_ const sigma_r[], _MT_ const sigma_i[],
                       TYPE(const_Mvec_ptr) Vproj,
                       TYPE(const_mvec_ptr) Rhs,
                       TYPE(const_mvec_ptr) x, TYPE(const_mvec_ptr) xi, TYPE(sdMat_ptr) xp, TYPE(sdMat_ptr) xpi,
                       TYPE(mvec_ptr) r,        TYPE(mvec_ptr) ri, TYPE(const_sdMat_ptr) rp, TYPE(const_sdMat_ptr) rpi,
                       _MT_  *nrms2, bool rc_variant, int *iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);

  if ((xp!=NULL)||(rp!=NULL))
  {
    *iflag=-99;
    return;
  }
#ifndef IS_COMPLEX
  if ((xpi!=NULL)||(rpi!=NULL)&&rc_variant)
  {
    *iflag=-99;
    return;
  }
#endif
  ST shifts[nvec];
  for (int i=0;i<nvec;i++)
  {
    shifts[i]=-(ST)sigma_r[i];
#ifdef IS_COMPLEX
    shifts[i]-=sigma_i[i]*st::cmplx_I();
#endif
  }

  TYPE(mvec_ptr) R=r, RI=ri;
  if (R==NULL)
  {
    const_map_ptr_t map;
    PHIST_CHK_IERR(SUBR(mvec_get_map)(x,&map,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_create)(&R,map,nvec,iflag),*iflag);
  }
  if (RI==NULL && rc_variant)
  {
    const_map_ptr_t map;
    PHIST_CHK_IERR(SUBR(mvec_get_map)(xi,&map,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_create)(&RI,map,nvec,iflag),*iflag);
  }

  // r = b-(A-sI)x
  // r = b-(A-sr)xr - si*xi
  //    -i[(A-sr)xi - si*xr]
  
  // r=b
  PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),Rhs,st::zero(),R,iflag),*iflag);
  // r=b-(A-sI)*x
  PHIST_CHK_IERR(SUBR(sparseMat_times_mvec_vadd_mvec)
    (-st::one(),A,shifts,x,st::one(),R,iflag),*iflag);
  if (rc_variant)
  {
    //ri=-(A-srI)xi
    PHIST_CHK_IERR(SUBR(sparseMat_times_mvec_vadd_mvec)
      (-st::one(),A,shifts,xi,st::zero(),RI,iflag),*iflag);
    //r=r-si*xi
    for (int i=0;i<nvec;i++) shifts[i]=-sigma_i[i];
    PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(shifts, xi,st::one(),R,iflag),*iflag);
    // ri+=si*xr
    for (int i=0;i<nvec;i++) shifts[i]=sigma_i[i];
    PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(shifts,x,st::one(),RI,iflag),*iflag);
  }
  // now compute the 2-norm of each column
#ifdef IS_COMPLEX
          ST tmp[nvec];
#else
          ST* tmp=nrms2;
#endif          
  PHIST_CHK_IERR(SUBR(private_dotProd)(R, RI, R, RI, nvec, tmp, NULL, rc_variant, iflag),*iflag);
#ifdef IS_COMPLEX
  for (int j=0;j<nvec;j++)
  {
    nrms2[j]=st::real(tmp[j]);
  }
#endif
  if (R!=r && R!=NULL)
  {
    PHIST_CHK_IERR(SUBR(mvec_delete)(R,iflag),*iflag);
  }
  if (RI!=ri && RI!=NULL)
  {
    PHIST_CHK_IERR(SUBR(mvec_delete)(RI,iflag),*iflag);
  }
  return;
}

// compute v'*w, where v and w may have an imaginary part even in real (S/D) case.
// in the complex case (C/Z), vi and wi are ignored, the result dots is complex and
// dotsi is ignored. In the real case, if vi and wi are NULL, they are assumed
// to be 0 and dotsi is ignored.
void SUBR(private_dotProd)(TYPE(const_mvec_ptr) v, TYPE(const_mvec_ptr) vi, TYPE(const_sdMat_ptr) vp, TYPE(const_sdMat_ptr) vpi,
                           TYPE(const_mvec_ptr) w, TYPE(const_mvec_ptr) wi, TYPE(const_sdMat_ptr) wp, TYPE(const_sdMat_ptr) wpi,
                           int nvec, _ST_ *dots, _MT_ *dotsi, bool rc_variant, int *iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);

  *iflag=0;
  PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(v,w,dots,iflag),*iflag);
  if (vp!=NULL || wp!=NULL)
  {
    *iflag=-99;
    return;
  }
#ifndef IS_COMPLEX
  if (rc_variant)
  {
    ST tmp[nvec];
    PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(vi,wi,tmp,iflag),*iflag);
    for (int j=0;j<nvec;j++)
    {
      dots[j]+=tmp[j];
    }
    if ( (v!=w || vi!=wi) && dotsi!=NULL)
    {
      PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(v,wi,dotsi,iflag),*iflag);      
      PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(vi,w,tmp,iflag),*iflag);      
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
    if (vpi!=NULL || wpi!=NULL)
    {
      *iflag=-99;
      return;
    }
  }
#endif
  return;
}

// print residual info. If we get nvec=0, the input arrays are ignored
// and a header is printed.
void SUBR(private_printResid)(int it, int nvec, _ST_ const* normR, 
        _MT_ const* normR0, _MT_ const* normB, int const* locked)
{
#include "phist_std_typedefs.hpp"
  const char* carp_label = "CARP_CG";
  if (nvec==0)
  {
    PHIST_SOUT(PHIST_INFO,"%s\tit\t||r||\t||r||/||b||\t||r||/||r0||\n",carp_label);
    //PHIST_SOUT(PHIST_INFO,"TEST1.\n"); 
  }
  else if (nvec==-1)
  {
    PHIST_SOUT(PHIST_INFO,"%s\tfinished after %d iterations.\n",carp_label,it);
   // PHIST_SOUT(PHIST_INFO,"TEST2.\n"); 
  }
  else if (PHIST_OUTLEV<=PHIST_INFO && nvec>2)
  {
    // PHIST_SOUT(PHIST_INFO,"TEST3.\n"); 
    // print maximum and minimum residual norms only
    int max_pos=0, min_pos=0;
    ST nrmR[2];
    MT nrmB[2],nrmR0[2];
    for (int i=1;i<nvec;i++)
    {
      if (st::real(normR[i])>st::real(normR[max_pos])) max_pos=i;
      if (st::real(normR[i])<st::real(normR[min_pos])) min_pos=i;
    }
    nrmR[0]=normR[min_pos]; nrmR[1]=normR[max_pos];
    nrmR0[0]=normR0[min_pos]; nrmR0[1]=normR0[max_pos];
    nrmB[0]=normB[min_pos]; nrmB[1]=normB[max_pos];

    PHIST_SOUT(PHIST_INFO,"min and max residuals:\n");
    SUBR(private_printResid)(it,2,nrmR,nrmR0,nrmB,(int const*)NULL);
  }
  else
  {
    MT tmp=mt::sqrt(st::real(normR[0]));
    std::string lock_str="";
    //PHIST_SOUT(PHIST_INFO,"TEST4 %e.\n",tmp/normB[0]); 
    if (locked!=NULL) lock_str=locked[0]?"(converged)":"";
    PHIST_SOUT(PHIST_INFO,"%s %d\t%e\t%e\t%e\t%s\n",carp_label,it,
          tmp,tmp/normB[0],tmp/normR0[0],lock_str.c_str());
    for (int j=1;j<nvec;j++)
    {
      if (locked!=NULL) lock_str=locked[j]?"(converged)":"";
      MT tmp=mt::sqrt(st::real(normR[j]));
      PHIST_SOUT(PHIST_INFO,"%s\t\t%e\t%e\t%e\t%s\n",carp_label,
             tmp,tmp/normB[j],tmp/normR0[j],lock_str.c_str());
    }
  }
}

void SUBR(destroy_vector_elements_at_random)(TYPE(mvec_ptr) V, double probability, int nloc, int *iflag) 
{
#include "phist_std_typedefs.hpp"
//std::cout << "Function used." <<std::endl;
  ST* V_val;
  lidx_t lda;
  int nvec;
  //ST tmp;
  double rand_prob;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nvec,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_extract_view)(V,&V_val,&lda,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_from_device)(V,iflag),*iflag);
  //std::cout << "nloc " << nloc <<std::endl;
  for (int i=0; i<nloc; i++)
  {
    //std::cout << "First loop." <<std::endl;
    for (int j=0;j<nvec; j++)
    {
      rand_prob=std::abs(mt::rand());
      //std::cout << "rand_prob "<< rand_prob <<std::endl;
      //std::cout << "Second loop." <<std::endl
      if (rand_prob>probability){
        //std::cout << "Element destroyed." <<std::endl;
       
        //tmp = V_val[i*lda+j];
        V_val[i*lda+j] =st::rand();
        //V_val[i*lda+j] =0;
        //std::cout << j << std::endl;
        //std::cout << "Changed Value from " << tmp << " to " << V_val[i*lda+j] << std::endl;
      }
    }
  }
  PHIST_CHK_IERR(SUBR(mvec_to_device)(V,iflag),*iflag);
  //std::cout <<" Destruction complete " << std::endl;
}

void SUBR(show_vector_elements)(TYPE(mvec_ptr) V, int *iflag) 
{
  #include "phist_std_typedefs.hpp"
ST* V_val;
lidx_t lda;
int nvec;
int nloc = 1;
PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nvec,iflag),*iflag);
PHIST_CHK_IERR(SUBR(mvec_extract_view)(V,&V_val,&lda,iflag),*iflag);
PHIST_CHK_IERR(SUBR(mvec_from_device)(V,iflag),*iflag);
for (int i=0; i<nloc; i++)
  {
    for (int j=0;j<nvec; j++)
    {
#ifdef PHIST_MVECS_ROW_MAJOR
        std::cout << "Element is "<< -V_val[i*lda+j] << std::endl;
#else
        std::cout << "Element is "<< -V_val[j*lda+i] << std::endl;
#endif
    }
  }
}
