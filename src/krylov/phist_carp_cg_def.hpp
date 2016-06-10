// in this implementation, there is no mixed real/complex arithmetic
// simply since we don't have a kernel lib interface and we also want
// to support kernel libs like epetra or our own fortran variant, which
// do not have complex arithmetic let alone mixed functionality.

//TODO - our pseudo complex dot product and MVM are probably very inefficient 
//      compared to using proper complex vectors, we should check if the kernel 
//      lib supports mixed real/complex arithmetic and use it if possible.

// We incorporate the case of a shifted and/or 'augmented' system
//
// |A-sigma[j]I Vproj ||x |   |b |
// | Vproj'       0   ||x'| = |b'|
//
// The small additional components x', b' etc. are denoted by a 'p' in the code, e.g. xp, bp,
// and represented as sdMats.

///////////////////////////////////////////////////////////////////////////////////////////////////
// private helper functions                                                                      //
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "phist_carp_cg_kernels_decl.hpp"


//! compute residual, note that this function only accepts a real right-hand side vector in the RC case
void SUBR(my_compResid)(TYPE(x_sparseMat) const* A, 
                       TYPE(const_mvec_ptr) Rhs,
                       TYPE(x_mvec) const* x, 
                       TYPE(x_mvec)* r,
                       _MT_  *nrms, int *iflag);

// pretty-print convergence history. nvec=0: print header. nvec=-1: print footer
void SUBR(my_printResid)(int it, int nvec, _ST_ const* normR, 
        _MT_ const* normR0, _MT_ const* normB, int const* locked);

// allocate CG vector blocks
void SUBR(my_carp_cgState_alloc)(TYPE(carp_cgState_ptr) S, int* iflag);
// allocate CG vector blocks
void SUBR(my_carp_cgState_dealloc)(TYPE(carp_cgState_ptr) S, int* iflag);

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
  
  phist_const_map_ptr map=NULL;
  PHIST_CHK_IERR(SUBR(sparseMat_get_row_map)(A,&map,iflag),*iflag);
  phist_lidx nloc;
  PHIST_CHK_IERR(phist_map_get_local_length(map,&nloc,iflag),*iflag);
  
  *state = new TYPE(carp_cgState);
  (*state)->iflag=-2;// not initialized (need to call reset())
  
  // note: this doesn't allocate memory yet. The iterate() function below
  // alocates and deallocates the vectors each time it is called because
  // the overhead may be big if many systems need to be solved as in FEAST.
  // A pipelining strategy for the linear systems with small memory consumption
  // can then be used to solve for ~1000 RHS using blocks of ~8.
  (*state)->p_=new TYPE(x_mvec);
  (*state)->q_=new TYPE(x_mvec);
  (*state)->r_=new TYPE(x_mvec);

  (*state)->A_=new TYPE(x_sparseMat);
  (*state)->A_->A_=A;
  (*state)->A_->Vproj_=Vproj;
  (*state)->A_->sigma_r_ = new MT[nvec];
  (*state)->A_->sigma_i_ = new MT[nvec];

  (*state)->conv = new int[nvec];
  
  (*state)->normR0_= new MT[nvec];
  (*state)->normB_= new MT[nvec];
  (*state)->normR= new MT[nvec];
  (*state)->normR_old= new MT[nvec];

  (*state)->beta_ =  new MT[nvec];
  (*state)->alpha_ = new ST[nvec];
  (*state)->alpha_i_ = new MT[nvec];
  (*state)->omega_=new MT[nvec]; // relaxation parameter for col j
                                 // moment just set it to 1, which
                                 // gave good results for Graphene
                                 // in the matlab tests.

  (*state)->rc_variant_=false;
  for (int i=0; i<nvec; i++)
  {
    (*state)->omega_[i]=mt::one();
    (*state)->A_->sigma_r_[i]=sigma_r[i];
    (*state)->A_->sigma_i_[i]=sigma_i[i];
#ifndef IS_COMPLEX
    if (sigma_i[i]!=mt::zero()) (*state)->rc_variant_=true;
#endif
  }
  
  if (!(*state)->rc_variant_)
  {
    delete [] (*state)->A_->sigma_i_;
    (*state)->A_->sigma_i_=NULL;
  }
  
  (*state)->nvec_=nvec;
  (*state)->nproj_=0;
  if(Vproj!=NULL)
  {
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(Vproj,&(*state)->nproj_,iflag),*iflag);
  }

  // setup the CARP kernel and get the required data structures:
  void* aux=NULL;
  if ((*state)->rc_variant_)
  {
#ifdef IS_COMPLEX
  *iflag=-1;
  PHIST_SOUT(PHIST_ERROR,"in complex arithmetic, the flag rc_variant_ should not be set to true\n"
                        "(file %s, line %d)\n",__FILE__,__LINE__);
  return;
#else
//    *iflag=PHIST_NOT_IMPLEMENTED;
//    PHIST_SOUT(PHIST_ERROR,"The real variant of CARP-CG with complex shifts is currently broken and we therefore abort here\n");
//    return;
    PHIST_CHK_IERR(SUBR(carp_setup_rc)(A,nvec,sigma_r,sigma_i,
      &aux, iflag),*iflag);
#endif
  }
  else
  {
    _ST_ *sigma;
#ifdef IS_COMPLEX
    sigma=new _ST_[nvec];
    for (int i=0; i<nvec; i++) sigma[i]=sigma_r[i]+sigma_i[i]*st::cmplx_I();
#else
    sigma=sigma_r;
#endif
    PHIST_CHK_IERR(SUBR(carp_setup)(A,nvec,sigma,
      &aux, iflag),*iflag);
#ifdef IS_COMPLEX
    delete [] sigma;
#endif
  }
  (*state)->aux_=aux;

  for (int j=0;j<nvec;j++)
  {
    (*state)->conv[j]=0;
    (*state)->normR0_[j]=-mt::one(); // not initialized
    (*state)->normB_[j]=-mt::one(); // not initialized
    (*state)->normR[j]=-mt::one(); // not initialized
  }
}

/////////////////////////////////////////////////////////////////////////////////////
//=================================================================================//
/////////////////////////////////////////////////////////////////////////////////////

void SUBR(my_carp_cgState_alloc)(TYPE(carp_cgState_ptr) S, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  
  int nvec=S->nvec_;
  int nproj=S->nproj_;
  const void* map=NULL;
  PHIST_CHK_IERR(SUBR(mvec_get_map)(S->b_,&map,iflag),*iflag);
  
  bool rc = (S->rc_variant_!=0);
    
  PHIST_CHK_IERR(S->q_->allocate(map,nvec,nproj,rc,iflag),*iflag);
  PHIST_CHK_IERR(S->r_->allocate(map,nvec,nproj,rc,iflag),*iflag);
  PHIST_CHK_IERR(S->p_->allocate(map,nvec,nproj,rc,iflag),*iflag);
}

void SUBR(my_carp_cgState_dealloc)(TYPE(carp_cgState_ptr) S, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);

  S->p_->deallocate();
  S->q_->deallocate();
  S->r_->deallocate();
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

    delete [] state->omega_;

    delete [] state->beta_;
    delete [] state->alpha_;
    delete [] state->alpha_i_;

    delete state->p_;
    delete state->q_;
    delete state->r_;

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
  }

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

////////////////////////////////////////////////////////////////////////////////
// settings for the algorithm                                                 //
////////////////////////////////////////////////////////////////////////////////


  // some internal settings 
  bool correction_step = false, correction_needed=false; 
  int cor_count = 0;
  double cor_tol = 0.99;

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

  // used for premature termination of the loop if requested
  int minConv=abortIfOneConverges? 1: S->nvec_;  

  int numConverged=0;
  
////////////////////////////////////////////////////////////////////////////////
// set some useful pointers for convenience                                   //
////////////////////////////////////////////////////////////////////////////////
  
  TYPE(mvec_ptr) bnul=NULL; // we can just pass in b=NULL if all entries for a carp_sweep
                            // are 0 (during CG iteration), for clarity we give it a name

  // get some pointers to ease the notation
  TYPE(x_sparseMat) *A=S->A_;
  int nvec=S->nvec_;
  TYPE(const_mvec_ptr) b=S->b_;
  // copy or view input vectors in our internal x_mvec type
  TYPE(x_mvec)* x=NULL;
  if (S->rc_variant_) // rc variant
  {
    PHIST_CHK_IERR(x=new TYPE(x_mvec)(X_r,X_i,S->nproj_,iflag),*iflag);
  }
  else
  {
    PHIST_CHK_IERR(x=new TYPE(x_mvec)(X_r,NULL,S->nproj_,iflag),*iflag);
  }

  TYPE(x_mvec)* r =S->r_;
  TYPE(x_mvec)* q =S->q_;
  TYPE(x_mvec)* p =S->p_;
    
  // CG (Lanczos) coefficients
  ST* alpha=S->alpha_;
  MT* alpha_i=S->alpha_i_;
  MT* beta=S->beta_;
    
  int *conv=S->conv;
    
////////////////////////////////////////////////////////////////////////////////
// check if input vectors have correct #cols                                  //
////////////////////////////////////////////////////////////////////////////////
   
  int nvecX,nvecB;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(b,&nvecB,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(x->v_,&nvecX,iflag),*iflag);
      
  if (nvec!=nvecB || nvec!=nvecX)
  {
    PHIST_SOUT(PHIST_ERROR,"input vectors to %s must have same num vectors as block size \n"
                           "passed to constructor (expected %d, found nvec(X)=%d, nvec(B)=%d instead)\n"
                           "(in %s, line %d)\n",
                           __FUNCTION__,nvec,nvecX,nvecB,__FILE__,__LINE__);
    *iflag=PHIST_INVALID_INPUT;
    return;      
  }

#ifndef IS_COMPLEX
  if (x->vi_!=NULL)
  {
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(x->vi_,&nvecX,iflag),*iflag);
    if (nvec!=nvecX)
    {
      PHIST_SOUT(PHIST_ERROR,"input vectors X_i to %s must have same num vectors as X_r and B\n",__FUNCTION__);
      *iflag=PHIST_INVALID_INPUT;
      return;
    }
  }
#endif    

////////////////////////////////////////////////////////////////////////////////
// initialization                                                             //
////////////////////////////////////////////////////////////////////////////////


  // allocate CG vectors: one per rhs. We allocate new memory each time
  // the iterate function is called and release it afterwards to reduce
  // the memory consumption between calls. This allows e.g. having an  
  // array of carp-cg objects that are used one by one for different shifts,
  // vector blocks or intervals in the spectrum.
  PHIST_CHK_IERR(SUBR(my_carp_cgState_alloc)(S,iflag),*iflag);

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
    PHIST_CHK_IERR(SUBR(my_compResid)(A, b, x, NULL, S->normR, iflag),*iflag);
    for (int j=0;j<nvec;j++)
    {
      S->normR0_[j]=std::sqrt(S->normR[j]);
      S->normR_old[j] = S->normR0_[j];
    }
  }

  S->iflag=1; // unconverged

  MT reltol2[nvec];

  for (int j=0;j<nvec;j++)
  {
    reltol2[j]=tol*tol*S->normR[j];
    reltol2[j]=std::max(reltol2[j],tol*tol*S->normB_[j]*S->normB_[j]);
  }

  // initial Kaczmarz/CARP sweep. Note that our function carp_sweep operates
  // in place, but we do not want to update x right now, we just want to
  // get a CG direction.

  //ALG r=b-OP*x
  PHIST_CHK_IERR(SUBR(x_mvec_add_mvec)(st::one(),x,st::zero(),r,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(x_carp_sweep)(A,b,r,S->aux_,S->omega_,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(x_mvec_add_mvec)(-st::one(),x,st::one(),r,iflag),*iflag);
  //ALG p=r
  PHIST_CHK_IERR(SUBR(x_mvec_add_mvec)(st::one(),r,st::zero(),p,iflag),*iflag);
  //r2_new = ||r||_2^2. For technical reasons we store it as complex type if IS_COMPLEX
  ST r2_new[nvec];
  MT r2_old[nvec];
  PHIST_CHK_IERR(SUBR(x_mvec_dot_mvec)(r,r,r2_new,NULL,iflag),*iflag);

  if (itprint>0)
  {
    // prints a header for the convergence history
    SUBR(my_printResid)(0,0,NULL,NULL,NULL,NULL);
    // print initial residual
    SUBR(my_printResid)(0,nvec,r2_new,S->normR0_,S->normB_,S->conv);
  }

  //ALG for (it=1..maxIter
  for (int it=1;it<maxIter; it++)
  {
    //ALG correction_needed=false
    correction_needed=false;

    //ALG q = OP*p

    //q=p-carp_sweep(A,sigma,B,bnul,p,omega);
    PHIST_CHK_IERR(SUBR(x_mvec_add_mvec)(st::one(),p,st::zero(),q,iflag),*iflag);
    // double carp sweep in place, updates q to carp_sweep(p)
    PHIST_CHK_IERR(SUBR(x_carp_sweep)(A,bnul,q, S->aux_,S->omega_,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(x_mvec_add_mvec)(st::one(),p,-st::one(),q,iflag),*iflag);

    // explicit computation of residual (helps to overcome soft errors, see 
    // `self-stabilizing CG' described in Sao & Vuduc, ScalA'14 proceedings. 
    //ALG if (correction_step)
    if (correction_step)
    {
      PHIST_SOUT(PHIST_INFO,"CARP_CG - correction step\n");
      cor_count++;

      //ALG r=b-OP*x
      
      //r=(dkswp(A,0,x)-I)x
      PHIST_CHK_IERR(SUBR(x_mvec_add_mvec)(st::one(),x,st::zero(),r,iflag),*iflag);

      // double carp sweep in place, updates r=dkswp(A-sI,omega,r)
      PHIST_CHK_IERR(SUBR(x_carp_sweep)(A, b,r,S->aux_,S->omega_,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(x_mvec_add_mvec)(-st::one(),x,st::one(),r,iflag),*iflag);
    }
    //ALG end if
    
    ////////////////////////////
    // update solution x      //
    ////////////////////////////
    
    //ALG alpha = (r'*r)/(p'*q);
    ST denom  [nvec];
    MT denom_i[nvec];
    PHIST_CHK_IERR(SUBR(x_mvec_dot_mvec)(r,p,alpha,alpha_i,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(x_mvec_dot_mvec)(p,q,denom,denom_i,iflag),*iflag);

    // stop updating x if the system is already converged (1-conv[j])
    // alpha = (r'p)/(p'q), may be complex
    for (int j=0;j<nvec;j++)
    {
      if (!S->rc_variant_)
      {
        alpha[j]=(ST)(1-conv[j])*(alpha[j]/denom[j]);
      }
#ifndef IS_COMPLEX
      else
      {
        CT tmp1(alpha[j],alpha_i[j]);
        tmp1*=(MT)(1-conv[j]);
        CT tmp2(denom[j],denom_i[j]);
        CT tmp3=tmp1/tmp2;
        alpha[j]=ct::real(tmp3);
        alpha_i[j]    = ct::imag(tmp3);
      }
#endif
    }
    
    //ALG x = x + alpha*p
    PHIST_CHK_IERR(SUBR(x_mvec_vadd_mvec)(alpha,alpha_i,p,st::one(),x,iflag),*iflag);
    
    // once in a while we check the actual residual of the linear system, b-Ax
    // and decide if we need a correction step and if we are converged
    // TODO - the check for correction steps is quite rare (every 10 iterations),
    //        but I suspect that Florian (who implemented it) did it this way because
    //        checking with the norm of r during CGMN was not working.
    //ALG if (mod(it,itcheck)==0)
    if ( it%itcheck == 0)
    {
      for (int j=0;j<nvec;j++)
      {
        S->normR_old[j] = (S->normR[j]);
      }
      //ALG ||r||_2 = ||b-Ax||_2
      PHIST_CHK_IERR(SUBR(my_compResid)(A,b,x,NULL,S->normR,iflag),*iflag);
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
        SUBR(my_printResid)(it, nvec, tmp, S->normR0_, S->normB_,S->conv);  
      }
      
      // check if sufficient systems have converged and exit the loop if so
      //ALG if ||r||_2/||r_0||_2<tol
      if (numConverged>=minConv)
      {
        // print footer
        SUBR(my_printResid)(it,-1,NULL,NULL,NULL,NULL);
        numSolved++; 
        //ALG break
        break;
      }
      //ALG end if

      // check if the next step should be a correction step
      //ALG if ||r||_2/||r_old||_2 > 0.99
      for (int i=0;i<nvec;i++)
      {
        //ALG correction_needed=true
        if (std::sqrt(*(S->normR))/std::sqrt(*(S->normR_old)) > cor_tol)
        {
          std::cout <<  std::sqrt(*(S->normR)) << " " << std::sqrt((*(S->normR_old))) << " "<< std::sqrt(*(S->normR))/std::sqrt((*(S->normR_old))) <<std::endl;
          std::cout << "Cstep(!) needed, because " << std::sqrt(*(S->normR))/std::sqrt((*(S->normR_old))) <<" > " << cor_tol <<std::endl;
                                              
          correction_needed = true;
        }
      }
      //ALG end if
    }
    //ALG end if
    
    // regular CG step?
    if (!correction_step)
    {
      //ALG r=r-alpha*q;
      ST min_alpha[nvec];
      MT min_alpha_i[nvec];
      for (int j=0;j<nvec;j++)
      {
        min_alpha[j]=-alpha[j];
        min_alpha_i[j]=-alpha_i[j];
      }
      PHIST_CHK_IERR(SUBR(x_mvec_vadd_mvec)(min_alpha,min_alpha_i,
                q,st::one(),r,iflag),*iflag);
    }
      
    //ALG if (correction_step)
    if (correction_step)
    {
      //ALG beta= r'q/p'q
      ST rq  [nvec];
      MT rq_i[nvec];
      ST pq  [nvec];
      MT pq_i[nvec];

      PHIST_CHK_IERR(SUBR(x_mvec_dot_mvec)(r,q,rq,rq_i,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(x_mvec_dot_mvec)(p,q,pq,pq_i,iflag),*iflag);

      for (int j=0;j<nvec;j++)
      {
#ifndef IS_COMPLEX
        CT tmp1(rq[j],rq_i[j]);
        CT tmp2(pq[j],pq_i[j]);
        CT tmp3=-tmp1/tmp2;
        beta[j]=ct::real(tmp3);
#else
        beta[j]=ct::real(rq[j]/pq[j]);
#endif       
      }
    }//ALG else
    else
    {
      //ALG beta = r'r/r_old'r_old
      for (int j=0;j<nvec;j++)
      {
        r2_old[j]=std::abs(r2_new[j]);
      }
      PHIST_CHK_IERR(SUBR(x_mvec_dot_mvec)(r,r,r2_new,NULL,iflag),*iflag);
      for (int j=0;j<nvec;j++)
      {
        beta[j]=std::abs(r2_new[j])/r2_old[j];
      }
    }//ALG end if

    //ALG p=r+beta*p
    //TODO - add kernel function to do this, mvec_vadd_mvec only allows multiple scalars for the added vector
    PHIST_CHK_IERR(SUBR(x_mvec_vscale)(p,beta,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(x_mvec_add_mvec)(st::one(),r,st::one(),p,iflag),*iflag);

    //ALG correction_step=correction_needed
    correction_step=correction_needed;
    correction_needed=false;
  }
  //ALG end for

    
  // take square root so that normR is indeed the norm, not the
  // squared 2-norm:
  for (int j=0;j<nvec;j++)
  {
    S->normR[j]=std::sqrt(S->normR[j]);
  }
  
  // to save memory, deallocate CG vector data
  PHIST_CHK_IERR(SUBR(my_carp_cgState_dealloc)(S,iflag),*iflag);

  // copy the solution from the x_mvec
  if (x->v_!=X_r) PHIST_CHK_IERR(SUBR(mvec_set_block)(X_r,x->v_,0,nvec-1,iflag),*iflag);
  if (x->vi_!=X_i&&S->rc_variant_) PHIST_CHK_IERR(SUBR(mvec_set_block)(X_i,x->vi_,0,nvec-1,iflag),*iflag);
  
  // this is just a wrapper object, the arguments X_r and X_i carry the solution
  delete x;
  *iflag=nvec-numSolved;
  return;
}

// compute residual r=b-(A-sI)x and ||r||_2^2 in nrms2. If r and ri are NULL, a temporary
// vector is used and discarded.
void SUBR(my_compResid)(TYPE(x_sparseMat) const* A,
                       TYPE(const_mvec_ptr) Rhs,
                       TYPE(x_mvec) const* x,
                       TYPE(x_mvec)* r,
                       _MT_  *nrms2, int *iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);

  TYPE(x_mvec) *R=NULL;
  
  bool  rc= rc_variant(A,x,x);
  bool aug=aug_variant(A,x,x);
  
  int nvec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(Rhs,&nvec,iflag),*iflag);

  int naug=0;
  if (aug)
  {
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(A->Vproj_,&naug,iflag),*iflag);
  }
  
  if (r) 
  {
    R=r;
  }
  else
  {
    phist_const_map_ptr map;
    PHIST_CHK_IERR(SUBR(mvec_get_map)(x->v_,&map,iflag),*iflag);
    R=new TYPE(x_mvec);
  PHIST_CHK_IERR(R->allocate(map,nvec,naug,rc,iflag),*iflag);
  }
  
  // r=-(A-sigma[j]I)x
  PHIST_CHK_IERR(SUBR(x_sparseMat_times_mvec)
    (-st::one(),A,x,st::zero(),R,iflag),*iflag);

  // r+=b
  PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),Rhs,st::one(),R->v_,iflag),*iflag);

  // now compute the 2-norm of each column
#ifdef IS_COMPLEX
          ST tmp[nvec];
#else
          ST* tmp=nrms2;
#endif          
  PHIST_CHK_IERR(SUBR(x_mvec_dot_mvec)(R, R, tmp, NULL,iflag),*iflag);
#ifdef IS_COMPLEX
  for (int j=0;j<nvec;j++)
  {
    nrms2[j]=st::real(tmp[j]);
  }
#endif
  if (R!=r && R!=NULL)
  {
    delete R;
  }
  return;
}

// print residual info. If we get nvec=0, the input arrays are ignored
// and a header is printed.
void SUBR(my_printResid)(int it, int nvec, _ST_ const* normR, 
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
    SUBR(my_printResid)(it,2,nrmR,nrmR0,nrmB,(int const*)NULL);
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

// kernel implementation for this data type
#include "phist_carp_cg_kernels_def.hpp"
