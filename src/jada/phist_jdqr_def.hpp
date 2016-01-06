//! tries to compute a given number of eigenpairs (num_eig) of 
//! an general operator A using the Jacobi-Davidson QR (JDQR)
//! method. We allow for generalized EVPs here, A*x=lambda*B*x,
//! where B should be positive definite.
//!
//! input arguments:
//!
//! A: pointer to the operator(-wrapper) A
//! B: pointer to the "mass-matrix" operator (if B==NULL, B=I is assumed)
//!
//! X:  Should be allocated with           
//!     at least <num_eigs> columns (num_eigs+1 in the real case         
//!     to avoid the situation where the last eigenvalue to con-         
//!     verge is a complex pair).                                        
//!                                                                      

//! is_cmplx,evals and resid should be pre-allocated arrays of size at least opts.numEigs
//! (in the real case opts.numEigs+1)
//! 
//! 
//! output arguments:
//!
//! num_eigs: number of converged eigenpairs
//! num_iters: number of iterations performed
//! X: latest eigenvector approximations
//! evals: latest eigenvalue approximations
//! resid: Ritz residuals
//! is_cmplx: in the real case, if is_cmplx[i]==is_cmplx[i+1]=1, then
//! a evals[i] +/- evals[i+1]*I forms a complex conjugate pair,
//! and the corresponding eigenvectors are X(:,i)+/-X(:,i+1)*I.
//! In the complex case is_cmplx is not referenced.
//! iflag: return code of the solver (0 on success, negative on error, positive on warning)
//!

// First a forward decleration of a helper function to compute the residual
void SUBR(computeResidual)(TYPE(const_op_ptr) B_op, TYPE(mvec_ptr) r_ptr,
        TYPE(mvec_ptr) Au_ptr, TYPE(mvec_ptr) u_ptr, TYPE(mvec_ptr) rtil_ptr,
        TYPE(mvec_ptr) Qv, TYPE(mvec_ptr) tmp, TYPE(sdMat_ptr) Theta,
        TYPE(sdMat_ptr) atil, TYPE(sdMat_ptr) *atilv, _MT_ *resid,
        int nv, int nconv, int* iflag);

// Now the actual main function
void SUBR(jdqr)(TYPE(const_op_ptr) A_op, TYPE(const_op_ptr) B_op,
                TYPE(mvec_ptr) X, TYPE(mvec_ptr) Qout, TYPE(sdMat_ptr) Rout,
                _ST_* evals, _MT_* resid, int* is_cmplx,
        phist_jadaOpts_t opts, int* num_eigs, int* num_iters,
        int* iflag)
  {
  PHIST_ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"

  int nconv=0; // number of converged eigenpairs

  // partial QR decomposition (converged eigenvectors form Q)
  mvec_ptr_t Qv=NULL;  // will be created as a view of the converged eigenvectors
  mvec_ptr_t Qtil=NULL; // for representing [Q u] with u the current approximation
  sdMat_ptr_t R=NULL;      // store R part of QR decomp
  sdMat_ptr_t Rv=NULL; // to view the first nconv x nconv block of R
  
  // current basis 
  mvec_ptr_t V=NULL, AV=NULL;
  // temporary storage
  mvec_ptr_t Vtmp=NULL;
  // current update
  mvec_ptr_t t=NULL;
  // some more (single) vectors we need:
  mvec_ptr_t u=NULL, Au=NULL, r=NULL, rtil=NULL;
#ifndef IS_COMPLEX
  // the above vectors may be complex even for a real   
  // matrix. These are the first columns to be used for 
  // real eigenvectors                                  
  // (in the complex case x_r aliases x and x_i==NULL). 
  mvec_ptr_t    u_r=NULL, Au_r=NULL, r_r=NULL, rtil_r=NULL, t_r=NULL;
#endif
  // these point either to u or u_r etc. to make live simpler further down
  mvec_ptr_t u_ptr=NULL, Au_ptr=NULL, r_ptr=NULL, rtil_ptr=NULL, t_ptr=NULL;
  // Q*s (temporary vector)
  sdMat_ptr_t atil=NULL, atilv=NULL;
  // matrix <V,AV>
  sdMat_ptr_t M=NULL;
  sdMat_ptr_t Mv=NULL; // to create views of parts of M
  // Schur-decomposition of M
  sdMat_ptr_t T=NULL, S=NULL;
  sdMat_ptr_t Tv=NULL,Sv=NULL; // to create views of parts of T and S

  // for extracting views of M, T and S (to call lapack etc.)
  lidx_t ldM,ldS,ldT;
  ST *M_raw, *S_raw, *T_raw;
  
  //! view of certain columns of V
  mvec_ptr_t Vv=NULL, Vm=NULL;
  //! view of certain columns of AV
  mvec_ptr_t AVv=NULL, AVm=NULL;

  //! Hessenberg matrix from the initial Arnoldi iteration
  sdMat_ptr_t H0=NULL;
  
  int expand=1;
  int solve=1;
  std::complex<MT> theta; // next eigenvalue to go for
  // theta as a sdMat_t (either 1x1 or 2x2 in real case)
  sdMat_ptr_t Theta=NULL;
  
  int numEigs=opts.numEigs;
  MT tol = opts.convTol;
  int maxIter=opts.maxIters;
  eigSort_t which=opts.which;
  eigExtr_t how=opts.how;
  int minBas = opts.minBas;
  int maxBas = opts.maxBas;
  linSolv_t innerSolvType = opts.innerSolvType;
  int innerSolvMaxBas = opts.innerSolvMaxBas;
  int innerSolvBlockSize = opts.innerSolvBlockSize;
  bool arno=(bool)opts.arno;
  CT initialShift=(ST)opts.initialShift_r
                 +(ST)opts.initialShift_i*st::cmplx_I();

  std::complex<MT> ev[maxBas+1];
  
  int i,it,m,mm;
  MT res_nrm; // residual norm
  MT prev_nrm; // previous norm, to detect wether we're still working on the same eigenvalue

  // set output format for floating point numbers
  std::cout << std::scientific << std::setprecision(4) << std::setw(8);
  
  if (how!=STANDARD)
  {
    PHIST_SOUT(PHIST_ERROR,"only Ritz extraction is implemented (jadaOpts.how=%s), found %s\n",
        eigExtr2str(STANDARD),eigExtr2str(how));
    *iflag=PHIST_NOT_IMPLEMENTED;
    return;
  }

  *num_iters=0;
  *num_eigs=0;
  res_nrm=1.0e20;// some random large value 
  const_comm_ptr_t comm;
  PHIST_CHK_IERR(phist_map_get_comm(A_op->range_map,&comm,iflag),*iflag);

#ifdef IS_COMPLEX
  const int nv_max=1;
#else
  const int nv_max=2;
#endif

  _ST_ sigma[nv_max];

  // In the real case we allow the basis to grow up to maxBas+1 so we can add a complex    
  // pair even if a restart would be required strictly speaking (nv_max=2)
  int maxVecs=maxBas+nv_max-1;

  // we need maxBas vectors to store the maximum size JaDa basis, but we add some temporary
  // storage so that we can store [V,Q], where Q contains the already converged eigenspace.
  PHIST_CHK_IERR(SUBR(mvec_create)(&V,A_op->domain_map,maxVecs+numEigs,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_create)(&AV,A_op->domain_map,maxVecs,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_create)(&Vtmp,A_op->domain_map,maxVecs,iflag),*iflag);

  PHIST_CHK_IERR(SUBR(mvec_create)(&u,A_op->domain_map,nv_max,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_create)(&Au,A_op->domain_map,nv_max,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_create)(&r,A_op->domain_map,nv_max,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_create)(&rtil,A_op->domain_map,nv_max,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_create)(&t,A_op->domain_map,nv_max,iflag),*iflag);
  
#ifndef IS_COMPLEX
  PHIST_CHK_IERR(SUBR(mvec_view_block)(u,&u_r,0,0,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_view_block)(Au,&Au_r,0,0,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_view_block)(r,&r_r,0,0,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_view_block)(rtil,&rtil_r,0,0,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_view_block)(t,&t_r,0,0,iflag),*iflag);
#endif
  PHIST_CHK_IERR(SUBR(sdMat_create)(&M,maxVecs,maxVecs,comm,iflag),*iflag);
  // these two are made bigger because they are used as temporary storage when 
  // orthogonalizing against [V Q], which is of dimension up to n x (maxBas+numEigs)
  PHIST_CHK_IERR(SUBR(sdMat_create)(&S,maxVecs+numEigs,maxVecs+numEigs,comm,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_create)(&T,maxVecs+numEigs,maxVecs+numEigs,comm,iflag),*iflag);

  if (Rout == NULL)
  {
  PHIST_CHK_IERR(SUBR(sdMat_create)(&R,numEigs+1,numEigs+1,comm,iflag),*iflag);
  }
  else
  {
    R = Rout;
  }
  PHIST_CHK_IERR(SUBR(sdMat_put_value)(R,st::zero(),iflag),*iflag);

  PHIST_CHK_IERR(SUBR(sdMat_create)(&atil,maxVecs,nv_max,comm,iflag),*iflag);

  // pointer to the data in M, S and T
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(M,&M_raw, &ldM,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(S,&S_raw, &ldS,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(T,&T_raw, &ldT,iflag),*iflag);

  // view the first minBas+1 columns of V and H as M(1:minBas+1,1:minBas)
  PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&Vv,0,minBas,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_view_block)(M,&H0, 0,minBas,0,minBas-1,iflag),*iflag);

  // setup inner solver
  TYPE(jadaCorrectionSolver_ptr) innerSolv = NULL;

  // print parameters passed in to the method
  PHIST_SOUT(PHIST_VERBOSE,"====================\n");
  PHIST_SOUT(PHIST_VERBOSE,"| JDQR parameters  |\n");
  PHIST_SOUT(PHIST_VERBOSE,"====================\n");
  PHIST_SOUT(PHIST_VERBOSE,"#eigs\t%d\n",numEigs);
  PHIST_SOUT(PHIST_VERBOSE,"which\t%s\n",eigSort2str(which));
  PHIST_SOUT(PHIST_VERBOSE,"tol\t%5.3g\n",tol);
  PHIST_SOUT(PHIST_VERBOSE,"#iter\t%d\n",maxIter);
  PHIST_SOUT(PHIST_VERBOSE,"minBas\t%d\n",minBas);
  PHIST_SOUT(PHIST_VERBOSE,"maxBas\t%d\n",maxBas);
  PHIST_SOUT(PHIST_VERBOSE,"innerSolvType\t%s\n",linSolv2str(innerSolvType));
  if (!arno)
  {
    PHIST_SOUT(PHIST_VERBOSE,"initialShift\t%4.2g%+4.2g\n",(MT)ct::real(initialShift),(MT)ct::imag(initialShift));
  }
  PHIST_SOUT(PHIST_VERBOSE,"====================\n");

  mvec_ptr_t v0=(mvec_ptr_t)opts.v0;
  //TODO - we promised to use the user-provided subspace v0, but here we
  //       use only the first vector as start vec when Arnoldi is used.
  if (v0==NULL)
  {
    PHIST_CHK_IERR(SUBR(mvec_view_block)(X,&v0,0,0,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_random)(v0,iflag),*iflag);
    MT nrmV0;
    PHIST_CHK_IERR(SUBR(mvec_normalize)(v0,&nrmV0,iflag),*iflag);
  }
  else if (arno)
  {
    // only use the first column of the given mvec as v0
    int nv0;
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(v0,&nv0,iflag),*iflag);
    if (nv0>1)
    {
      mvec_ptr_t V0=v0;
      PHIST_CHK_IERR(SUBR(mvec_view_block)(V0,&v0,0,0,iflag),*iflag);
    }
  }

  if (arno)
  {
    PHIST_SOUT(1,"%d steps of Arnoldi as start-up\n",minBas);
    PHIST_CHK_IERR(SUBR(simple_arnoldi)(A_op,NULL,v0,Vv,NULL,NULL,H0,minBas,iflag),*iflag);

#if PHIST_OUTLEV>=PHIST_DEBUG
    PHIST_DEB("initial H (by Arnoldi)\n");
    PHIST_CHK_IERR(SUBR(sdMat_print)(H0,iflag),*iflag);
#endif

    PHIST_CHK_IERR(SUBR(mvec_view_block)(AV,&AVv,0,minBas-1,iflag),*iflag);
    // set Vv=V(:,1:minBas)
    PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&Vv,0,minBas-1,iflag),*iflag);
    // compute A*V for the first minBas columns which we have now using the Arnoldi relation
    //PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(st::one(),Vv,H0,st::zero(),AVv,iflag),*iflag);
    // fails if orthogonal basis was extended by random vectors
    PHIST_CHK_IERR(A_op->apply(st::one(),A_op->A,Vv,st::zero(),AVv,iflag),*iflag);

    // compute initial projection M = V'AV
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(M,&Mv,0,minBas-1,0,minBas-1,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),Vv,AVv,st::zero(),Mv,iflag),*iflag);

    m=minBas; // number of vectors in current basis V
    PHIST_CHK_IERR(SUBR(sdMat_delete)(H0,iflag),*iflag);
    expand=0;
  }
  else
  {
    // start Jacobi-Davidson directly with a given shift
#ifdef IS_COMPLEX
    t_ptr=t;
#else
    t_ptr=t_r;
#endif

    int nv0;
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(v0,&nv0,iflag),*iflag);
    t_ptr = v0;
    m=0;
    expand=nv0;
  }

  // delete the views (not the data) v0 and H0
  if (v0!=opts.v0)
  {
    PHIST_CHK_IERR(SUBR(mvec_delete)(v0,iflag),*iflag);
  }

  PHIST_SOUT(PHIST_INFO,"Start Jacobi-Davidson\n");

  it=0; // total number of JaDa iterations
  mm=0; //count iterations per eigenvalue/pair of complex ev's
  if (!arno) mm--; // to get the inner conv tol right we need this hack
  int nv=1; // number of vectors in current update (1 or 2 in this implementation, 2 for
            // complex eigenvectors of a real matrix)
  
  //ALG while (nconv<numEigs && it < maxIter)
  while (nconv<numEigs && it < maxIter)
  {
    //ALG if (expand)
    if (expand)
    {
      PHIST_DEB("EXPAND basis by %d vector(s)\n",nv);
      int m0=m;
      m+=expand;
      mm++;
      it++;
      *num_iters=it;
  
      // create views of V(:,1:m-1), V(:,m) and AV likewise
      
      // V temporarily gets extra storage to store [V Q] for the orthogonalization
      int ncVQ=m0+nconv;// number of columns of [V,Q]
      PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&Vm,m0,m-1,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(mvec_view_block)(AV,&AVm,m0,m-1,iflag),*iflag);

      if (m0>0)
      {
        PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&Vv,0,ncVQ-1,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(mvec_view_block)(AV,&AVv,0,m0-1,iflag),*iflag);
        // orthogonalize t against V(:,0:m-1) and the converged eigenspace Q.
        // We use T and S as temporary storage here
        PHIST_CHK_IERR(SUBR(sdMat_view_block)(T,&Tv,0,nv-1,0,nv-1,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(sdMat_view_block)(S,&Sv,0,ncVQ-1,0,nv-1,iflag),*iflag);
        if (nconv>0)
          {
          PHIST_CHK_IERR(SUBR(mvec_set_block)(Vv,Qv,m0,ncVQ-1,iflag),*iflag);
          }
        int rankV;
        //ALG orthogonalize t against [Q V(:,0:m-1)]
        PHIST_CHK_NEG_IERR(SUBR(orthog)(Vv,t_ptr,B_op,Tv,Sv,3,&rankV,iflag),*iflag);
        // reset the view of V without the Q.
        PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&Vv,0,m0-1,iflag),*iflag);
      }
      //ALG set V(:,m)=t
      PHIST_CHK_IERR(SUBR(mvec_set_block)(V,t_ptr,m0,m-1,iflag),*iflag);
      // ALG compute AV(:,m) = A*V(:,m)
      PHIST_CHK_IERR(A_op->apply(st::one(),A_op->A,Vm,st::zero(),AVm,iflag),*iflag);

      // Galerkin for non-Hermitian A
      PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&Vv,0,m-1,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(mvec_view_block)(AV,&AVv,0,m-1,iflag),*iflag);
      //ALG compute M = V'A*V
      PHIST_CHK_IERR(SUBR(sdMat_view_block)(M,&Mv,0,m-1,0,m-1,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),Vv,AVv,st::zero(),Mv,iflag),*iflag);

      PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&Vv,0,m-1,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(mvec_view_block)(AV,&AVv,0,m-1,iflag),*iflag);
      PHIST_DEB("basis size is now m=%d\n",m);
#ifdef TESTING
      // check orthogonality of [V Q]
      sdMat_ptr_t tmp1=NULL,tmp2=NULL;
      ST *tmp1_raw, *tmp2_raw;
      lidx_t ld1,ld2;
      PHIST_CHK_IERR(SUBR(sdMat_create)(&tmp1,m,m,comm,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_extract_view)(tmp1,&tmp1_raw,&ld1,iflag),*iflag);
        int ncV;
        PHIST_CHK_IERR(SUBR(mvec_num_vectors)(Vv,&ncV,iflag),*iflag);
        PHIST_SOUT(PHIST_VERBOSE,"#vectors in V: %d\n",ncV);
      PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),Vv,Vv,st::zero(),tmp1,iflag),*iflag);
      if (nconv>0)
      {
        int ncQ;
        PHIST_CHK_IERR(SUBR(mvec_num_vectors)(Qv,&ncQ,iflag),*iflag);
        PHIST_SOUT(PHIST_VERBOSE,"#vectors in Q: %d\n",ncQ);
        PHIST_CHK_IERR(SUBR(sdMat_create)(&tmp2,m,nconv,comm,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(sdMat_extract_view)(tmp2,&tmp2_raw,&ld2,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),Vv,Qv,st::zero(),tmp2,iflag),*iflag);
      }
      MT err1=mt::zero(), err2=mt::zero();
      for (int i=0;i<m;i++)
      {
        for (int j=0;j<m;j++)
        {
          if (i!=j)
          {
            err1+=std::abs(tmp1_raw[j*ld1+i]);
          }
        }
        err1+=mt::one()-std::abs(tmp1_raw[i*ld1+i]);// diagonal element should be 1
        if (nconv>0)
        {
          for (int j=0;j<nconv;j++)
          {
            err2+=std::abs(tmp2_raw[j*ld2+i]);
          }
        }
      }
      PHIST_SOUT(PHIST_VERBOSE,"orthogonality monitors:\n");
      PHIST_SOUT(PHIST_VERBOSE,"||V'V-I||=%e\n",err1);
      PHIST_SOUT(PHIST_VERBOSE,"||V'Q||=%e\n",err2);
      PHIST_CHK_IERR(SUBR(sdMat_delete)(tmp1,iflag),*iflag);
      if (nconv>0)
      {
        PHIST_CHK_IERR(SUBR(sdMat_delete)(tmp2,iflag),*iflag);
      }
#endif
    }//ALG else
    else
    {
      //ALG expand=1
      expand=1;// the expand flag is just used at start-up
               // for starting with a vector or some Arnoldi steps.
    }
    //ALG end if
      
    // now do a Schur decomposition of M into S and T, with the
    // Ritz values on the diagonal of T.

    // copy T=M. Note that as we just 'view' the upper left blocks,
    // the pointers M_raw, T_raw, S_raw are still correct
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(M,&Mv,0,m-1,0,m-1,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(T,&Tv,0,m-1,0,m-1,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_get_block)(Mv,Tv,0,m-1,0,m-1,iflag),*iflag);

    // if a restart is impending, shift the next <minBas> Ritz values
    // to the top left. Shift the next <nsort> of them to the 
    // top left (next shift(s) to choose) in sorted order.
    //int nselect = m==maxBas-1? minBas: 1;
    //int nsort = 1;
    
    // TODO - check out how much sorting is required, in our matlab prototype
    // we have to sort the entire thing here because of the way we select the 
    // converged eigenvalues. In JDQR Sleijpen just sorts as stated above.
    // Unnecessary sorting should be avoided because it may change the value
    // of ill-conditioned Ritz values.
    int nselect=m;
    int nsort=nselect;
    //ALG sorted Schur decomposition
    PHIST_CHK_IERR(SUBR(SchurDecomp)
        (T_raw,ldT,S_raw,ldS,m,nselect,nsort,which,tol,ev,iflag),*iflag);
    int ev_pos=0;
    theta=ev[ev_pos];

#if PHIST_OUTLEV>=PHIST_DEBUG
  PHIST_DEB("it=%d, m=%d, mm=%d\n",it,m,mm);
  PHIST_DEB("projected matrix\n");
  PHIST_CHK_IERR(SUBR(sdMat_print)(Mv,iflag),*iflag);
  PHIST_DEB("sorted Schur form\n");
  PHIST_CHK_IERR(SUBR(sdMat_print)(Tv,iflag),*iflag);
#endif

    nv=1; // nv==2 indicates complex eigenvalue of real matrix, we
          // will use this flag later on.
    // set some pointers
    u_ptr=u;
    Au_ptr=Au;
    r_ptr=r;
    rtil_ptr=rtil;
    t_ptr=t;
#ifndef IS_COMPLEX
    if (mt::abs(ct::imag(theta))>mt::eps())
    {
      PHIST_DEB("complex eigenvalue in real arithmetic\n");
      nv++; // get imaginary part
    }
    else
    {
      // real matrix, real Ritz theta
      u_ptr=u_r;
      Au_ptr=Au_r;
      r_ptr=r_r;
      rtil_ptr=rtil_r;
      t_ptr=t_r;
    }
#endif

    // get the diagonal block (1x1 or 2x2) corresponding to theta
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(T,&Theta,0,nv-1,0,nv-1,iflag),*iflag);

    // get Ritz vector corresponding to theta
    // (for complex RV in real case, two columns are extracted)
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(S,&Sv,0,m-1,0,nv-1,iflag),*iflag);
    //u=V*s;
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(st::one(),Vv,Sv,st::zero(),u_ptr,iflag),*iflag);
    //Au=AV*s;
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(st::one(),AVv,Sv,st::zero(),Au_ptr,iflag),*iflag);

    prev_nrm=res_nrm;

    PHIST_CHK_IERR(SUBR(computeResidual)(B_op, r_ptr, Au_ptr, u_ptr,
        rtil_ptr, Qv, Vtmp, Theta, atil, &atilv, &res_nrm, nv, nconv, iflag),*iflag);

    PHIST_SOUT(PHIST_INFO,"JDQR Iter %d\tdim(V)=%d\ttheta=%8.4g%+8.4gi\t\tr_est=%8.4e\n",it,m,ct::real(theta),ct::imag(theta),res_nrm);
    // deflate converged eigenpairs
    if (res_nrm<=tol)
    {
      mm=0;// counts iterations for next theta
      // number of converged eigenpairs: +2 for complex conjugate pairs in real case
      int nq0=nconv;
      nconv=nconv+nv;
      *num_eigs=nconv;// tell the user how many we return
      for (int j=nq0;j<nconv;j++)
      {
        resid[j]=res_nrm; // will be returned to the user
      }

      // view first nconv columns of X as 'Q'
      PHIST_CHK_IERR(SUBR(mvec_view_block)(X,&Qv,0,nconv-1,iflag),*iflag);

       //Q=[Q,u];
       PHIST_CHK_IERR(SUBR(mvec_set_block)(Qv,u_ptr,nq0,nconv-1,iflag),*iflag);
       //R=[R, atil;
       //   0, theta];
       if (nq0>0)
       {
         PHIST_CHK_IERR(SUBR(sdMat_set_block)(R,atilv,0,nq0-1,nq0,nconv-1,iflag),*iflag);
       }
       PHIST_CHK_IERR(SUBR(sdMat_set_block)(R,Theta,nq0,nconv-1,nq0,nconv-1,iflag),*iflag);

      PHIST_SOUT(PHIST_VERBOSE,"eigenvalue %d (%8.4g%+8.4gi) is converged.\n",
          nconv,ct::real(theta),ct::imag(theta));
      if (nconv>=numEigs)
      {
        PHIST_SOUT(PHIST_VERBOSE,"stopping JDQR loop\n");
        solve=false;
        break;
      }

      //S=S(:,2:m); (select remaining Ritz vectors)
      PHIST_CHK_IERR(SUBR(sdMat_view_block)(S,&Sv,0,m-1,nv,m-1,iflag),*iflag);
      //M=T(2:m,2:m);
      // let Mblock point to the first m-nv x m-nv block of M
      PHIST_CHK_IERR(SUBR(sdMat_view_block)(M,&Mv,0,m-nv-1,0,m-nv-1,iflag),*iflag);
      // copy T(nv+1:m,nv+1:m) into it
      PHIST_CHK_IERR(SUBR(sdMat_get_block)(T,Mv,nv,m-1,nv,m-1,iflag),*iflag);
      m=m-nv;
      //V=V*S;
      // It is not allowed to alias V in mvec_times_sdMat, so we use a temporary vector
      mvec_ptr_t v_tmp=NULL;
      PHIST_CHK_IERR(SUBR(mvec_view_block)(Vtmp,&v_tmp,0,m-1,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(st::one(),Vv,Sv,st::zero(),v_tmp,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(mvec_set_block)(V,v_tmp,0,m-1,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&Vv,0,m-1,iflag),*iflag);
      //AV=AV*S;
      PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(st::one(),AVv,Sv,st::zero(),v_tmp,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(mvec_set_block)(AV,v_tmp,0,m-1,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(mvec_view_block)(AV,&AVv,0,m-1,iflag),*iflag);
    
      // delete the temporary view (not the data, of course)
      PHIST_CHK_IERR(SUBR(mvec_delete)(v_tmp,iflag),*iflag);

      // select next target ev. This need not be the next according to
      // the sort criterion because the Schur-form is not completely sorted
      // (TODO: is it desirable to do them in order?)

      ev_pos+=nv; // skip conjugate ev for real case w/ cmplx theta
      theta=ev[ev_pos];
      nv=1;
#ifndef IS_COMPLEX
      if (mt::abs(ct::imag(theta))>mt::eps())
      {
        nv++;
        u_ptr=u;
        Au_ptr=Au;
        r_ptr=r;
        rtil_ptr=rtil;
        t_ptr=t;
      }
      else
      {
        u_ptr=u_r;
        Au_ptr=Au_r;
        r_ptr=r_r;
        rtil_ptr=rtil_r;
        t_ptr=t_r;
      }
#endif

      //u=V(:,1);
      PHIST_CHK_IERR(SUBR(mvec_get_block)(Vv,u_ptr,0,nv-1,iflag),*iflag);
      //Au=AV(:,1);
      PHIST_CHK_IERR(SUBR(mvec_get_block)(AVv,Au_ptr,0,nv-1,iflag),*iflag);

      // get the diagonal block (1x1 or 2x2) corresponding to theta
      PHIST_CHK_IERR(SUBR(sdMat_view_block)(T,&Theta,ev_pos,ev_pos+nv-1,ev_pos,ev_pos+nv-1,iflag),*iflag);

      PHIST_CHK_IERR(SUBR(computeResidual)(B_op, r_ptr, Au_ptr, u_ptr,
          rtil_ptr, Qv, Vtmp, Theta, atil, &atilv, &res_nrm, nv, nconv, iflag),*iflag);

      // again, sort the largest Ritz value to the top
//    nselect = m==maxBas-1? minBas: 1;
//    nsort = 1;
      nselect=m;
      nsort=nselect;
    // set T=M
      PHIST_CHK_IERR(SUBR(sdMat_set_block)(T,Mv,0,m-1,0,m-1,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(SchurDecomp)
          (T_raw,ldT,S_raw,ldS,m,nselect,nsort,which,tol,ev,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_view_block)(T,&Tv,0,m-1,0,m-1,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_view_block)(S,&Sv,0,m-1,0,m-1,iflag),*iflag);
    } //if (deflate)

    if (res_nrm>100*prev_nrm)
    {
      // we probably aren't working on the same
      // eigenvalue anymore, so we reset the convergence
      // tolerance
    }

    // restart if necessary
    //ALG if (m>=maxBas restart
    if (m>=maxBas)
    {
      int m0=m;
      m=minBas;

      PHIST_SOUT(PHIST_VERBOSE,"restart JDQR, shrink basis from %d to %d vectors\n",m0,m);
    
      //S=S(:,1:m);
      PHIST_CHK_IERR(SUBR(sdMat_view_block)(S,&Sv,0,m0-1,0,m-1,iflag),*iflag);
      //M=T(1:m,1:m);
      PHIST_CHK_IERR(SUBR(sdMat_view_block)(M,&Mv,0,m-1,0,m-1,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_get_block)(T,Mv,0,m-1,0,m-1,iflag),*iflag);
      //V=V*S;
      mvec_ptr_t v_tmp=NULL;
      PHIST_CHK_IERR(SUBR(mvec_view_block)(Vtmp,&v_tmp,0,m-1,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(st::one(),Vv,Sv,st::zero(),v_tmp,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(mvec_set_block)(V,v_tmp,0,m-1,iflag),*iflag);
      //AV=AV*S;
      PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(st::one(),AVv,Sv,st::zero(),v_tmp,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(mvec_set_block)(AV,v_tmp,0,m-1,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(mvec_delete)(v_tmp,iflag),*iflag);
      //V=V(:,1:mmin);
      PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&Vv,0,m-1,iflag),*iflag);
      //AV=AV(:,1:mmin);
      PHIST_CHK_IERR(SUBR(mvec_view_block)(AV,&AVv,0,m-1,iflag),*iflag);
      solve=true;
      expand=nv;
    }
    //ALG end if
    //ALG if (solve)
    if (solve)
    {
      if (nv>1)
      {
        PHIST_SOUT(PHIST_VERBOSE,"complex eigenvalue of real matrix, using real GMRES right now.\n");
      }
      // maintain orthogonality against
      // all converged vectors (Q) and the
      // new one u:
      //Qtil=[Q,u];
      PHIST_CHK_IERR(SUBR(mvec_view_block)(X,&Qtil,0,nconv+nv-1,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(mvec_set_block)(Qtil,u_ptr,nconv,nconv+nv-1,iflag),*iflag);
      CT actual_shift=theta;
      if (m<minBas && res_nrm>1.0e-2)
      {
        PHIST_SOUT(PHIST_VERBOSE,"start-up step with fixed sigma=%4.2f%+4.2f\n",
        (MT)ct::real(initialShift),(MT)ct::imag(initialShift));
        actual_shift=(CT)initialShift; // start-up without Arnoldi
      }
#ifdef IS_COMPLEX
      sigma[0] = actual_shift;
#else
      sigma[0] = ct::real(actual_shift);
      sigma[1] = ct::real(actual_shift);
#endif
      // solve approximately 
      // (I-uu')(A-theta*I)(I-uu')*t=-r
      // to get t \orth u (u ^= Qtil here)
      // In the real case with complex Ritz value theta,
      // we solve the update equation in real arithmetic with
      // real shift. This results in a linear system with two
      // right-hand sides.

      // the block case is not implemented here:
      // (I-uu')A(I-uu')*t - (I-uu')t*Theta=-r

      // 1/2^mm, but at most the outer tol as conv tol for GMRES
      MT innerTol[2];
      innerTol[0] = std::max(tol,mt::one()/((MT)(pow(2.0,mm+1))));
      innerTol[1] = std::max(tol,mt::one()/((MT)(pow(2.0,mm+1))));
      PHIST_SOUT(PHIST_VERBOSE,"inner conv tol: %g\n",innerTol[0]);

      opts.innerSolvBlockSize=nv;
      PHIST_CHK_IERR(SUBR(jadaCorrectionSolver_create)(&innerSolv, opts,
          A_op->domain_map, iflag), *iflag);


      // allow at most the given number of iterations
      int nIt=opts.innerSolvMaxIters;
      PHIST_CHK_NEG_IERR(SUBR(jadaCorrectionSolver_run)(innerSolv, A_op, B_op, Qtil, NULL, 
                                                      sigma, rtil_ptr, NULL, 
                                                      innerTol, nIt, t_ptr, opts.innerSolvRobust, false, iflag), *iflag);
      expand=nv;
      PHIST_CHK_IERR(SUBR(jadaCorrectionSolver_delete)(innerSolv,iflag),*iflag);
      innerSolv=NULL;
    }//ALG else
    else
    {
      expand=0;
      solve=true;
    }
    //ALG end if
  }
  //ALG end while
    
  // TODO Sleijpen in his jdqr code does a refinement step at this      
  // point taking the information from the current V into Q and R       
  // (might be a nice feature). It also allows us to return unconverged 
  // eigenvectors so the user can restart from them.                    
  // In that case the value of the residual norm in resid should be     
  // considered an upper bound (or updated, of course).                 
  PHIST_TOUCH(Rv); // we would need it for this
  
  // some sanity checks - the user provides only num_eigs+1 slots
  nconv = std::min(numEigs+1,nconv);
  
  if (nconv>0)
    {
    // copy at most num_eigs+1 converged eigenvalues into the user
    // provided array
    ST* R_raw=NULL;
    lidx_t ldR;
    PHIST_CHK_IERR(SUBR(sdMat_extract_view)(R,&R_raw,&ldR,iflag),*iflag);
    i=0;
    while (i<nconv)
    {
      evals[i] = R_raw[i*ldR+i];
#ifndef IS_COMPLEX
      if (i<nconv-1)
      {
        if (st::abs(R_raw[i*ldR+i+1])>mt::eps())
        {
          is_cmplx[i]=1;
          // sqrt((T^2)/4-D) gives the imaginary part of the eigenpair of
          // a 2x2 matrix, with T the trace and D the determinant. In our
          // case A11=A22 => im(lambda)=im(sqrt(A12*A21))
          evals[i+1]=st::sqrt(-R_raw[i*ldR+i+1]*R_raw[(i+1)*ldR+i]);
          is_cmplx[i+1]=1;
          i++;
        }
        else
        {
            is_cmplx[i]=0;
        }
      }
      else
      {
          is_cmplx[i]=0;
      }
#else
      PHIST_TOUCH(is_cmplx);
#endif    
      i++;
    }
  
    // compute eigenvectors. Given the Schur form R, first compute all its 
    // eigenvectors in S. Then compute the eigenvectors of A as Q*S.
    const char* side="R";
    const char* howmny="A";
    int m_out;
    MT* work=new MT[5*nconv];
    blas_idx_t ildS=static_cast<blas_idx_t>(ldS);
    blas_idx_t ildR=static_cast<blas_idx_t>(ldR);
    if (ildS<0 || ildR<0)
    {
    *iflag=PHIST_INTEGER_OVERFLOW;
    return;
    }
#ifdef IS_COMPLEX
    MT* rwork = work+4*nconv;
    PHIST_CHK_IERR(PREFIX(TREVC)((blas_char_t*)side, (blas_char_t*)howmny, NULL, 
    &nconv, (mt::blas_cmplx_t*)R_raw, &ildR, NULL, &ildS, 
    (mt::blas_cmplx_t*)S_raw, &ildS, &nconv, &m_out, (mt::blas_cmplx_t*)work, 
    rwork, iflag),*iflag);
#else
    PHIST_CHK_IERR(PREFIX(TREVC)((blas_char_t*)side, (blas_char_t*)howmny, NULL, 
    &nconv,R_raw, &ildR, NULL, &ildS, S_raw, &ildS, &nconv, &m_out, work, 
    iflag),*iflag);
#endif  
    delete [] work;
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(S,&Sv,0,nconv-1,0,nconv-1,iflag),*iflag);
    TYPE(mvec_ptr) Qcopy=NULL;
    if (Qout == NULL)
    {
    PHIST_CHK_IERR(SUBR(mvec_view_block)(Vtmp,&Qcopy,0,nconv-1,iflag),*iflag);
    }
    else
    {
      PHIST_CHK_IERR(SUBR(mvec_view_block)(Qout,&Qcopy,0,nconv-1,iflag),*iflag);
    }
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),Qv,st::zero(),Qcopy,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(st::one(),Qcopy,Sv,st::zero(),Qv,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_delete)(Qcopy,iflag),*iflag);
    }// any eigenpairs converged?

  PHIST_CHK_IERR(SUBR(mvec_delete)(V,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_delete)(AV,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_delete)(Vtmp,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_delete)(u,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_delete)(Au,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_delete)(r,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_delete)(rtil,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_delete)(t,iflag),*iflag);
  
  PHIST_CHK_IERR(SUBR(sdMat_delete)(M,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(S,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(T,iflag),*iflag);
  if (Rout == NULL)
  {
  PHIST_CHK_IERR(SUBR(sdMat_delete)(R,iflag),*iflag);
  }
  
  // we also call delete for the views to avoid small memory leaks
#ifndef IS_COMPLEX
  PHIST_CHK_IERR(SUBR(mvec_delete)(u_r,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_delete)(Au_r,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_delete)(r_r,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_delete)(rtil_r,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_delete)(t_r,iflag),*iflag);
#endif

  PHIST_CHK_IERR(SUBR(mvec_delete)(Vv,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_delete)(AVv, iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_delete)(Vm,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_delete)(AVm, iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_delete)(Qtil, iflag),*iflag);
  
  PHIST_CHK_IERR(SUBR(sdMat_delete)(Mv, iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(Tv, iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(Sv, iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(Theta, iflag),*iflag);
  // these may still be NULL if nothing converged
  if (nconv>0)
    {
    PHIST_CHK_IERR(SUBR(mvec_delete)(Qv, iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_delete)(atilv, iflag),*iflag);
    }
  return;
  }

void SUBR(computeResidual)(TYPE(const_op_ptr) B_op, TYPE(mvec_ptr) r_ptr,
        TYPE(mvec_ptr) Au_ptr, TYPE(mvec_ptr) u_ptr, TYPE(mvec_ptr) rtil_ptr,
        TYPE(mvec_ptr) Qv, TYPE(mvec_ptr) tmp, TYPE(sdMat_ptr) Theta,
        TYPE(sdMat_ptr) atil, TYPE(sdMat_ptr) *atilv, _MT_ *resid,
        int nv, int nconv, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  // r=Au-theta*u; ([r_r, r_i] = [Au_r, Au_i] - [u_r, u_i]*Theta in the real case with 
  // complex theta)

  MT nrm[2];

  // pointer to temporary storage for B*Q*atil if B is defined
  mvec_ptr_t tmp_ptr=NULL;

  // first: r=Au
  PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),Au_ptr,st::zero(),r_ptr,iflag),*iflag);

  if (B_op != NULL)
  {
    // rt = B * u
    PHIST_CHK_IERR(B_op->apply(st::one(),B_op->A,u_ptr,st::zero(),rtil_ptr,iflag),*iflag);

    // update r = r - rt*Theta
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(-st::one(),rtil_ptr,Theta,st::one(),r_ptr,iflag),*iflag);
  }
  else
  {
    // update r = r - u*Theta
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(-st::one(),u_ptr,Theta,st::one(),r_ptr,iflag),*iflag);
  }

  // set rtil=r
  PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),r_ptr,st::zero(),rtil_ptr,iflag),*iflag);

  // project out already converged eigenvectors
  // TODO - we could use our orthog routine here instead
  if (nconv>0)
  {
    // view next ~a, a temporary vector to compute ~a=Q'*r
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(atil,atilv,0,nconv-1,0,nv-1,iflag),*iflag);

    //atil = Q'*r;
    PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),Qv,r_ptr,st::zero(),*atilv,iflag),*iflag);

    if (B_op != NULL)
    {
      PHIST_CHK_IERR(SUBR(mvec_view_block)(tmp,&tmp_ptr,0,nv-1,iflag),*iflag);
      // tmp = Q*atil
      PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(st::one(),Qv,*atilv,st::zero(),tmp_ptr,iflag),*iflag);

      //rtil = r-B*Q*atil;
      PHIST_CHK_IERR(B_op->apply(-st::one(),B_op->A,tmp_ptr,st::one(),rtil_ptr,iflag),*iflag);

      PHIST_CHK_IERR(SUBR(mvec_delete)(tmp_ptr,iflag),*iflag);
    }
    else
    {
      //rtil = r-Q*atil;
      PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(-st::one(),Qv,*atilv,st::one(),rtil_ptr,iflag),*iflag);
    }
  }

  //nrm=norm(rtil);
  // real case with complex r: ||v+iw||=sqrt((v+iw).'*(v-iw))=sqrt(v'v+w'w).
  // in the complex case we pass in a 're-interpret cast' of nrm as complex,
  // which should be fine (imaginary part will be 0).
  nrm[1]=mt::zero();
  PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(rtil_ptr,rtil_ptr,(ST*)nrm,iflag),*iflag);
  *resid=mt::sqrt(nrm[0]+nrm[1]);
}
