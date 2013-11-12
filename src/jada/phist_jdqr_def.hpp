//! tries to compute a given number of eigenpairs (num_eig) of 
//! an general operator A using the Jacobi-Davidson QR (JDQR)
//! method. We allow for generalized EVPs here, A*x=lambda*B*x,
//! where B should be positive definite.
//!
//! input arguments:
//!
//! A: pointer to the operator(-wrapper) A
//! B: pointer to the "mass-matrix" operator (if B==NULL, B=I is assumed)
//! X:  start vectors for the num_eigs desired                           
//!     eigenvectors (random vectors if you don't know better)           
//!     (currently ignored on input). Should be allocated with           
//!     at least <num_eigs> columns (num_eigs+1 in the real case         
//!     to avoid the situation where the last eigenvalue to con-         
//!     verge is a complex pair).                                        
//!                                                                      
//! which - decide at which end of the spectrum to look for eigenvalues
//! tol: convergence tolerance
//! num_eigs: number of desired eigenpairs
//! num_iters: maximum number of iterations allowed
//! is_cmplx,evals and resid should be pre-allocated arrays of size at least num_eigs
//! (in the real case num_eigs+1, cf. comment for 'X' above).
//! 
//! minBas: start up from a basis consisting of minBas vectors (using Arnoldi)
//! maxBas: when the basis reaches <maxBas> vectors, restart from <minBas> vectors.
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
//! ierr: return code of the solver (0 on success, negative on error, positive on warning)
//!
void SUBR(jdqr)(TYPE(const_op_ptr) A_op, TYPE(const_op_ptr) B_op,
        TYPE(mvec_ptr) X, _ST_* evals, _MT_* resid, int* is_cmplx,
        eigSort_t which, _MT_ tol, int* num_eigs, int* num_iters,
        int minBas, int maxBas,
        int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
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
#ifndef _IS_COMPLEX_
  // the above vectors may be complex even for a real   
  // matrix. These are the first columns to be used for 
  // real eigenvectors                                  
  // (in the complex case x_r aliases x and x_i==NULL). 
  mvec_ptr_t    u_r=NULL, Au_r=NULL, r_r=NULL, rtil_r=NULL, t_r=NULL;
#endif
  // these point either to u or u_r etc. to make live simpler further down
  mvec_ptr_t u_ptr=NULL, Au_ptr=NULL, r_ptr=NULL, rtil_ptr=NULL,t_ptr=NULL;
  // Q*s (temporary vector)
  sdMat_ptr_t atil=NULL, atilv=NULL;
  // matrix <V,AV>
  sdMat_ptr_t M=NULL;
  sdMat_ptr_t Mv=NULL; // to create views of parts of M
  // Schur-decomposition of M
  sdMat_ptr_t T=NULL, S=NULL;
  sdMat_ptr_t Tv=NULL,Sv=NULL; // to create views of parts of T and S

  // for extracting views of M, T and S (to call lapack etc.)
  int ldM,ldS,ldT;
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
  sdMat_ptr_t shift = NULL;
  // theta as a sdMat_t (either 1x1 or 2x2 in real case)
  sdMat_ptr_t Theta=NULL;
  std::complex<MT> ev[maxBas];
  
  int numEigs=*num_eigs;
  int maxIter=*num_iters;
  int i,it,m,mm;
  MT nrm[2]; // for computing residual norms

  // set output format for floating point numbers
  std::cout << std::scientific << std::setprecision(4) << std::setw(8);

  *num_iters=0;
  
  const_comm_ptr_t comm;
  PHIST_CHK_IERR(phist_map_get_comm(A_op->range_map,&comm,ierr),*ierr);

  // we need maxBas vectors to store the maximum size JaDa basis, but we add some temporary
  // storage so that we can store [V,Q], where Q contains the already converged eigenspace.
  PHIST_CHK_IERR(SUBR(mvec_create)(&V,A_op->domain_map,maxBas+numEigs,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_create)(&AV,A_op->domain_map,maxBas,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_create)(&Vtmp,A_op->domain_map,maxBas,ierr),*ierr);

#ifdef _IS_COMPLEX_
  int nv_max=1;
#else
  int nv_max=2;
#endif
  PHIST_CHK_IERR(SUBR(mvec_create)(&u,A_op->domain_map,nv_max,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_create)(&Au,A_op->domain_map,nv_max,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_create)(&r,A_op->domain_map,nv_max,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_create)(&rtil,A_op->domain_map,nv_max,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_create)(&t,A_op->domain_map,nv_max,ierr),*ierr);
  
#ifndef _IS_COMPLEX_
  u_r=NULL;
  PHIST_CHK_IERR(SUBR(mvec_view_block)(u,&u_r,0,0,ierr),*ierr);

  Au_r=NULL;
  PHIST_CHK_IERR(SUBR(mvec_view_block)(Au,&Au_r,0,0,ierr),*ierr);

  r_r=NULL;
  PHIST_CHK_IERR(SUBR(mvec_view_block)(r,&r_r,0,0,ierr),*ierr);

  rtil_r=NULL;
  PHIST_CHK_IERR(SUBR(mvec_view_block)(rtil,&rtil_r,0,0,ierr),*ierr);

  t_r=NULL;
  PHIST_CHK_IERR(SUBR(mvec_view_block)(t,&t_r,0,0,ierr),*ierr);
#endif
  PHIST_CHK_IERR(SUBR(sdMat_create)(&shift,1,1,NULL,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_create)(&M,maxBas,maxBas,comm,ierr),*ierr);
  // these two are made bigger because they are used as temporary storage when 
  // orthogonalizing against [V Q], which is of dimension up to n x (maxBas+numEigs)
  PHIST_CHK_IERR(SUBR(sdMat_create)(&S,maxBas+numEigs,maxBas+numEigs,comm,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_create)(&T,maxBas+numEigs,maxBas+numEigs,comm,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_create)(&R,*num_eigs,*num_eigs,comm,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_put_value)(R,st::zero(),ierr),*ierr);

  PHIST_CHK_IERR(SUBR(sdMat_create)(&atil,maxBas,nv_max,comm,ierr),*ierr);

  // pointer to the data in M, S and T
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(M,&M_raw, &ldM,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(S,&S_raw, &ldS,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(T,&T_raw, &ldT,ierr),*ierr);

  // view the first minBas+1 columns of V and H as M(1:minBas+1,1:minBas)
  PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&Vv,0,minBas,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_view_block)(M,&H0, 0,minBas,0,minBas-1,ierr),*ierr);

  // print parameters passed in to the method
  PHIST_OUT(PHIST_VERBOSE,"====================");
  PHIST_OUT(PHIST_VERBOSE,"| JDQR parameters  |");
  PHIST_OUT(PHIST_VERBOSE,"====================");
  PHIST_OUT(PHIST_VERBOSE,"#eigs\t%d",*num_eigs);
  PHIST_OUT(PHIST_VERBOSE,"which\t%s",eigSort2str(which));
  PHIST_OUT(PHIST_VERBOSE,"tol\t%4.2g",tol);
  PHIST_OUT(PHIST_VERBOSE,"#iter\t%d",maxIter);
  PHIST_OUT(PHIST_VERBOSE,"minBas\t%d",minBas);
  PHIST_OUT(PHIST_VERBOSE,"maxBas\t%d",maxBas);
  PHIST_OUT(PHIST_VERBOSE,"====================");


  PHIST_OUT(1,"%d steps of Arnoldi as start-up",minBas);
  
  // start by filling the first minBas vectors using Arnoldi
  //TODO - we promised to use the user input in X, but here we
  // use only the first vector as start vec. We could do block
  // Arnoldi instead, or use a linear combination of the X columns.
  // We're also assuming that v0 is normalized here.
  mvec_ptr_t v0=NULL;
  PHIST_CHK_IERR(SUBR(mvec_view_block)(X,&v0,0,0,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(simple_arnoldi)(A_op,NULL,v0,Vv,NULL,H0,minBas,ierr),*ierr);

#if PHIST_OUTLEV>=PHIST_DEBUG
  PHIST_DEB("initial H (by Arnoldi)");
  PHIST_CHK_IERR(SUBR(sdMat_print)(H0,ierr),*ierr);
#endif
  // delete the view (not the data) v0
  PHIST_CHK_IERR(SUBR(mvec_delete)(v0,ierr),*ierr);
  v0=NULL; // always important to nullify pointers in phist...
  PHIST_CHK_IERR(SUBR(mvec_view_block)(AV,&AVv,0,minBas-1,ierr),*ierr);
  // set Vv=V(:,1:minBas)
  PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&Vv,0,minBas-1,ierr),*ierr);
  // compute A*V for the first minBas columns which we have now using the Arnoldi relation
  //PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(st::one(),Vv,H0,st::zero(),AVv,ierr),*ierr);
  // fails if orthogonal basis was extended by random vectors
  PHIST_CHK_IERR(A_op->apply(st::one(),A_op->A,Vv,st::zero(),AVv,ierr),*ierr);

  // compute initial projection M = V'AV
  PHIST_CHK_IERR(SUBR(sdMat_view_block)(M,&Mv,0,minBas-1,0,minBas-1,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),Vv,AVv,st::zero(),Mv,ierr),*ierr);

  PHIST_OUT(PHIST_INFO,"Start Jacobi-Davidson");

  it=0; // total number of JaDa iterations
  mm=0; //count iterations per eigenvalue/pair of complex ev's
  m=minBas; // number of vectors in current basis V
  expand=0;
  int nv=1; // number of vectors in current update (1 or 2 in this implementation, 2 for
            // complex eigenvectors of a real matrix)
  
  while (nconv<numEigs && it < maxIter)
    {
    //TODO - thoroughly check the switch from 1-based m in MATLAB to 0-based here
    if (expand)
      {
      PHIST_DEB("EXPAND basis by %d vector(s)",nv);
      int m0=m;
      m+=expand;
      mm++;
      it++;
      // create views of V(:,1:m-1), V(:,m) and AV likewise
      
      // V temporarily gets extra storage to store [V Q] for the orthogonalization
      int ncVQ=m0+nconv;// number of columns of [V,Q]
    PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&Vv,0,ncVQ-1,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(mvec_view_block)(AV,&AVv,0,m0-1,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&Vm,m0,m-1,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(mvec_view_block)(AV,&AVm,m0,m-1,ierr),*ierr);
      // orthogonalize t against V(:,0:m-1) and the converged eigenspace Q.
      // We use T and S as temporary storage here
      PHIST_CHK_IERR(SUBR(sdMat_view_block)(T,&Tv,0,nv-1,0,nv-1,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(sdMat_view_block)(S,&Sv,0,ncVQ-1,0,nv-1,ierr),*ierr);
      if (nconv>0)
        {
        PHIST_CHK_IERR(SUBR(mvec_set_block)(Vv,Qv,m0,ncVQ-1,ierr),*ierr);
        }
      PHIST_CHK_NEG_IERR(SUBR(orthog)(Vv,t_ptr,Tv,Sv,3,ierr),*ierr);
      // reset the view of V without the Q.
      PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&Vv,0,m0-1,ierr),*ierr);
      
      // set V(:,m)=t
      PHIST_CHK_IERR(SUBR(mvec_set_block)(V,t,m0,m-1,ierr),*ierr);
      // compute AV(:,m) = A*V(:,m)
      PHIST_CHK_IERR(A_op->apply(st::one(),A_op->A,Vm,st::zero(),AVm,ierr),*ierr);
      // Galerkin for non-Hermitian A
      // TODO - is it maybe more efficient to just compute V'AV completely? It
      //        means more ops and more data transfer but fewer messages.
      // M(1:m-1,m)=V(:,1:m-1)'*AV(:,m);
      PHIST_CHK_IERR(SUBR(sdMat_view_block)(M,&Mv,0,m0-1,m0,m-1,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),Vv,AVm,st::zero(),Mv,ierr),*ierr);
      // M(m,1:m-1)=V(:,m)'*AV(:,1:m-1);
      PHIST_CHK_IERR(SUBR(sdMat_view_block)(M,&Mv,m0,m-1,0,m0-1,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),Vm,AVv,st::zero(),Mv,ierr),*ierr);
      // M(m,m)=V(:,m)'*AV(:,m) note that we could do this using a dot product in the
      // single-vector case, but we want to extend it to a block variant so we use the more
      // general dense matmul.
      PHIST_CHK_IERR(SUBR(sdMat_view_block)(M,&Mv,m0,m-1,m0,m-1,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),Vm,AVm,st::zero(),Mv,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&Vv,0,m-1,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(mvec_view_block)(AV,&AVv,0,m-1,ierr),*ierr);
      PHIST_DEB("basis size is now m=%d",m+1);
#ifdef TESTING
      // check orthogonality of [V Q]
      sdMat_ptr_t tmp1=NULL,tmp2=NULL;
      ST *tmp1_raw, *tmp2_raw;
      int ld1,ld2;
      PHIST_CHK_IERR(SUBR(sdMat_create)(&tmp1,m,m,NULL,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(sdMat_extract_view)(tmp1,&tmp1_raw,&ld1,ierr),*ierr);
        int ncV;
        PHIST_CHK_IERR(SUBR(mvec_num_vectors)(Vv,&ncV,ierr),*ierr);
        PHIST_OUT(PHIST_VERBOSE,"#vectors in V: %d",ncV);
      PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),Vv,Vv,st::zero(),tmp1,ierr),*ierr);
      if (nconv>0)
        {
        int ncQ;
        PHIST_CHK_IERR(SUBR(mvec_num_vectors)(Qv,&ncQ,ierr),*ierr);
        PHIST_OUT(PHIST_VERBOSE,"#vectors in Q: %d",ncQ);
        PHIST_CHK_IERR(SUBR(sdMat_create)(&tmp2,m,nconv,NULL,ierr),*ierr);
        PHIST_CHK_IERR(SUBR(sdMat_extract_view)(tmp2,&tmp2_raw,&ld2,ierr),*ierr);
        PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),Vv,Qv,st::zero(),tmp2,ierr),*ierr);
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
      PHIST_OUT(PHIST_VERBOSE,"orthogonality monitors:");
      PHIST_OUT(PHIST_VERBOSE,"||V'V-I||=%e",err1);
      PHIST_OUT(PHIST_VERBOSE,"||V'Q||=%e",err2);
      PHIST_CHK_IERR(SUBR(sdMat_delete)(tmp1,ierr),*ierr);
      if (nconv>0)
        {
        PHIST_CHK_IERR(SUBR(sdMat_delete)(tmp2,ierr),*ierr);
        }
#endif
      }
    else
      {
      expand=1;//TODO - handle nv correctly
      }
      
    // now do a Schur decomposition of M into S and T, with the
    // Ritz values on the diagonal of T.

    // copy T=M. Note that as we just 'view' the upper left blocks,
    // the pointers M_raw, T_raw, S_raw are still correct
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(M,&Mv,0,m-1,0,m-1,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(T,&Tv,0,m-1,0,m-1,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_get_block)(Mv,Tv,0,m-1,0,m-1,ierr),*ierr);

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
    PHIST_CHK_IERR(SUBR(SchurDecomp)
        (T_raw,ldT,S_raw,ldS,m,nselect,nsort,which,ev,ierr),*ierr);
    int ev_pos=0;
    theta=ev[ev_pos];

#if PHIST_OUTLEV>=PHIST_DEBUG
  PHIST_DEB("it=%d, m=%d, mm=%d",it,m,mm);
  PHIST_DEB("projected matrix");
  PHIST_CHK_IERR(SUBR(sdMat_print)(Mv,ierr),*ierr);
  PHIST_DEB("sorted Schur form");
  PHIST_CHK_IERR(SUBR(sdMat_print)(Tv,ierr),*ierr);
#endif

    nv=1; // nv==2 indicates complex eigenvalue of real matrix, we
          // will use this flag later on.
    // set some pointers
    u_ptr=u;
    Au_ptr=Au;
    r_ptr=r;
    rtil_ptr=rtil;
    t_ptr=t;
#ifndef _IS_COMPLEX_
    if (mt::abs(ct::imag(theta))>mt::eps())
      {
      PHIST_DEB("complex eigenvalue in real arithmetic");
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
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(T,&Theta,0,nv-1,0,nv-1,ierr),*ierr);

    // get Ritz vector corresponding to theta
    // (for complex RV in real case, two columns are extracted)
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(S,&Sv,0,m-1,0,nv-1,ierr),*ierr);
    //u=V*s;
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(st::one(),Vv,Sv,st::zero(),u_ptr,ierr),*ierr);
    //Au=AV*s;
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(st::one(),AVv,Sv,st::zero(),Au_ptr,ierr),*ierr);

// part CMP_RESID: compute residual and its norm (TODO: make this a helper function)
//{
    // r=Au-theta*u; ([r_r, r_i] = [Au_r, Au_i] - [u_r, u_i]*Theta in the real case with 
    // complex theta)
    
    // first: r=Au
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),Au_ptr,st::zero(),r_ptr,ierr),*ierr);

    // update r = r - U*Theta
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(-st::one(),u_ptr,Theta,st::one(),r_ptr,ierr),*ierr);

    // set rtil=r
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),r_ptr,st::zero(),rtil_ptr,ierr),*ierr);

    // project out already converged eigenvectors
    // TODO - we could use our orthog routine here instead
    if (nconv>0)
      {
      //atil = Q'*r;
      PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),Qv,r_ptr,st::zero(),atilv,ierr),*ierr);
    
      //rtil = r-Q*atil;
      PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(-st::one(),Qv,atilv,st::one(),rtil_ptr,ierr),*ierr);
      }

    //nrm=norm(rtil);
    // real case with complex r: ||v+iw||=sqrt((v+iw).'*(v-iw))=sqrt(v'v+w'w)
    nrm[1]=mt::zero();
    // in the complex case we pass in a 're-interpret cast' of nrm as complex,
    // which should be fine (imaginary part will be 0)
    PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(rtil_ptr,rtil_ptr,(ST*)nrm,ierr),*ierr);
    nrm[0]=mt::sqrt(nrm[0]*nrm[0]+nrm[1]*nrm[1]);
//}
    PHIST_OUT(PHIST_INFO,"JDQR Iter %d\tdim(V)=%d\ttheta=%8.4g%+8.4gi\t\tr_est=%8.4e\n",it,m,ct::real(theta),ct::imag(theta),nrm[0]);
  // deflate converged eigenpairs
  while (nrm[0]<=tol)
    {
    mm=0;// counts iterations for next theta
    // number of converged eigenpairs: +2 for complex conjugate pairs in real case
    int nq0=nconv;
    nconv=nconv+nv;
    for (int j=nq0;j<nconv;j++)
      {
      resid[j]=nrm[0]; // will be returned to the user
      }

    // view first nconv columns of X as 'Q'
    PHIST_CHK_IERR(SUBR(mvec_view_block)(X,&Qv,0,nconv-1,ierr),*ierr);

     //Q=[Q,u];
     PHIST_CHK_IERR(SUBR(mvec_set_block)(Qv,u_ptr,nq0,nconv-1,ierr),*ierr);
     //R=[R, atil;
     //   0, theta];
     if (nq0>0)
       {
       PHIST_CHK_IERR(SUBR(sdMat_set_block)(R,atilv,0,nq0-1,nq0,nconv-1,ierr),*ierr);
       }
     PHIST_CHK_IERR(SUBR(sdMat_set_block)(R,Theta,nq0,nconv-1,nq0,nconv-1,ierr),*ierr);

     // view next ~a, a temporary vector to compute ~a=Q'*r
     PHIST_CHK_IERR(SUBR(sdMat_view_block)(atil,&atilv,0,nconv-1,0,nv-1,ierr),*ierr);

    PHIST_OUT(PHIST_VERBOSE,"eigenvalue %d (%8.4g%+8.4gi) converged.",nconv,ct::real(theta),ct::imag(theta));
    if (nconv>=numEigs)
      {
      PHIST_OUT(PHIST_VERBOSE,"stopping JDQR loop");
      solve=false;
      break;
      }

    //S=S(:,2:m); (select remaining Ritz vectors)
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(S,&Sv,0,m-1,nv,m,ierr),*ierr);
    //M=T(2:m,2:m);
    // let Mblock point to the first m-nv x m-nv block of M
      PHIST_CHK_IERR(SUBR(sdMat_view_block)(M,&Mv,0,m-nv-1,0,m-nv-1,ierr),*ierr);
    // copy T(nv+1:m,nv+1:m) into it
    PHIST_CHK_IERR(SUBR(sdMat_get_block)(T,Mv,nv,m-1,nv,m-1,ierr),*ierr);
    int m0=m;
    m=m-nv;
    it=it-1;
    //V=V*S;
    // It is not allowed to alias V in mvec_times_sdMat, so we use a temporary vector
    mvec_ptr_t v_tmp=NULL;
    PHIST_CHK_IERR(SUBR(mvec_view_block)(Vtmp,&v_tmp,0,m0-1,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(st::one(),Vv,Sv,st::zero(),v_tmp,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_set_block)(Vv,v_tmp,0,m0-1,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&Vv,0,m-1,ierr),*ierr);
    //AV=AV*S;
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(st::one(),AVv,Sv,st::zero(),v_tmp,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_set_block)(AVv,v_tmp,0,m0-1,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_view_block)(AV,&AVv,0,m-1,ierr),*ierr);

    //u=V(:,1);
    PHIST_CHK_IERR(SUBR(mvec_get_block)(Vv,u_ptr,0,nv-1,ierr),*ierr);
    //Au=AV(:,1);
    PHIST_CHK_IERR(SUBR(mvec_get_block)(AVv,Au_ptr,0,nv-1,ierr),*ierr);

    // select next target ev. This need not be the next according to
    // the sort criterion because the Schur-form is not completely sorted
    // (TODO: is it desirable to do them in order?)

    ev_pos+=nv; // skip conjugate ev for real case w/ cmplx theta
                // TODO: we would have to add it ot Q anyway, wouldn't we?
    theta=ev[ev_pos];
    nv=1;
#ifndef _IS_COMPLEX_
    if (mt::abs(ct::imag(theta))>mt::eps())
      {
      nv++;
      }
#endif
    // get the diagonal block (1x1 or 2x2) corresponding to theta
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(T,&Theta,ev_pos,ev_pos+nv-1,ev_pos,ev_pos+nv-1,ierr),*ierr);

//TODO: same as part CMP_RESID up there, put it in an aux function
//{
    // r=Au-theta*u; ([r_r, r_i] = [Au_r, Au_i] - [u_r, u_i]*Theta in the real case with 
    // complex theta)
    
    // first: r=Au
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),Au_ptr,st::zero(),r_ptr,ierr),*ierr);

    // update r = r - U*Theta
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(-st::one(),u_ptr,Theta,st::one(),r_ptr,ierr),*ierr);

    // set rtil=r
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),r_ptr,st::zero(),rtil_ptr,ierr),*ierr);

    //atil = Q'*r;
    PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),Qv,r_ptr,st::zero(),atilv,ierr),*ierr);

    //rtil = r-Q*atil;
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(-st::one(),Qv,atilv,st::one(),r,ierr),*ierr);
      
    //nrm=norm(rtil);
    // real case with complex r: ||v+iw||=sqrt((v+iw).'*(v-iw))=sqrt(v'v+w'w)
    nrm[1]=mt::zero();
    PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(rtil,rtil_ptr,(ST*)nrm,ierr),*ierr);
    nrm[0]=mt::sqrt(nrm[0]*nrm[0]+nrm[1]*nrm[1]);
//}
    // again, sort the largest Ritz value to the top
//    nselect = m==maxBas-1? minBas: 1;
//    nsort = 1;
    nselect=m;
    nsort=nselect;
    // set T=M
    PHIST_CHK_IERR(SUBR(sdMat_set_block)(T,Mv,0,m-1,0,m-1,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(SchurDecomp)
        (T_raw,ldT,S_raw,ldS,m,nselect,nsort,which,ev,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(T,&Tv,0,m-1,0,m-1,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(S,&Sv,0,m-1,0,m-1,ierr),*ierr);
    } //while (deflate)

  // restart if necessary
  if (m>=maxBas)
    {
    PHIST_OUT(PHIST_VERBOSE,"restart JDQR");
    int m0=m;
    m=minBas;
    
    //S=S(:,1:m);
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(S,&Sv,0,m0-1,0,m-1,ierr),*ierr);
    //M=T(1:m,1:m);
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(M,&Mv,0,m-1,0,m-1,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_get_block)(T,Mv,0,m-1,0,m-1,ierr),*ierr);
    //V=V*S;
    mvec_ptr_t v_tmp=NULL;
    PHIST_CHK_IERR(SUBR(mvec_view_block)(Vtmp,&v_tmp,0,m-1,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(st::one(),Vv,Sv,st::zero(),v_tmp,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_set_block)(Vv,v_tmp,0,m-1,ierr),*ierr);
    //AV=AV*S;
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(st::one(),AVv,Sv,st::zero(),v_tmp,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_set_block)(AVv,v_tmp,0,m-1,ierr),*ierr);

    //V=V(:,1:mmin);
    PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&Vv,0,m-1,ierr),*ierr);
    //AV=AV(:,1:mmin);
    PHIST_CHK_IERR(SUBR(mvec_view_block)(AV,&AVv,0,m-1,ierr),*ierr);
    solve=true;
    expand=nv;
    }
  if (solve)
    {
    if (nv>1)
      {
      PHIST_OUT(PHIST_ERROR,"case real A with complex eigs not (fully) implemented");
      PHIST_CHK_IERR(*ierr=-99,*ierr);
      }
    // maintain orthogonality against
    // all converged vectors (Q) and the
    // new one u:
    //Qtil=[Q,u];
    PHIST_CHK_IERR(SUBR(mvec_view_block)(X,&Qtil,0,nconv+nv-1,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_set_block)(Qtil,u_ptr,nconv,nconv+nv-1,ierr),*ierr);
    /*
    if (nrm[0]<switchTol)
      {
      PHIST_OUT(1,"using RQI");
      shift=theta;
      }
    else
      {
      PHIST_OUT(1,"using SI");
      shift=target;
      }
    */
    
#ifdef _IS_COMPLEX_
//    shift=theta;
    PHIST_CHK_IERR(SUBR(sdMat_put_value)(shift,theta,ierr),*ierr);
#else
  if (nv>1)
    {
    PHIST_OUT(PHIST_ERROR,"real case with complex eig not implemented (file %s, line %d)",__FILE__,__LINE__);
    PHIST_CHK_IERR(*ierr=-99,*ierr);
    }
  PHIST_CHK_IERR(SUBR(sdMat_put_value)(shift,ct::real(theta),ierr),*ierr);
#endif    
    // solve approximately 
    // (I-uu')(A-theta*I)(I-uu')*t=-r
    // to get t \orth u (u ^= Qtil here)

    // the block case is not implemented yet:
    // (I-uu')A(I-uu')*t - (I-uu')t*Theta=-r
    op_ptr_t jada_op=NULL;
    if (B_op!=NULL)
      {
      PHIST_OUT(PHIST_WARNING,"case B!=I not implemented (file %s, line %d)",__FILE__,__LINE__);
      }
      
    PHIST_CHK_IERR(SUBR(jadaOp_create)(A_op,NULL,Qtil,NULL,shift,NULL,&jada_op,ierr),*ierr);

    // 1/2^mm, but at most the outer tol as conv tol for GMRES
    MT innerTol = std::max(tol,mt::one()/((MT)(2<<mm)));
    PHIST_OUT(PHIST_VERBOSE,"inner conv tol: %g",innerTol);

    // allow at most 25 iterations (TODO: make these settings available to the user)
    int nIt=25;
    int maxKSpace=25; // maximum size of Krylov subspace created
    
    //TODO - scaling and preconditioning
    
    // set t=0 as initial guess (TODO, better start vector)
    PHIST_CHK_IERR(SUBR(mvec_put_value)(t_ptr,st::zero(),ierr),*ierr);

    int variant=0; //0:block GMRES, 1: pseudo-BGMRES
    SUBR(bgmres)(jada_op,t_ptr,rtil_ptr,innerTol,&nIt,maxKSpace,variant,NULL,ierr);
      
    expand=true;
    }
  else
    {
    expand=0;
    solve=true;
    }// if solve
  }// while loop
    
  *num_iters=it;
  
  // TODO Sleijpen in his jdqr code does a refinement step at this      
  // point taking the information from the current V into Q and R       
  // (might be a nice feature). It also allows us to return unconverged 
  // eigenvectors so the user can restart from them.                    
  // In that case the value of the residual norm in resid should be     
  // considered an upper bound (or updated, of course).                 
  TOUCH(Rv); // we would need it for this
  
  // some sanity checks - the user provides only num_eigs+1 slots
  nconv = std::min(*num_eigs+1,nconv);
  *num_eigs=nconv;// tell the user how many we return
  
  // copy at most num_eigs+1 converged eigenvalues into the user
  // provided array
  ST* R_raw=NULL;
  lidx_t ldR;
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(R,&R_raw,&ldR,ierr),*ierr);

  i=0;
  while (i<nconv)
    {
    evals[i] = R_raw[i*ldR+i];
#ifndef _IS_COMPLEX_
    if (i<nconv-1)
      {
      if (st::abs(R_raw[i*ldR+i+1])>mt::eps())
        {
        is_cmplx[i]=1;
        // sqrt((T^2)/4-D) gives the imaginary part of the eigenpair of
        // a 2x2 matrix, with T the trace and D the determinant. In our
        // case A11=A22 and the thing simplifies to im(lambda)=A12*A21.
        evals[i+1]=T_raw[i*ldR+i+1]*T_raw[(i+1)*ldR+i];
        is_cmplx[i+1]=1;
        i++;
        }
      }
#else
    TOUCH(is_cmplx);
#endif    
    i++;
    }
  
  // TODO - compute eigenvectors 
  // * Jordan decomposition: R*S = S*J, with J in Jordan form (diagonal with possibly a 1 on 
  // the super diagonal and the eigenvalues on the diagonal)
  // * X = Q*S
  // (which LAPACK routines? XTREVC to get eigenvectors?)
  PHIST_OUT(PHIST_WARNING,"warning: returning only a basis of the computed eigenspace");

  // free memory
  PHIST_CHK_IERR(SUBR(mvec_delete)(V,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_delete)(AV,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_delete)(Vtmp,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_delete)(u,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_delete)(Au,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_delete)(r,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_delete)(rtil,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_delete)(t,ierr),*ierr);
  
  PHIST_CHK_IERR(SUBR(sdMat_delete)(M,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(S,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(T,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(R,ierr),*ierr);
  
  // TODO - we should also call delete for the views to avoid small memory leaks
  
  return;
  }


