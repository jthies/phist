
//! a simple Arnoldi process to start up the JaDa iteration.
//! Given a minimum basis size m, compute V(:,1:m), 
//! H(1:m+1,1:m) such that A*V(:,1:m) = V(:,1:m+1)*H(1:m+1,1:m)
//! input: v0 with ||v0||_2=1, V and H allocated with m+1 columns
//! and nloc resp. m+1 rows.
//!
//! TODO - block Arnoldi, cf. matlab/krylov/arnoldi.m for a prototype.
//! TODO - we may want to include a check if any Ritz values have already
//!        converged in Arnoldi, but this requires some additional programming
//!        effort which we leave for later.
void SUBR(arnoldi)(TYPE(const_op_ptr) op, TYPE(const_mvec_ptr) v0,
        TYPE(mvec_ptr) V, TYPE(sdMat_ptr) H, int m, int* ierr);

 //                                                                                             
 // this function does an in-place Schur decomposition of T(1:m,1:m) into T and S.              
 // The Ritz values appear on the diagonal of T (for the real case there may be 2x2 blocks      
 // for complex conjugate pairs). The nselect and nsort flags indicate in which order they      
 // should appear:                                                                              
 // nselect=nsort=0: unsorted                                                                   
 // nselect>0: the first <nselect> Ritz values (if the last one is a complex conjugate pair     
 //            <nselect+1>) appear in any order in the upper left corner of T                    
 // 0<nsort<nselect: in addition to moving a cluster of <nselect> eigenvalues, sort the         
 //             first <nsort> of them in the upper left corner.                                 
 //                                                                                             
 // Example: if the Schur form is                                                               
 //                                                                                             
 //     a x x x x                                                                               
 //     0 b c x x, and the eigenvalues are |a|<|e|<|lambda([b,c;d,b])|<|f|,                     
 //     0 d b x x                                                                               
 //     0 0 0 e x                                                                               
 //     0 0 0 0 f                                                                               
 //                                                                                             
 // we get for nselect=3, nsort=1, which=LM:                                                    
 //                                                                                             
 //   f x x x x                                                                                 
 //     b c x x                                                                                 
 //     d b x x                                                                                 
 //         e x                                                                                 
 //           a                                                                                 
 //                                                                                             
 // so we guarantee that the largest one <nsort> is in the upper left corner,                   
 // and that the 3 largest ones appear first in any order.                                      
 //                                                                                             
 void SUBR(SchurDecomp)(_ST_* T, int ldT, _ST_* S, int ldS,
         int m, int nselect, int nsort, eigSort_t which, 
         std::complex<_MT_>* ev, int *ierr);

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
  // the above vectors may be complex even for a real
  // matrix. These are real- and imaginary parts, respectively
  // (in the complex case x_r aliases x and x_i==NULL)
  mvec_ptr_t    u_r=NULL,       u_i=NULL,
                Au_r=NULL,      Au_i=NULL,
                r_r=NULL,       r_i=NULL,
                rtil_r=NULL,    rtil_i=NULL,
                t_r=NULL,       t_i=NULL;
  // these point either to u or u_r etc. to make live simpler further down
  mvec_ptr_t u_ptr=NULL, Au_ptr=NULL, r_ptr=NULL, rtil_ptr=NULL,t_ptr=NULL;
  // Q*s (temporary vector)
  sdMat_ptr_t atil=NULL, atil_r=NULL, atil_i=NULL, atil_ptr;
  // matrix <V,AV>
  sdMat_ptr_t M=NULL;
  sdMat_ptr_t Mv=NULL; // to create views of parts of M
  // Schur-decomposition of M
  sdMat_ptr_t T, S=NULL;
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
  ST shift;
  // theta as a sdMat_t (either 1x1 or 2x2 in real case)
  sdMat_ptr_t Theta=NULL;
  std::complex<MT> ev[maxBas];
  
  int numEigs=*num_eigs;
  int maxIter=*num_iters;
  int i,it,m,mm;
  MT nrm[2]; // for computing residual norms

  *num_iters=0;
  
  const_comm_ptr_t comm;
  PHIST_CHK_IERR(phist_map_get_comm(A_op->range_map,&comm,ierr),*ierr);

  PHIST_CHK_IERR(SUBR(mvec_create)(&V,A_op->domain_map,maxBas,ierr),*ierr);
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
  
  u_r=u; u_i=NULL;
  Au_r=Au; Au_i=NULL;
  r_r=r; r_i=NULL;
  rtil_r=rtil; rtil_i=NULL;
  t_r=t; t_i=NULL;

#ifndef _IS_COMPLEX_
  PHIST_CHK_IERR(SUBR(mvec_view_block)(u,&u_r,0,0,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_view_block)(u,&u_i,1,1,ierr),*ierr);

  PHIST_CHK_IERR(SUBR(mvec_view_block)(Au,&Au_r,0,0,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_view_block)(Au,&Au_i,1,1,ierr),*ierr);

  PHIST_CHK_IERR(SUBR(mvec_view_block)(r,&r_r,0,0,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_view_block)(r,&r_i,1,1,ierr),*ierr);

  PHIST_CHK_IERR(SUBR(mvec_view_block)(rtil,&rtil_r,0,0,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_view_block)(rtil,&rtil_i,1,1,ierr),*ierr);  

  PHIST_CHK_IERR(SUBR(mvec_view_block)(t,&t_r,0,0,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_view_block)(t,&t_i,1,1,ierr),*ierr);  
#endif
  PHIST_CHK_IERR(SUBR(sdMat_create)(&M,maxBas,maxBas,comm,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_create)(&S,maxBas,maxBas,comm,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_create)(&T,maxBas,maxBas,comm,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_create)(&R,*num_eigs,*num_eigs,comm,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_put_value)(R,st::zero(),ierr),*ierr);

  // pointer to the data in M, S and T
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(M,&M_raw, &ldM,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(S,&S_raw, &ldS,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(T,&T_raw, &ldT,ierr),*ierr);

  // view the first minBas+1 columns of V and H as M(1:minBas+1,1:minBas)
  PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&Vv,0,minBas,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_view_block)(M,&H0, 0,minBas,0,minBas-1,ierr),*ierr);

  PHIST_OUT(1,"%d steps of Arnoldi as start-up",minBas);
  
  // start by filling the first minBas vectors using Arnoldi
  //TODO - we promised to use the user input in X, here we just
  // use a random vector instead.
  mvec_ptr_t v0=NULL;
  PHIST_CHK_IERR(SUBR(mvec_view_block)(Vtmp,&v0,0,0,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_random)(v0,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(arnoldi)(A_op,v0,Vv,H0,minBas,ierr),*ierr);

#if PHIST_OUTLEV>PHIST_DEBUG
  PHIST_DEB("initial H (by Arnoldi)");
  PHIST_CHK_IERR(SUBR(sdMat_print)(H0,ierr),*ierr);
#endif
  // delete the view (not the data) v0
  PHIST_CHK_IERR(SUBR(mvec_delete)(v0,ierr),*ierr);
  v0=NULL; // always important to nullify pointers in phist...
  PHIST_CHK_IERR(SUBR(mvec_view_block)(AV,&AVv,0,minBas-1,ierr),*ierr);
  // compute A*V for the first minBas columns which we have now using the Arnoldi relation
  PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(st::one(),Vv,H0,st::zero(),AVv,ierr),*ierr);

  // set Vv=V(:,1:minBas)
  PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&Vv,0,minBas-1,ierr),*ierr);

  PHIST_OUT(1,"Jacobi-Davidson");
  PHIST_OUT(1,"%s\t%s\t%s\t\t%s\n","iter","m","approx","resid");

  it=0; // total number of JaDa iterations
  mm=0; //count iterations per eigenvalue/pair of complex ev's
  m=minBas-1;// location of last valid vector in V
  expand=0;
  int nv=1; // number of vectors in current update (1 or 2 in this implementation, 2 for
            // complex eigenvectors of a real matrix.
  
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
      PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&Vv,0,m0,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(mvec_view_block)(AV,&AVv,0,m0,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&Vm,m0+1,m,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(mvec_view_block)(AV,&AVm,m0+1,m,ierr),*ierr);
      // orthogonalize t against V(:,0:m-1). We use T and S as temporary storage here
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(T,&Tv,0,nv-1,0,nv-1,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(S,&Sv,0,m0,0,nv-1,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(orthog)(Vv,t_ptr,Tv,Sv,2,ierr),*ierr);
      // set V(:,m)=t
      PHIST_CHK_IERR(SUBR(mvec_set_block)(V,t,m0+1,m,ierr),*ierr);
      // compute AV(:,m) = A*V(:,m)
      PHIST_CHK_IERR(A_op->apply(st::one(),A_op->A,Vm,st::zero(),AVm,ierr),*ierr);
      // Galerkin for non-Hermitian A
      // TODO - is it maybe more efficient to just compute V'AV completely? It
      //        means more ops and more data transfer but fewer messages.
      // M(1:m-1,m)=V(:,1:m-1)'*AV(:,m);
      PHIST_CHK_IERR(SUBR(sdMat_view_block)(M,&Mv,0,m0,m0+1,m,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),Vv,AVm,st::zero(),Mv,ierr),*ierr);
      // M(m,1:m-1)=V(:,m)'*AV(:,1:m-1);
      PHIST_CHK_IERR(SUBR(sdMat_view_block)(M,&Mv,m0+1,m,0,m0,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),Vm,AVv,st::zero(),Mv,ierr),*ierr);
      // M(m,m)=V(:,m)'*AV(:,m) note that we could do this using a dot product in the
      // single-vector case, but we want to extend it to a block variant so we use the more
      // general dense matmul.
      PHIST_CHK_IERR(SUBR(sdMat_view_block)(M,&Mv,m0+1,m,m0+1,m,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),Vm,AVm,st::zero(),Mv,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&Vv,0,m,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(mvec_view_block)(AV,&AVv,0,m,ierr),*ierr);
      PHIST_DEB("basis size is now m=%d",m+1);
      }
    else
      {
      expand=1;//TODO - handle nv correctly
      }
      
    // now do a Schur decomposition of M into S and T, with the
    // Ritz values on the diagonal of T.

    // copy T=M. Note that as we just 'view' the upper left blocks,
    // the pointers M_raw, T_raw, S_raw are still correct
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(M,&Mv,0,m,0,m,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(T,&Tv,0,m,0,m,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_get_block)(Mv,Tv,0,m,0,m,ierr),*ierr);

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

#if PHIST_OUTLEV>PHIST_DEBUG
  PHIST_DEB("it=%d, m=%d, mm=%d",it,m,mm);
  PHIST_DEB("sorted Schur form");
  PHIST_CHK_IERR(SUBR(sdMat_print)(Tv,ierr),*ierr);
#endif

    nv=1; // nv==2 indicates complex eigenvalue of real matrix, we
          // will use this flag later on.
    // set some pointers
    u_ptr=u;
    Au_ptr=Au;
    atil_ptr=atil;
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
      atil_ptr=atil_r;
      u_ptr=u_r;
      Au_ptr=Au_r;
      r_ptr=r_r;
      rtil_ptr=rtil_r;
      }
#endif

    // get the diagonal block (1x1 or 2x2) corresponding to theta
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(T,&Theta,0,nv-1,0,nv-1,ierr),*ierr);

    // get Ritz vector corresponding to theta
    // (for complex RV in real case, two columns are extracted)
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(S,&Sv,0,m,0,nv-1,ierr),*ierr);
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
    if (nconv>0)
      {
      //atil = Q'*r;
      PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),Qv,r,st::zero(),atil,ierr),*ierr);
    
      //rtil = r-Q*atil;
      PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(-st::one(),Qv,atil,st::one(),r,ierr),*ierr);
      }

    //nrm=norm(rtil);
    // real case with complex r: ||v+iw||=sqrt((v+iw).'*(v-iw))=sqrt(v'v+w'w)
    nrm[1]=mt::zero();
    // in the complex case we pass in a 're-interpret cast' of nrm as complex,
    // which should be fine (imaginary part will be 0)
    PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(rtil_ptr,rtil_ptr,(ST*)nrm,ierr),*ierr);
    nrm[0]=mt::sqrt(nrm[0]*nrm[0]+nrm[1]*nrm[1]);
//}
    PHIST_OUT(1,"%d\t%d\t%8.4g%+8.4gi\t\t%8.4g\n",it,m,ct::real(theta),ct::imag(theta),nrm[0]);
  // deflate converged eigenpairs
  while (nrm[0]<=tol)
    {
    mm=0;// counts iterations for next theta
    // number of converged eigenpairs: +2 for complex conjugate pairs in real case
    int nq0=nconv;
    nconv=nconv+nv;

    // view first nconv columns of X as 'Q'
    PHIST_CHK_IERR(SUBR(mvec_view_block)(X,&Qv,0,nconv-1,ierr),*ierr);

     //Q=[Q,u];
     PHIST_CHK_IERR(SUBR(mvec_set_block)(Qv,u_ptr,nq0,nconv-1,ierr),*ierr);
     //R=[R, atil;
     //   0, theta];
     PHIST_CHK_IERR(SUBR(sdMat_set_block)(R,atil_ptr,0,m,nq0,nconv-1,ierr),*ierr);
     PHIST_CHK_IERR(SUBR(sdMat_set_block)(R,Theta,nq0,nconv-1,nq0,nconv-1,ierr),*ierr);

    PHIST_OUT(1,"eigenvalue %d (%8.4g%+8.4gi) is converged.",nconv,ct::real(theta),ct::imag(theta));
    if (nconv>=numEigs)
      {
      solve=false;
      break;
      }

    //S=S(:,2:m); (select remaining Ritz vectors)
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(S,&Sv,0,m-1,nv,m,ierr),*ierr);
    //M=T(2:m,2:m);
    // let Mblock point to the first m-nv x m-nv block of M
      PHIST_CHK_IERR(SUBR(sdMat_view_block)(M,&Mv,0,m-nv-1,0,m-nv-1,ierr),*ierr);
    // copy T(nv+1:m,nv+1:m) into it
    PHIST_CHK_IERR(SUBR(sdMat_get_block)(T,Mv,nv,m,nv,m,ierr),*ierr);
    m=m-nv;
    //V=V*S;
    // It is not allowed to alias V in mvec_times_sdMat, so we use a temporary vector
    mvec_ptr_t v_tmp=NULL;
    PHIST_CHK_IERR(SUBR(mvec_view_block)(Vtmp,&v_tmp,0,m,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(st::one(),Vv,Sv,st::zero(),v_tmp,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_set_block)(Vv,v_tmp,0,m,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&Vv,0,m,ierr),*ierr);
    //AV=AV*S;
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(st::one(),AVv,Sv,st::zero(),v_tmp,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_set_block)(AVv,v_tmp,0,m,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_view_block)(AV,&AVv,0,m,ierr),*ierr);

    //u=V(:,1);
    PHIST_CHK_IERR(SUBR(mvec_get_block)(Vv,u,0,nv,ierr),*ierr);
    //Au=AV(:,1);
    PHIST_CHK_IERR(SUBR(mvec_get_block)(AVv,Au,0,nv,ierr),*ierr);

    // select next target ev. This need not be the next according to
    // the sort criterion because the Schur-form is not completely sorted
    // (TODO: is it desirable to do them in order?)

    ev_pos+=nv; // skip conjugate ev for real case w/ cmplx theta
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
    PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),Qv,r_ptr,st::zero(),atil,ierr),*ierr);
    
    //rtil = r-Q*atil;
      PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(-st::one(),Qv,atil_ptr,st::one(),r,ierr),*ierr);

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
    PHIST_CHK_IERR(SUBR(sdMat_set_block)(T,Mv,0,m,0,m,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(SchurDecomp)
        (T_raw,ldT,S_raw,ldS,m,nselect,nsort,which,ev,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(T,&Tv,0,m,0,m,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(S,&Sv,0,m,0,m,ierr),*ierr);
    } //while (deflate)

  // restart if necessary
  if (m>=maxBas)
    {
    PHIST_OUT(1,"restart JDQR");
    int m0=m;
    m=minBas-1;
    
    //S=S(:,1:m);
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(S,&Sv,0,m0-1,0,m-1,ierr),*ierr);
    //M=T(1:m,1:m);
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(M,&Mv,0,m-1,0,m-1,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_get_block)(T,Mv,0,m-1,0,m-1,ierr),*ierr);
    //V=V*S;
    mvec_ptr_t v_tmp=NULL;
    PHIST_CHK_IERR(SUBR(mvec_view_block)(Vtmp,&v_tmp,0,m,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(st::one(),Vv,Sv,st::zero(),v_tmp,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_set_block)(Vv,v_tmp,0,m,ierr),*ierr);
    //AV=AV*S;
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(st::one(),AVv,Sv,st::zero(),v_tmp,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_set_block)(AVv,v_tmp,0,m,ierr),*ierr);

    //V=V(:,1:mmin);
    PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&Vv,0,m,ierr),*ierr);
    //AV=AV(:,1:mmin);
    PHIST_CHK_IERR(SUBR(mvec_view_block)(AV,&AVv,0,m,ierr),*ierr);
    solve=true;
    expand=nv;
    }
  if (solve)
    {
    if (nv>1)
      {
      PHIST_OUT(PHIST_ERROR,"case real A with complex eigs not (fully) implemented");
      PHIST_CHK_IERR(-99,*ierr);
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
    shift=theta;
#else
    shift=ct::real(theta);
#endif    
    // solve approximately 
    // (I-uu')(A-theta*I)(I-uu')*t=-r
    // to get t \orth u (u ^= Qtil here)

    // the block case is not implemented yet:
    // (I-uu')A(I-uu')*t - (I-uu')t*Theta=-r
    op_ptr_t jada_op;
    PHIST_CHK_IERR(SUBR(jadaOp_create)(A_op,B_op,Qtil,shift,&jada_op,ierr),*ierr);
    
    // 1/2^mm, but at most the outer tol as conv tol for GMRES
    MT innerTol = std::max(tol,mt::one()/((MT)(2<<mm)));
    PHIST_OUT(1,"inner conv tol: %g",innerTol);

    // allow at most 25 iterations (TODO: make these settings available to the user)
    int nIt=25;
    int maxKSpace=25; // maximum size of Krylov subspace created
    
    //TODO - scaling and preconditioning
    
    // set t=-r as initial guess
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),rtil,st::zero(),u,ierr),*ierr);

    SUBR(bgmres)(jada_op,t_ptr,rtil_ptr,innerTol,&nIt,maxKSpace,ierr);
      
    expand=true;

//    PHIST_OUT(1,"(Qtil,t)=%f",num2str(norm(Qtil'*t)));
//    PHIST_OUT(1,"(A-sI)t+rtil=%f",norm(A*t-shift*t+rtil));         
    }
  else
    {
    expand=0;
    solve=true;
    }// if solve
  }// while loop
    
  *num_iters=it;
  // copy at most num_eigs+1 converged eigenvalues into the user
  // provided array
  /* TODO
  for (i=n-1;i>=0, inext<*num_eigs;i--)
    {
    if (converged[i])
      {
      evals[inext] = falphas[i];
      resid[inext] = r_est[i];
      inext++;
      }
    }
  */
  
  if (nconv!=*num_eigs) *num_eigs=nconv;
  
  // TODO - compute eigenvectors

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



void SUBR(arnoldi)(TYPE(const_op_ptr) op, TYPE(const_mvec_ptr) v0,
        TYPE(mvec_ptr) V, TYPE(sdMat_ptr) H, int m, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  TYPE(mvec_ptr) v=NULL,av=NULL,vprev=NULL;
  TYPE(sdMat_ptr) R1=NULL,R2=NULL;
  PHIST_CHK_IERR(SUBR(mvec_set_block)(V,v0,0,0,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&v,0,0,ierr),*ierr);
  
  for (int i=0;i<m;i++)
    {
    PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&vprev,0,i,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&av,i+1,i+1,ierr),*ierr);
    PHIST_CHK_IERR(op->apply(st::one(),op->A,v,st::zero(),av,ierr),*ierr);
    // orthogonalize, Q*R1 = W - V*R2
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(H,&R2,0,i,i,i,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(H,&R1,i+1,i+1,i,i,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(orthog)(vprev,av,R1,R2,2,ierr),*ierr);
    v=av;
    }
  }


 // this function does an in-place Schur decomposition of T(1:m,1:m) into T and S.              
 // The Ritz values appear on the diagonal of T (for the real case there may be 2x2 blocks      
 // for complex conjugate pairs). The nselect and nsort flags indicate in which order they      
 // should appear:                                                                              
 // nselect=nsort=0: unsorted                                                                   
 // nselect>0: the first <nselect> Ritz values (if the last one is a complex conjugate pair     
 //            <nselect+1> appear in any order in the upper left corner of T                    
 // 0<nsort<nselect: in addition to moving a cluster of <nselect> eigenvalues, sort the         
 //             first <nsort> of them in the upper left corner.                                 
void SUBR(SchurDecomp)(_ST_* T, int ldT, _ST_* S, int ldS,
         int m, int nselect, int nsort, eigSort_t which, 
         std::complex<_MT_>* ev, int *ierr)
   {
   ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
   // this is for XGEES (computing the Schur form)
   int lwork = std::max(20*m,2*nselect*(m-nselect));
                          // min required workspace is 3*m 
                          // for GEES and 2*m*(m-nsort) for 
                          // TRSEN with condition estimate,
                          // so this should be enough for  
                          // good performance of GEES as.  
   ST work[lwork];
   // real and imag part of ritz values
   MT ev_r[m];   // in the complex case this is used as RWORK
   MT ev_i[m];

   const char *jobvs="V"; // compute the ritz vectors in S
   const char *sort="N";  // do not sort Ritz values (we do that later
                          // because gees only accepts the simple select
                          // function which does not compare the Ritz values)
  int sdim;

  // can select at most m Ritz values (and at least 0)
  nselect=std::max(0,std::min(nselect,m));
  // can sort at most nselect Ritz values (and at least 0)
  nsort=std::max(0,std::min(nsort,nselect));

#ifdef _IS_COMPLEX_
     PHIST_CHK_IERR(PREFIX(GEES)(jobvs,sort,NULL,&m,(blas_cmplx_t*)T,&ldT,
         &sdim,(blas_cmplx_t*)ev,(blas_cmplx_t*)S,&ldS,(blas_cmplx_t*)work,&lwork,ev_r,NULL,ierr),*ierr);
#else
     PHIST_CHK_IERR(PREFIX(GEES)(jobvs,sort,NULL,&m,T,&ldT,
         &sdim,ev_r,ev_i,S,&ldS,work,&lwork,NULL,ierr),*ierr);
     for (int i=0;i<m;i++)
       {
       ev[i]=std::complex<MT>(ev_r[i],ev_i[i]);
       }
#endif


   if (nselect==0) return;

   // find indices for the first howMany eigenvalues. A pair of complex conjugate
   // eigs is counted as a single one because we will skip solving the update equation
   // in that case. howMany is adjusted to include the pairs on output, for instance,
   // if howMany=1 on input but the first eig encountered is a complex conjugate pair,
   // the 2x2 block is shifted to the upper left of T and howMany=2 on output.
   int idx[m];

   // sort all eigenvalues according to 'which'.
   PHIST_CHK_IERR(SortEig(ev,m,idx,which,ierr),*ierr);

   // permute the first <nsort> eigenvalues according to idx
   // to the top left, taking the vectors along
   int select[m];
   for (int i=0;i<m;i++) select[i]=0;
   for (int i=0;i<nselect;i++) select[std::abs(idx[i])]=1;

  // call lapack routine to reorder Schur form
  const char *job="N"; // indicates wether we want condition estimates
                       // for [E]igenvalues, the invariant [S]ubspace or [B]oth
                       // (or [N]one, just sort)
  MT S_cond;
  MT sep;

#ifdef _IS_COMPLEX_
   PHIST_CHK_IERR(PREFIX(TRSEN)(job,jobvs,select,&m,(blas_cmplx_t*)T,&ldT,(blas_cmplx_t*)S,&ldS,(blas_cmplx_t*)ev,&nsort,
        &S_cond, &sep, (blas_cmplx_t*)work, &lwork, ierr),*ierr);
#else
   int liwork=nselect*(m-nselect);
   int iwork[liwork];
   PHIST_CHK_IERR(PREFIX(TRSEN)(job,jobvs,select,&m,T,&ldT,S,&ldS,ev_r,ev_i,&nsort,
        &S_cond, &sep, work, &lwork, iwork, &liwork, ierr),*ierr);   
   for (int i=0;i<m;i++)
     {
     ev[i]=std::complex<MT>(ev_r[i],ev_i[i]);
     }
#endif   

   for (int i=0;i<m;i++) select[i]=0;

   // sort all eigenvalues according to 'which'.
   PHIST_CHK_IERR(SortEig(ev,nselect,idx,which,ierr),*ierr);

   for (int i=0;i<nsort;i++)
     {
     // sort next candidate to top. Keep the select array in this loop
     // so that previoously sorted ones aren't moved back down.
     select[std::abs(idx[i])]=1;
#ifdef _IS_COMPLEX_
     PHIST_CHK_IERR(PREFIX(TRSEN)(job,jobvs,select,&m,(blas_cmplx_t*)T,&ldT,(blas_cmplx_t*)S,&ldS,
        (blas_cmplx_t*)ev,&nsort,&S_cond, &sep, (blas_cmplx_t*)work, &lwork, ierr),*ierr);
#else
     PHIST_CHK_IERR(PREFIX(TRSEN)(job,jobvs,select,&m,T,&ldT,S,&ldS,ev_r,ev_i,&nsort,
          &S_cond, &sep, work, &lwork, iwork, &liwork, ierr),*ierr);   
#endif
     }

#ifndef _IS_COMPLEX_
   if (nsort>0)
     {
     for (int i=0;i<m;i++)
       {
       ev[i]=std::complex<MT>(ev_r[i],ev_i[i]);
       }
     }
#endif
   }

