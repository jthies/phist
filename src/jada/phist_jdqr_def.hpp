
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
 //            <nselect+1> appear in any order in the upper left corner of T                    
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
         int m, int nselect, int nsort, std::complex<MT>* ev, int *ierr);

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
        TYPE(mvec_ptr) X, eigSort_t which,
        _ST_* evals, _MT_* resid, int* is_cmplx,
        _MT_ tol, int* num_eigs, int* num_iters,
        int minBas, int maxBas,
        int* ierr)
  {
#include "phist_std_typedefs.hpp"

  int nconv=0; // number of converged eigenpairs

  // partial QR decomposition (converged eigenvectors form Q)
  mvec_ptr_t Q=NULL;  // will be created as a view of the converged eigenvectors
  sdMat_ptr_t R;      // store R part of QR decomp
  sdMat_ptr_t Rblock; // to view the first nconv x nconv block of R
  
  // current basis 
  mvec_ptr_t V, AV;
  // current update
  mvec_ptr_t t;
  // some more (single) vectors we need:
  mvec_ptr_t u, Au, r, rtil;
  // the above vectors may be complex even for a real
  // matrix. These are real- and imaginary parts, respectively
  // (in the complex case x_r aliases x and x_i==NULL)
  mvec_ptr_t u_r, u_i, Au_r, Au_i, r_r, r_i, rtil_r, rtil_i;
  // these point either to u or u_r etc. to make live simpler further down
  mvec_ptr_t u_ptr, Au_ptr, r_ptr, rtil_ptr;
  // Q*s (temporary vector)
  sdMat_ptr_t atil, atil_r, atil_i;
  // matrix (V,AV)
  sdMat_ptr_t M;
  sdMat_ptr_t Mblock; // to create views of parts of M
  // Schur-decomposition of M
  sdMat_ptr_t T, S;
  sdMat_ptr_t Tblock,Sblock; // to create views of parts of T and S

  // for extracting views of M, T and S (to call lapack etc.)
  int ldM,ldS,ldT;
  ST *M_ptr, *S_ptr, *T_ptr; 
  
  //! view of certain columns of V
  mvec_ptr_t Vblock, Vm;
  //! view of certain columns of AV
  mvec_ptr_t AVblock, AVm;

  //! Hessenberg matrix from the initial Arnoldi iteration
  sdMat_ptr_t H0;
  
  bool expand=true;
  bool solve=true;
  std::complex<MT> theta; // next eigenvalue to go for
  // theta as a sdMat_t (either 1x1 or 2x2 in real case)
  sdMat_ptr_t Theta;
  std::complex<MT> ev[*maxBas];
  
  int maxIter=*num_iters;
  int i,it,m,mm;
  int converged[maxIter];
  MT r_est[maxIter];
  MT nrm[2];

  *num_iters=0;
  
  const_comm_ptr_t comm;
  PHIST_CHK_IERR(phist_map_get_comm(op->range_map_,&comm,ierr),*ierr);

  PHIST_CHK_IERR(SUBR(mvec_create)(&V,op->domain_map_,maxBas,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_create)(&AV,op->domain_map_,maxBas,ierr),*ierr);

#ifdef IS_COMPLEX
  PHIST_CHK_IERR(SUBR(mvec_create)(&u,op->domain_map_,1,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_create)(&Au,op->domain_map_,1,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_create)(&r,op->domain_map_,1,ierr),*ierr);
  
  u_r=u; u_i=NULL;
  Au_r=Au; Au_i=NULL;
  r_r=r; r_i=NULL;
  rtil_r=rtil; rtil_i=NULL;
#else
  PHIST_CHK_IERR(SUBR(mvec_create)(&u,op->domain_map_,2,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_create)(&Au,op->domain_map_,2,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_create)(&r,op->domain_map_,2,ierr),*ierr);
  
  PHIST_CHK_IERR(SUBR(mvec_view_block)(u,&u_r,0,0,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_view_block)(u,&u_i,1,1,ierr),*ierr);

  PHIST_CHK_IERR(SUBR(mvec_view_block)(Au,&Au_r,0,0,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_view_block)(Au,&Au_i,1,1,ierr),*ierr);

  PHIST_CHK_IERR(SUBR(mvec_view_block)(r,&r_r,0,0,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_view_block)(r,&r_i,1,1,ierr),*ierr);

  PHIST_CHK_IERR(SUBR(mvec_view_block)(rtil,&rtil_r,0,0,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_view_block)(rtil,&rtil_i,1,1,ierr),*ierr);  
#endif
  PHIST_CHK_IERR(SUBR(sdMat_create)(&M,maxBas,maxBas,comm,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_create)(&S,maxBas,maxBas,comm,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_create)(&T,maxBas,maxBas,comm,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_create)(&R,*num_eigs,*num_eigs,comm,ierr),*ierr);

  // pointer to the data in M, S and T
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(M,&M_ptr, &ldM,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(S,&S_ptr, &ldS,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(T,&T_ptr, &ldT,ierr),*ierr);

  PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&Vblock,0,minBas,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_view_block)(M,&H0, 0,minBas,0,minBas,ierr),*ierr);

  PHIST_OUT(1,"%d steps of Arnoldi as start-up",minBas);
  
  // start by filling the first minBas vectors using Arnoldi
  PHIST_CHK_IERR(SUBR(arnoldi)(op,v0,Vblock,H0,minBas,ierr),*ierr);

  PHIST_CHK_IERR(SUBR(mvec_view_block)(AV,&AVblock,0,minBas-1,ierr),*ierr);
  // compute A*V for the first minBas columns which we have now using the Arnoldi relation
  PHIST_CHK_IERR(mvec_times_sdMat)(st::one(),Vblock,H0,st::zero(),AVblock,ierr),*ierr);
  m=minBas-1;
  expand=false;
  
  PHIST_OUT(1,"Jacobi-Davidson");
  PHIST_OUT(1,"%s\t%s\t%s\t\t%s\n","iter","m","approx","resid");

  it=0;
  int mm=0;
  while (nconv<numEigs && it < maxIter)
    {
    //TODO - thoroughly check the switch from 1-based m in MATLAB to 0-based here
    if (expand)
      {
      m++;
      mm++;
      it++;
      // create views of V(:,1:m-1), V(:,m) and AV likewise
      PHIST_CHK_IERR(SUBR(mvec_extract_view)(V,&Vblock,0,m-1,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(mvec_extract_view)(AV,&AVblock,0,m-1,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(mvec_extract_view)(V,&Vm,m,m,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(mvec_extract_view)(AV,&AVm,m,m,ierr),*ierr);
      // orthogonalize t against V(:,0:m-1)
      PHIST_CHK_IERR(SUBR(orthog)(Vblock,t,ierr),*ierr);
      // set V(:,m)=t
      PHIST_CHK_IERR(SUBR(mvec_set_block)(V,t,m,m,ierr),*ierr);
      // compute AV(:,m) = A*V(:,m)
      PHIST_CHK_IERR(op->apply(st::one(),op->A_,Vm,st::zero(),AVm,ierr),*ierr);
      // Galerkin for non-Hermitian A
      // TODO - is it maybe more efficient to just compute V'AV completely? It
      //        means more ops and more data transfer but fewer messages.
      // M(1:m-1,m)=V(:,1:m-1)'*AV(:,m);
      PHIST_CHK_IERR(SUBR(sdMat_view_block)(M,&Mblock,0,m-1,m,m,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),Vblock,AVm,st::zero,Mblock,ierr),*ierr);
      // M(m,1:m-1)=V(:,m)'*AV(:,1:m-1);
      PHIST_CHK_IERR(SUBR(sdMat_view_block)(M,&Mblock,m,m,0,m-1,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),Vm,AVblock,st::zero,Mblock,ierr),*ierr);
      // M(m,m)=V(:,m)'*AV(:,m) note that we could do this using a dot product in the
      // single-vector case, but we want to extend it to a block variant so we use the more
      // general dense matmul.
      PHIST_CHK_IERR(SUBR(sdMat_view_block)(M,&Mblock,m,m,m,m,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),Vm,AVm,st::zero,Mblock,ierr),*ierr);
      }
    else
      {
      expand=true;
      }
      
    // now do a Schur decomposition of M into S and T, with the
    // Ritz values on the diagonal of T.

    // copy T=M. Note that as we just 'view' the upper left blocks,
    // the pointers M_ptr, T_ptr, S_ptr are still correct
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(M,&Mblock,0,m,0,m,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(T,&Tblock,0,m,0,m,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),Mblock,st::zero(),Tblock,ierr),*ierr);

    // if a restart is impending, shift the next <minBas> Ritz values
    // to the top left. Shift the next <nsort>=1 one of them to the 
    // top left (next shift to choose).
    int nselect = m==maxBas-1? minBas: 1;
    int nsort = 1;
    PHIST_CHK_IERR(SUBR(SchurDecomp)
        (T_ptr,ldT,S_ptr,ldS,m,nselect,nsort,ev,ierr),*ierr);
    theta=ev[0];

    int i0=0;
    int i1=0; // i1==1 indicates complex eigenvalue of real matrix, we
              // will use this flag later on.
    // set some pointers
    u_ptr=u;
    Au_ptr=Au;
    r_ptr=r;
    rtil_ptr=rtil;
#ifndef IS_COMPLEX
    if (mt::abs(ct::imag(theta))>mt::eps())
      {
      i1++; // get imaginary part
      }
    else
      {
      // real matrix, real Ritz theta
      u_ptr=u_r;
      Au_ptr=Au_r;
      r_ptr=r_r;
      rtil_ptr=rtil_r;
      // TODO Theta as 1x1 matrix?
      }
#endif

    // get the diagonal block (1x1 or 2x2) corresponding to theta
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(T,&Theta,i0,i1,i0,i1,ierr),*ierr);

    // get Ritz vector corresponding to theta
    // (for complex RV in real case, two columns are extracted)
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(S,&Sblock,0,m,i0,i1,ierr),*ierr);

    //u=V*s;
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(st::one(),V,Sblock,st::zero(),u_ptr,ierr),*ierr);
    //Au=AV*s;
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(st::one(),AV,Sblock,st::zero(),Au_ptr,ierr),*ierr);

    // r=Au-theta*u; ([r_r, r_i] = [Au_r, Au_i] - [u_r, u_i]*Theta in the real case with 
    // complex theta)
    
    // first: r=Au
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),Au_ptr,st::zero(),r_ptr,ierr),*ierr);

    // update r = r - U*Theta
    PHIST_CHK_IERR(SUBR(mvec_times_sdMAt)(-st::one(),u_ptr,Theta,st::one(),r_ptr,ierr),*ierr);

    // set rtil=r
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),r_ptr,st::zero(),rtil_ptr,ierr),*ierr);

    if (nrowsQ>0)
      {
      //atil = Q'*r;
      PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),Q,r,st::zero(),atil,ierr),*ierr);
    
      //rtil = r-Q*atil;
      PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(-st::one(),Q,atil,st::one(),r,ierr),*ierr);
      }

    //nrm=norm(rtil);
    // real case with complex r: ||v+iw||=sqrt((v+iw).'*(v-iw))=sqrt(v'v+w'w)
    nrm[1]=mt::zero();
    PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(rtil,rtil,nrmv,ierr),*ierr);
    nrm[0]=mt::sqrt(nrm[0]*nrm[0]+nrm[1]*nrm[1]);
    PHIST_OUT(1,"%d\t%d\t%8.4g %8.4gi\t\t%8.4g\n",it,m,st::real(theta),st::imag(theta),nrm);

  // deflate converged eigenpairs
  while (nrm<=tol)
    {
//    R=[R, atil;
//       0, theta];
// ...
//    Q=[Q,u];
// ...
    k=k+1;
    mm=0;
    PHIST_OUT(1,"eigenvalue %d (%g+%gi) is converged.",k,st::real(theta),st::imag(theta));
    if (k==numEigs)
      {
      solve=false;
      break;
      }

    if (abs(st::imag(theta))>tol)
       {
       t=imag(u/sign(max(u)));
       if (norm(t)>tol)
         {
         t=orthog(Qschur,t); 
         expand=(norm(t)>sqrt(tol)); 
         }
       }
    // TROET! This requires theta to be at pos 0
    S=S(:,2:m);
    M=T(2:m,2:m);
    m=m-1;
    V=V*S;
    AV=AV*S;
    u=V(:,1);
    Au=AV(:,1);
    theta=M(1,1);
    r=Au-theta*u;
    atil=Q'*r;
    rtil=r-Q*atil;
    nrm=norm(rtil);
    // here we sort again
    [S,T]=SortSchur(M,target,m==mmax,mmin);
    } //while (deflate)

  // restart if necessary
  if (m>=mmax)
    {
    PHIST_OUT(1,"restart JDQR");
    m=mmin;
    S=S(:,1:m);
    M=T(1:m,1:m);
    V=V*S;
    AV=AV*S;
    V=V(:,1:mmin);
    AV=AV(:,1:mmin);
    solve=true;
    expand=true;
    }
  if (solve)
    {
    // maintain orthogonality against
    // all converged vectors (Q) and the
    // new one:
    Qtil=[Q,u];
    if (nrm<switchTol)
      {
      PHIST_OUT(1,"using RQI");
      shift=theta;
      }
    else
      {
      PHIST_OUT(1,"using SI");
      shift=target;
      }
    // solve approximately 
    // (I-uu')(A-theta*I)(I-uu')*t=-r
    // to get t \orth u (u ^= Qtil here)

      op=comp_jada_op(A,shift,speye(n),Qtil);
      //precOp=precOp.compute(A-shift*speye(n),precOp);
      t0=zeros(size(rtil));
      lsOpts.tol=max(tol,1/2.^max(1,mm-1));

      PHIST_OUT(1,"inner conv tol: %f",lsTol);

      [t,flag,relres,iter,resvec] = ...
        lsFun(op,rtil,t0,lsOpts,precOp);
      
    expand=true;
//    PHIST_OUT(1,"(Qtil,t)=%f",num2str(norm(Qtil'*t)));
//    PHIST_OUT(1,"(A-sI)t+rtil=%f",norm(A*t-shift*t+rtil));         
    }
  else
    {
    expand=false;
    solve=true;
    }// if solve
  }// while loop
    
  *num_iters=it;
  inext=0;
  // copy at most num_eigs converged eigenvalues into the user
  // provided array
  for (i=n-1;i>=0, inext<*num_eigs;i--)
    {
    if (converged[i])
      {
      evals[inext] = falphas[i];
      resid[inext] = r_est[i];
      inext++;
      }
    }
  if (nconv<*num_eigs) *num_eigs=nconv;
  
  // TODO - compute eigenvectors (need to store the whole basis V for that)
  
  PHIST_CHK_IERR(SUBR(mvec_delete)(vold,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_delete)(vnew,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(S,ierr),*ierr);
  
  free(alphas);
  free(betas);
  free(falphas);
  free(fbetas);
  free(work);

  return;
  }



void SUBR(arnoldi)(TYPE(const_op_ptr) op, TYPE(const_mvec_ptr) v0,
        TYPE(mvec_ptr) V, TYPE(sdMat_ptr) H, int m, int* ierr)
  {
#include "phist_std_typedefs.hpp"
  TYPE(mvec_ptr) v,av,vprev;
  TYPE(sdMat_ptr) R1,R2;
  PHIST_CHK_IERR(SUBR(mvec_set_block)(V,v0,0,0,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_view)(V,v,0,0,ierr),*ierr);
  
  for (int i=0;i<m;i++)
    {
    PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&vprev,0,i,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&av,i+1,i+1,ierr),*ierr);
    PHIST_CHK_IERR(op->apply(st::one(),op->A_,v,st::zero(),av,ierr),*ierr);
    // orthogonalize
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(H,&R1,0,i,i,i,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(H,&R2,i+1,i+1,i,i,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(orthog)(vprev,av,R1,R2,ierr),*ierr);
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
         std::complex<MT>* ev, int *ierr)
   {
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


  // can select at most m Ritz values (and at least 0)
  nselect=std::max(0,std::min(nselect,m));
  // can sort at most nselect Ritz values (and at least 0)
  nsort=std::max(0,std::min(nsort,nselect));

#ifdef IS_COMPLEX
     PHIST_CHK_IERR(PREFIX(GEES)(jobvs,sort,NULL,&m,T_ptr,&ldT,
         &sdim,ev,S_ptr,&ldS,work,&lwork,ev_r,NULL,ierr),*ierr);
#else
     PHIST_CHK_IERR(PREFIX(GEES)(jobvs,sort,NULL,&m,T_ptr,&ldT,
         &sdim,ev_r,ev_i,S_ptr,&ldS,work,&lwork,NULL,ierr),*ierr);
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
  MT S[nselect];
  MT* sep[nselect];

#ifdef IS_COMPLEX
   PHIST_CHK_IERR(PREFIX(TRSEN)(job,jobvs,select,&m,T_ptr,&ldT,S_ptr,&ldS,ev,&nsort,
        S, sep, work, &lwork, ierr),*ierr);
#else
   int liwork=nselect*(m-nselect);
   int iwork[liwork];
   PHIST_CHK_IERR(PREFIX(TRSEN)(job,jobvs,select,&m,T_ptr,&ldT,S_ptr,&ldS,ev_r,ev_i,&nsort,
        S, sep, work, &lwork, iwork, &liwork, ierr),*ierr);   
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
#ifdef IS_COMPLEX
     PHIST_CHK_IERR(PREFIX(TRSEN)(job,jobvs,select,&m,T_ptr,&ldT,S_ptr,&ldS,ev,&nsort,
          S, sep, work, &lwork, ierr),*ierr);
#else
     PHIST_CHK_IERR(PREFIX(TRSEN)(job,jobvs,select,&m,T_ptr,&ldT,S_ptr,&ldS,ev_r,ev_i,&nsort,
          S, sep, work, &lwork, iwork, &liwork, ierr),*ierr);   
#endif   
     }

#ifndef IS_COMPLEX
   if (nsort>0)
     {
     for (int i=0;i<m;i++)
       {
       ev[i]=std::complex<MT>(ev_r[i],ev_i[i]);
       }
     }
#endif
   }

