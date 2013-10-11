
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
void _SUBR_(arnoldi)(TYPE(const_op_ptr) op, TYPE(const_mvec_ptr) v0,
        TYPE(mvec_ptr) V, TYPE(sdMat_ptr) H, int m, int* ierr);

 // this function does a Schur decomposition of M(1:m,1:m) into T and S. The Ritz values
 // appear on the diagonal of T (for the real case there may be 2x2 blocks for complex
 // conjugate pairs). The nsort flag indicates in which order they should appear:
 // 0: unsorted
 // >0: first <nsort> Ritz values according to the 'which' flag (for instance the ones with
 // largest magnitude, smallest real part etc.) in the top left corner
 void SUBR(SchurDecomp)(_ST_ M, int ldM, _ST_* T, int ldT, _ST_* S, int ldS,
         int m, int nsort, int *ierr);

//! tries to compute a given number of eigenpairs (num_eig) of 
//! an non-hermitian operator A using the Jacobi-Davidson
//! method. 
//!
//! input arguments:
//!
//! A: pointer to the operator(-wrapper) A
//! X:  start vectors for the num_eigs desired
//!     eigenvectors (random vectors if you don't know better)
//! tol: convergence tolerance
//! num_eigs: number of desired eigenpairs
//! num_eigs: maximum number of iterations allowed
//! evals and resid should be pre-allocated arrays of size at least num_eigs
//! 
//! output arguments:
//!
//! num_eigs: number of converged eigenpairs
//! num_iters: number of iterations performed
//! X: latest eigenvector approximations
//! evals: latest eigenvalue approximations
//! resid: Ritz residuals
//! ierr: return code of the solver (0 on success, negative on error, positive on warning)
//!
//! TODO even for real matrices we may get pairs of complex eigenvalues/-vectors. These can be
//! stored as re/im in the given (real) mvec_ptr X, but we must indicate somehow to the caller
//! which candidates are complex conjugate pairs and which are real values/vectors.
void _SUBR_(jada)(_TYPE_(const_op_ptr) op, 
        _TYPE_(mvec_ptr) X,
        _ST_* evals, _MT_* resid, eigSort_t which,
        _MT_ tol,int* num_iters, int* num_eigs,
        int minBas, int maxBas,
        int* ierr)
  {
#include "phist_std_typedefs.hpp"

  // partial QR decomposition (converged eigenvectors form Q)
  mvec_ptr_t Q=NULL; // will be created as a view of the converged eigenvectors
  sdMat_ptr_t R=NULL;// will be allocated when needed.
  
  // current basis 
  mvec_ptr_t V, AV;
  // current update
  mvec_ptr_t t;
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
  
  int nconv=0; // number of converged eigenpairs
  int maxIter=*num_iters;
  int i,it,n,m;
  ST alpha, beta;
  int converged[maxIter];
  MT r_est[maxIter];
  MT nrm;

  *num_iters=0;
  
  const_comm_ptr_t comm;
  PHIST_CHK_IERR(phist_map_get_comm(op->range_map_,&comm,ierr),*ierr);

  PHIST_CHK_IERR(_SUBR_(mvec_create)(&V,op->domain_map_,maxBas,ierr),*ierr);
  PHIST_CHK_IERR(_SUBR_(mvec_create)(&AV,op->domain_map_,maxBas,ierr),*ierr);

  PHIST_CHK_IERR(SUBR(sdMat_create)(&M,maxBas,maxBas,comm,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_create)(&S,maxBas,maxBas,comm,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_create)(&T,maxBas,maxBas,comm,ierr),*ierr);

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
      PHIST_CHK_IERR(op->apply(st::one(),op->A_,Vm,st::zero(),AVm);
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

    // unless a restart is impending, sort only the next Ritz value to the top left
    int sort = m==maxBas-1? minBas: 1;
    PHIST_CHK_IERR(SUBR(SchurDecomp)(M_ptr,ldM,T_ptr,ldT,S_ptr,ldS,m,sort,ierr),*ierr);

    theta=T_ptr[0]; // first diagonal entry of T is next EV to go for.
                    // In the real case it may be a complex eigenvalue,
                    // but we just take the real part then.
    
    // get Ritz vector corresponding to theta
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(S,&Sblock,0,m,0,0,ierr),*ierr);
    
    //u=V*s;
    //...
    
    //Au=AV*s;
    //...

    //r=Au-theta*u;
    //...

    //atil = Q'*r;
    //...
    
    //rtil = r-Q*atil;
    //...
    
    //nrm=norm(rtil);
    //...
    
  PHIST_OUT(1,"%d\t%d\t%8.4g %8.4gi\t\t%8.4g\n",it,m,st::real(theta),st::imag(theta),nrm);

  // deflate converged eigenpairs
  while (nrm<=tol)
    {
    R=[R, atil;
       zeros(1,k), theta];
    Q=[Q,u];
    k=k+1;
    mm=0;
    if (verbose)
      disp(['eigenvalue ',int2str(k),' (',num2str(theta),') is converged.']);
      %more on;
    end %if
    if (debug)
      disp('updated Schur form (R):');
      R
    end
    if (k==numEigs)
      SOLVE=false;
      break;
    end

%    if (abs(imag(theta))>tol)
%      t=imag(u/sign(max(u)));
%      if (norm(t)>tol)
%        t=orthog(Qschur,t); 
%        EXPAND=(norm(t)>sqrt(tol)); 
%      end
%    end
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

    [S,T]=SortSchur(M,target,m==mmax,mmin);

  end % while (deflate)

  % restart if necessary
  if m>=mmax
    if (verbose)
      disp('restart JDQR');
    end
    m=mmin;
    S=S(:,1:m);
    M=T(1:m,1:m);
    V=V*S;
    AV=AV*S;
    V=V(:,1:mmin);
    AV=AV(:,1:mmin);
    SOLVE=true;
    EXPAND=true;
  end
  if (solve)
    {
    // maintain orthogonality against
    // all converged vectors (Q) and the
    // new one:
    Qtil=[Q,u];
    if (nrm<switchTol)
      if (verbose)
        disp('using RQI');
      end
      shift=theta;
    else
      if (verbose)
        disp('using SI')
      end
      shift=tau;
    end
    % solve approximately 
    % (I-uu')(A-theta*I)(I-uu')*t=-r
    % to get t \orth u (u ^= Qtil here)
    if (strcmp(lsFun,'direct'))
      A_aug = [A-shift*speye(n), Qtil;
                Qtil'        , zeros(size(Qtil,2))];
      b_aug = [rtil;zeros(size(Qtil,2),size(rtil,2))];
      t_aug= A_aug\b_aug;
      t=t_aug(1:size(Qtil,1),:);
    else
      op=comp_jada_op(A,shift,speye(n),Qtil);
      if (isfield(printOpts,'indent'))
        printOpts.indent=printOpts.indent+1;
      else
        printOpts.indent=1;
      end
      % TODO - is this necessary?
      precOp=precOp.compute(A-shift*speye(n),precOp);
      t0=zeros(size(rtil));
      lsOpts.tol=max(tol,1/2.^max(1,mm-1));
      if (verbose)
        disp(['inner conv tol: ',num2str(lsOpts.tol)]);
      end
      [t,flag,relres,iter,resvec] = ...
        lsFun(op,rtil,t0,lsOpts,precOp);
      
      printOpts.indent=printOpts.indent-1;
    end % direct or user-supplied
    EXPAND=true;    
    if (verbose)
      disp(['(Qtil,t)=',num2str(norm(Qtil'*t))]);
      disp(['(A-sI)t+rtil=',num2str(norm(A*t-shift*t+rtil))]);         
        }
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
  
  PHIST_CHK_IERR(_SUBR_(mvec_delete)(vold,ierr),*ierr);
  PHIST_CHK_IERR(_SUBR_(mvec_delete)(vnew,ierr),*ierr);
  PHIST_CHK_IERR(_SUBR_(sdMat_delete)(S,ierr),*ierr);
  
  free(alphas);
  free(betas);
  free(falphas);
  free(fbetas);
  free(work);

  return;
  }



void _SUBR_(arnoldi)(TYPE(const_op_ptr) op, TYPE(const_mvec_ptr) v0,
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


 // this function does a Schur decomposition of M(1:m,1:m) into T and S. The Ritz values
 // appear on the diagonal of T (for the real case there may be 2x2 blocks for complex
 // conjugate pairs). The sort flag indicates in which order they should appear:
 // 0: unsorted
 // >0: first <sort> Ritz values according to the 'which' flag (for instance the ones with
 // largest magnitude, smallest real part etc.) in the top left corner
 void SUBR(SchurDecomp)(_ST_ M, int ldM, _ST_* T, int ldT, _ST_* S, int ldS,
         int m, int nsort, eigSort_t which, int *ierr)
   {
#include "phist_std_typedefs.hpp"
   // this is for XGEES (computing the Schur form)
   int lwork = std::max(20*m,2*nsort*(m-nsort));
                          // min required workspace is 3*m 
                          // for GEES and 2*m*(m-nsort) for 
                          // TRSEN with condition estimate,
                          // so this should be enough for  
                          // good performance of GEES as.  
   ST work[lwork];
   // real and imag part of ritz values
   MT ev_r[m];   // in the complex case this is used as RWORK
   MT ev_i[m];
   std::complex<MT> ev[m];

   const char *jobvs="V"; // compute the ritz vectors in S
   const char *sort="N";  // do not sort Ritz values (we do that later
                          // because gees only accepts the simple select
                          // function which does not compare the Ritz values)

 #ifdef IS_COMPLEX
     PHIST_CHK_IERR(PREFIX(GEES)(jobvs,sort,NULL,&m,T_ptr,&ldT,
         &sdim,ev,S_ptr,&ldS,work,&lwork,ev_r,NULL,ierr),*ierr);
 #else
     PHIST_CHK_IERR(PREFIX(GEES)(jobvs,sort,NULL,&m,T_ptr,&ldT,
         &sdim,ev_r,ev_i,S_ptr,&ldS,work,&lwork,NULL,ierr),*ierr);
#endif

     if (!nsort) return;

#ifdef IS_COMPLEX
     for (int i=0;i<m;i++)
       {
       ev[i]=std::complex<MT>(ev_r[i],ev_i[i]);
       }
#endif

   // find indices for the first howMany eigenvalues. A pair of complex conjugate
   // eigs is counted as a single one because we will skip solving the update equation
   // in that case. howMany is adjusted to include the pairs on output, for instance,
   // if howMany=1 on input but the first eig encountered is a complex conjugate pair,
   // the 2x2 block is shifted to the upper left of T and howMany=2 on output.
   int idx[m];

   // sort all eigenvalues according to 'which'.
   PHIST_CHK_IERR(SortEig(ev,idx,which,ierr),*ierr);

   // permute the first <nsort> eigenvalues according to idx
   // to the top left, taking the vectors along
   int select[m];
   for (int i=0;i<m;i++) select[i]=0;
   for (int i=0;i<nsort;i++) select[idx[i]]=1;

  // call lapack routine to reorder Schur form
  const char *job="N"; // indicates wether we want condition estimates
                       // for [E]igenvalues, the invariant [S]ubspace or [B]oth
                       // (or [N]one, just sort)
  MT S[nsort];
  MT* sep[nsort];

#ifdef IS_COMPLEX
   PHIST_CHK_IERR(PREFIX(TRSEN)(job,jobvs,select,&m,T_ptr,&ldT,S_ptr,&ldS,ev,&nsort,
        S, sep, work, &lwork, ierr),*ierr);
#else
   int liwork=nsort*(m-nsort);
   int iwork[liwork];
   PHIST_CHK_IERR(PREFIX(TRSEN)(job,jobvs,select,&m,T_ptr,&ldT,S_ptr,&ldS,ev_r,ev_i,&nsort,
        S, sep, work, &lwork, iwork, &liwork, ierr),*ierr);   
#endif   
   }

