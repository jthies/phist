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
         int m, int nselect, int nsort, eigSort_t which, _MT_ tol, 
         void* v_ev, int *ierr)
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

  MT ev_r[m];   // in the complex case this is used as RWORK
#ifndef IS_COMPLEX
  // real and imag part of ritz values
  MT ev_i[m];
#endif
  const char *jobvs="V"; // compute the ritz vectors in S
  const char *sort="N";  // do not sort Ritz values (we do that later
                          // because gees only accepts the simple select
                          // function which does not compare the Ritz values)
  int sdim;
  CT* ev = (CT*)v_ev;

  // can select at most m Ritz values (and at least 0)
  nselect=std::max(0,std::min(nselect,m));
  // can sort at most nselect Ritz values (and at least 0)
  nsort=std::max(0,std::min(nsort,nselect));

  PHIST_DEB("m=%d, nselect=%d, nsort=%d",m,nselect,nsort);

#ifdef IS_COMPLEX
  PHIST_DEB("call complex %cGEES",st::type_char());
  PHIST_CHK_IERR(PREFIX(GEES)(jobvs,sort,NULL,&m,(blas_cmplx_t*)T,&ldT,
         &sdim,(blas_cmplx_t*)ev,(blas_cmplx_t*)S,&ldS,(blas_cmplx_t*)work,&lwork,ev_r,NULL,ierr),*ierr);
#else
  PHIST_DEB("call real %cGEES",st::type_char());
  PHIST_CHK_IERR(PREFIX(GEES)(jobvs,sort,NULL,&m,T,&ldT,
         &sdim,ev_r,ev_i,S,&ldS,work,&lwork,NULL,ierr),*ierr);
  for (int i=0;i<m;i++)
    {
    ev[i]=std::complex<MT>(ev_r[i],ev_i[i]);
    }
#endif

#if PHIST_OUTLEV>=PHIST_DEBUG
PHIST_OUT(0,"eigenvalues of unsorted Schur form:");
for (int i=0;i<m;i++)
  {
  PHIST_OUT(0,"%d\t%16.8g%+16.8gi",i,ct::real(ev[i]),ct::imag(ev[i]));
  }
#endif

  if (nselect<=0) return;
  if (nsort>nselect || nsort<0)
    {
    PHIST_OUT(PHIST_WARNING,"nselect=%d>=nsort=%d, or nsort>=0 "
                             "not satisfied, returning      unsorted Schur form",
                             nselect,nsort);
    *ierr=1;
    return;
    }

  // find indices for the first howMany eigenvalues. A pair of complex conjugate
  // eigs is counted as a single one because we will skip solving the update equation
  // in that case. howMany is adjusted to include the pairs on output, for instance,
  // if howMany=1 on input but the first eig encountered is a complex conjugate pair,
  // the 2x2 block is shifted to the upper left of T and howMany=2 on output.
  int idx[m];

  // permute the first <nselect> eigenvalues according to idx
  // to the top left, taking the vectors along
  int select[m];
  int nsorted=0;

  for (int i=0;i<m;i++) select[i]=0;

  // call lapack routine to reorder Schur form
  const char *job="N"; // indicates wether we want condition estimates
                       // for [E]igenvalues, the invariant [S]ubspace or [B]oth
                       // (or [N]one, just sort)
  const char *compq = "V"; // indicates wether we want to update the schur vectors
                           // [V]: update schur vectors in Q
                           // [N]: don't update schur vectors in Q
  MT S_cond;
  MT sep;

#ifndef IS_COMPLEX
  int liwork=m*m;
  int iwork[liwork];
#endif  
  if (nselect<m)
    {
    PHIST_DEB("initial sort step, nselect=%d",nselect);
    // sort all eigenvalues according to 'which'.
    PHIST_CHK_IERR(SortEig(ev,m,idx,which,tol,ierr),*ierr);
    for (int i=0;i<nselect;i++) 
      {
      select[std::abs(idx[i])]=1;
      }
#ifdef IS_COMPLEX
    PHIST_CHK_IERR(PREFIX(TRSEN)(job,compq,select,&m,(blas_cmplx_t*)T,&ldT,(blas_cmplx_t*)S,&ldS,(blas_cmplx_t*)ev,&nsorted,
          &S_cond, &sep, (blas_cmplx_t*)work, &lwork, ierr),*ierr);
#else
    PHIST_CHK_IERR(PREFIX(TRSEN)(job,compq,select,&m,T,&ldT,S,&ldS,ev_r,ev_i,&nsorted,
         &S_cond, &sep, work, &lwork, iwork, &liwork, ierr),*ierr);   
    for (int i=0;i<m;i++)
      {
      ev[i]=std::complex<MT>(ev_r[i],ev_i[i]);
      }
#endif   
    PHIST_DEB("nsorted=%d",nsorted);
    }//nselect<m
  
  if (nselect==1) return; // the one (or two for complex pairs) selected eigenvalue
                          // according to 'which' is already 'sorted'
  
  int i=0;
  while (i<nsort)
     {
     PHIST_DEB("sort step %d, nsorted=%d [%d]",i,nsorted,nsort);
     // sort the next few eigenvalues (up to nselect)
     PHIST_CHK_IERR(SortEig(ev+i,nselect-i,idx+i,which,tol,ierr),*ierr);

     // sort next candidate to top.
     for (int j=0;j<i;j++) select[j]=1;
     for (int j=i;j<m;j++) select[j]=0;
     select[std::abs(idx[i])+i]=1; // sort next one to the top
                                   // note that the index returned by 
                                   // SortEig misses an offset i because
                                   // we pass in ev+i
     int nsorted_before=nsorted;
#ifdef IS_COMPLEX
     PHIST_CHK_IERR(PREFIX(TRSEN)(job,compq,select,&m,(blas_cmplx_t*)T,&ldT,(blas_cmplx_t*)S,&ldS,
        (blas_cmplx_t*)ev,&nsorted,&S_cond, &sep, (blas_cmplx_t*)work, &lwork, ierr),*ierr);
#else
     PHIST_CHK_IERR(PREFIX(TRSEN)(job,compq,select,&m,T,&ldT,S,&ldS,ev_r,ev_i,&nsorted,
          &S_cond, &sep, work, &lwork, iwork, &liwork, ierr),*ierr);   
    for (int j=0;j<m;j++)
      {
      ev[j]=std::complex<MT>(ev_r[j],ev_i[j]);
      }
#endif
    i+= std::max(nsorted-nsorted_before,1);
    }//while

#if PHIST_OUTLEV>=PHIST_DEBUG
PHIST_OUT(0,"eigenvalues of sorted Schur form:");
for (int i=0;i<m;i++)
  {
  PHIST_OUT(0,"%d\t%16.8g%+16.8gi",i,ct::real(ev[i]),ct::imag(ev[i]));
  }
#endif


  }


// reorder multiple eigenvalues in a given (partial) schur decomposition by the smallest residual norm of the unprojected problem
void SUBR(ReorderPartialSchurDecomp)(_ST_* T, int ldT, _ST_* S, int ldS,
      int m, int nselected, _MT_ tol, _MT_* resNorm, void* ev, TYPE(sdMat_ptr) *transFormation, int *ierr)
{
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *ierr = 0;
}

