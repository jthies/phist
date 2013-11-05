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

   // sort all eigenvalues according to 'which'.
   PHIST_CHK_IERR(SortEig(ev,m,idx,which,ierr),*ierr);

   // permute the first <nselect> eigenvalues according to idx
   // to the top left, taking the vectors along
   int select[m];
   for (int i=0;i<m;i++) select[i]=0;

  // call lapack routine to reorder Schur form
  const char *job="N"; // indicates wether we want condition estimates
                       // for [E]igenvalues, the invariant [S]ubspace or [B]oth
                       // (or [N]one, just sort)
  MT S_cond;
  MT sep;

  int liwork=std::max(1,nselect*(m-nselect));
  int iwork[liwork];
  
  if (nselect<m)
    {
    PHIST_DEB("first sort step, nselect=%d",nselect);
    for (int i=0;i<nselect;i++) 
      {
      select[std::abs(idx[i])]=1;
      }
#ifdef _IS_COMPLEX_
    PHIST_CHK_IERR(PREFIX(TRSEN)(job,jobvs,select,&m,(blas_cmplx_t*)T,&ldT,(blas_cmplx_t*)S,&ldS,(blas_cmplx_t*)ev,&nsort,
          &S_cond, &sep, (blas_cmplx_t*)work, &lwork, ierr),*ierr);
#else
    PHIST_CHK_IERR(PREFIX(TRSEN)(job,jobvs,select,&m,T,&ldT,S,&ldS,ev_r,ev_i,&nsort,
         &S_cond, &sep, work, &lwork, iwork, &liwork, ierr),*ierr);   
    for (int i=0;i<m;i++)
      {
      ev[i]=std::complex<MT>(ev_r[i],ev_i[i]);
      }
#endif   

    // sort the first nselect eigenvalues according to 'which'.
    PHIST_CHK_IERR(SortEig(ev,nselect,idx,which,ierr),*ierr);
    }//nselect<m

  for (int i=0;i<m;i++) select[i]=0;

  for (int i=0;i<nsort;i++)
     {
     PHIST_DEB("sort step, %d [%d]",i,nsort);
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

