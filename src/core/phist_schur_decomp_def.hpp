/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
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
         int m, int nselect, int nsort, phist_EeigSort which, _MT_ tol, 
         void* v_ev, int *iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag = 0;
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN_SMALLDETERMINISTIC(ComputeTask)

  MT ev_r[m];   // in the complex case this is used as RWORK
#ifndef IS_COMPLEX
  // real and imag part of ritz values
  MT ev_i[m];
#endif
  phist_blas_char jobvs='V'; // compute the ritz vectors in S
  phist_blas_char sort='N';  // do not sort Ritz values (we do that later
                          // because gees only accepts the simple select
                          // function which does not compare the Ritz values)
  int sdim;
  CT* ev = (CT*)v_ev;

  // can select at most m Ritz values (and at least 0)
  nselect=std::max(0,std::min(nselect,m));
  // can sort at most nselect Ritz values (and at least 0)
  nsort=std::max(0,std::min(nsort,nselect));

  PHIST_DEB("m=%d, nselect=%d, nsort=%d\n",m,nselect,nsort);

// prohibit parallel execution to assure identical results on different procs
#pragma omp parallel
  {
#pragma omp master
    {
#ifdef IS_COMPLEX
      PHIST_DEB("call complex %cGEES\n",st::type_char());
      *iflag=PHIST_TG_PREFIX(GEES)(SDMAT_FLAG,jobvs,sort,NULL,m,(blas_cmplx*)T,ldT,
             &sdim,(blas_cmplx*)ev,(blas_cmplx*)S,ldS);
#else
      PHIST_DEB("call real %cGEES\n",st::type_char());
      *iflag=PHIST_TG_PREFIX(GEES)(SDMAT_FLAG,jobvs,sort,NULL,m,T,ldT,
            &sdim,ev_r,ev_i,S,ldS);
      for (int i=0;i<m;i++)
      {
        ev[i]=std::complex<MT>(ev_r[i],ev_i[i]);
      }
#endif
    }
  }
  if (*iflag>0)
  {
    PHIST_SOUT(PHIST_WARNING,"Lapack routine xGEES returned positive value INfo=%d, input size was m=%d\n",*iflag,m);
  }
  PHIST_CHK_NEG_IERR(PHIST_TOUCH("GEES"),*iflag);

#if PHIST_OUTLEV>=PHIST_DEBUG
//PHIST_OUT(0,"eigenvalues of unsorted Schur form:\n");
//for (int i=0;i<m;i++)
//{
  //PHIST_OUT(0,"%d\t%16.8g%+16.8gi\n",i,ct::real(ev[i]),ct::imag(ev[i]));
//}
#endif

  if (nselect<=0) return;
  if (nsort>nselect || nsort<0)
  {
    PHIST_OUT(PHIST_WARNING,"nselect=%d>=nsort=%d, or nsort>=0 "
                             "not satisfied, returning      unsorted Schur form\n",
                             nselect,nsort);
    *iflag=1;
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
  phist_blas_char job='N'; // indicates wether we want condition estimates
                       // for [E]igenvalues, the invariant [S]ubspace or [B]oth
                       // (or [N]one, just sort)
  phist_blas_char compq = 'V'; // indicates wether we want to update the schur vectors
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
    PHIST_DEB("initial sort step, nselect=%d\n",nselect);
    // sort all eigenvalues according to 'which'.
    PHIST_CHK_IERR(SortEig(ev,m,idx,which,tol,iflag),*iflag);
    for (int i=0;i<nselect;i++) 
      select[std::abs(idx[i])]=1;
// prohibit parallel execution to assure identical results on different procs
#pragma omp parallel
    {
#pragma omp master
      {
#ifdef IS_COMPLEX
        *iflag=PHIST_TG_PREFIX(TRSEN)(SDMAT_FLAG,job,compq,select,m,(blas_cmplx*)T,ldT,(blas_cmplx*)S,ldS,(blas_cmplx*)ev,&nsorted,
              &S_cond, &sep);
#else
        *iflag=PHIST_TG_PREFIX(TRSEN)(SDMAT_FLAG,job,compq,select,m,T,ldT,S,ldS,ev_r,ev_i,&nsorted,
              &S_cond, &sep);
        for (int i=0;i<m;i++)
        {
          ev[i]=std::complex<MT>(ev_r[i],ev_i[i]);
        }
#endif   
      }
    }
    PHIST_CHK_IERR(;,*iflag);
    PHIST_DEB("nsorted=%d\n",nsorted);
    // *POSSIBLE PROBLEM*
    // if we select part of a complex conjugate eigenpair, trsen just increases nselect by one
    if( nsorted > nselect )
    {
      PHIST_DEB("detected nsorted > nselect, try to continue with nsort <- nselect and nselect <- nsorted\n");
      nsort = nselect;
      nselect = nsorted;
    }
  }//nselect<m
  
  if (nselect==1) return; // the one (or two for complex pairs) selected eigenvalue
                          // according to 'which' is already 'sorted'

// prohibit parallel execution to assure identical results on different procs
#pragma omp parallel
  {
#pragma omp master
    {
      int i=0;
      while (i<nsort)
      {
        PHIST_DEB("sort step %d, nsorted=%d [%d]\n",i,nsorted,nsort);
        // sort the next few eigenvalues (up to nselect)
        SortEig(ev+i,nselect-i,idx+i,which,tol,iflag);
        if( *iflag != 0 )
          break;

        // sort next candidate to top.
        for (int j=0;j<i;j++) select[j]=1;
        for (int j=i;j<m;j++) select[j]=0;
        select[std::abs(idx[i])+i]=1; // sort next one to the top
        // note that the index returned by 
        // SortEig misses an offset i because
        // we pass in ev+i
        int nsorted_before=nsorted;
#ifdef IS_COMPLEX
        PHIST_TG_PREFIX(TRSEN)(SDMAT_FLAG,job,compq,select,m,(blas_cmplx*)T,ldT,(blas_cmplx*)S,ldS,
              (blas_cmplx*)ev,&nsorted,&S_cond, &sep);
#else
        *iflag=PHIST_TG_PREFIX(TRSEN)(SDMAT_FLAG,job,compq,select,m,T,ldT,S,ldS,ev_r,ev_i,&nsorted,
              &S_cond, &sep);
        for (int j=0;j<m;j++)
        {
          ev[j]=std::complex<MT>(ev_r[j],ev_i[j]);
        }
#endif
        if( *iflag != 0 )
          break;
        i+= std::max(nsorted-nsorted_before,1);
      }//while
#if PHIST_OUTLEV>=PHIST_DEBUG
      //PHIST_OUT(0,"eigenvalues of sorted Schur form:\n");
      //for (int i=0;i<m;i++)
      //{
        //PHIST_OUT(0,"%d\t%16.8g%+16.8gi\n",i,ct::real(ev[i]),ct::imag(ev[i]));
      //}
#endif
    }
  }
PHIST_TASK_END(iflag)
}

// generalized Schur Decomposition, (S,T)->(~S,~T,VS,WS) such that
// (T,S) = ( VS*~S*WS^T, VS*~T*WS^T ), with ~S upper Schur and ~T upper triangular. The
// generalized Eigenvalues are returned in _CT_* ev[i]. As the lapack routine (XGGES) returns 
// ev=alpha/beta, beta may be found to be 0. In this case, we set *iflag=1 and ev[i]=0.
void SUBR(GenSchurDecomp)(_ST_* S, int ldS, _ST_* T, int ldT,
                          _ST_* VS, int ldVS, _ST_* WS, int ldWS,
                          int m, int nselect, int nsort, phist_EeigSort which, _MT_ tol,
                          void* v_ev, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  // this is for XGGES (computing the Schur form)
  int lwork = std::max(8*m+16,1024);
                          // min required workspace is 8*m+16
                          // for GGES 
  
  MT rwork[lwork];      // used as work in the real case
  int clwork=2*m;       // dimension of work in the complex case
  ST alpha[m], alphai[2*m], beta[m]; // in the complex case, alphai is used as WORK

  phist_blas_char jobvsl='V'; // compute the left Ritz vectors in VS
  phist_blas_char jobvsr='V'; // compute the right Ritz vectors in WR
  phist_blas_char sort='N';  // do not sort Ritz values (we do that later
                          // because gees only accepts the simple select
                          // function which does not compare the Ritz values)
  int sdim;
  CT* ev = (CT*)v_ev;
//  bool some_beta_zero=false;

  // can select at most m Ritz values (and at least 0)
  nselect=std::max(0,std::min(nselect,m));
  // can sort at most nselect Ritz values (and at least 0)
  nsort=std::max(0,std::min(nsort,nselect));

  PHIST_DEB("m=%d, nselect=%d, nsort=%d\n",m,nselect,nsort);

// prohibit parallel execution to assure identical results on different procs
#pragma omp parallel
  {
#pragma omp master
    {
#ifdef IS_COMPLEX
      PHIST_DEB("call complex %cGGES\n",st::type_char());
      *iflag=PHIST_TG_PREFIX(GGES)(SDMAT_FLAG,jobvsl,jobvsr,sort,NULL,m,(blas_cmplx*)S,ldS,
             (blas_cmplx*)T, ldT, &sdim,(blas_cmplx*)alpha, (blas_cmplx*)beta,
             (blas_cmplx*)VS,ldVS, (blas_cmplx*)WS, ldWS);
#else
      PHIST_DEB("call real %cGGES\n",st::type_char());
      *iflag=PHIST_TG_PREFIX(GGES)(SDMAT_FLAG,jobvsl,jobvsr,sort,NULL,m,S,ldS,
             T, ldT, &sdim, alpha, alphai, beta,
             VS,ldVS, WS, ldWS);

      for (int i=0;i<m;i++)
      {
        ev[i]=ct::zero();
        if (std::abs(beta[i])>mt::eps())
        {
          ev[i]=std::complex<MT>(alpha[i],alphai[i])/beta[i];
        }
//        else
//        {
//          some_beta_zero=true;
//        }
      }
#endif
    }
  }
  PHIST_CHK_IERR(;,*iflag);

#if PHIST_OUTLEV>=PHIST_DEBUG
//PHIST_OUT(0,"eigenvalues of unsorted Schur form:\n");
//for (int i=0;i<m;i++)
//{
  //PHIST_OUT(0,"%d\t%16.8g%+16.8gi\n",i,ct::real(ev[i]),ct::imag(ev[i]));
//}
#endif

  if (nselect<=0)
  { 
    *iflag=0;
    return;
  }
  if (nsort>nselect || nsort<0)
  {
    PHIST_OUT(PHIST_WARNING,"nselect=%d>=nsort=%d, or nsort>=0 "
                             "not satisfied, returning      unsorted Schur form\n",
                             nselect,nsort);
    *iflag=1;
    return;
  }

  //The function to sort a generalized Schur form in lapack is called XTGSEN

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
  MT pl, pr;// not used
  MT dif[2];
  int wantq=1,wantz=1;

  if (nselect<m)
  {
    PHIST_DEB("initial sort step, nselect=%d\n",nselect);
    // sort all eigenvalues according to 'which'.
    PHIST_CHK_IERR(SortEig(ev,m,idx,which,tol,iflag),*iflag);
    for (int i=0;i<nselect;i++) 
      select[std::abs(idx[i])]=1;
// prohibit parallel execution to assure identical results on different procs
#pragma omp parallel
    {
#pragma omp master
      {
        int ijob=0; // do not retrieve info on conditioning or deflating subspaces, maybe we could use it somehow?
#ifdef IS_COMPLEX
        int clwork=1;
        *iflag=PHIST_TG_PREFIX(TGSEN)(SDMAT_FLAG, ijob, wantq, wantz, select, m,
                (blas_cmplx*)S,ldS,(blas_cmplx*)T,ldT,(blas_cmplx*)alpha,(blas_cmplx*)beta,
                (blas_cmplx*)VS,ldVS,(blas_cmplx*)WS,ldWS,&m,
                &pl,&pr,dif);
#else
        *iflag=PHIST_TG_PREFIX(TGSEN)(SDMAT_FLAG, ijob, wantq, wantz, select, m,
                S,ldS,T,ldT,alpha,alphai,beta,
                VS,ldVS,WS,ldWS,&m,
                &pl,&pr,dif);


        for (int i=0;i<m;i++)
        {
          if (beta[i]!=st::zero())
          {
            ev[i]=std::complex<MT>(alpha[i],alphai[i])/beta[i];
          }
          else
          {
            ev[i]=st::zero();
          }
        }
#endif
      }
    }
    PHIST_CHK_IERR(;,*iflag);

    // TODO  - what about nsorted?
    PHIST_DEB("nsorted=%d\n",nsorted);
    // *POSSIBLE PROBLEM*
    // if we select part of a complex conjugate eigenpair, trsen just increases nselect by one
    if( nsorted > nselect )
    {
      PHIST_DEB("detected nsorted > nselect, try to continue with nsort <- nselect and nselect <- nsorted\n");
      nsort = nselect;
      nselect = nsorted;
    }
  }//nselect<m

  if (nselect==1) return; // the one (or two for complex pairs) selected eigenvalue
                          // according to 'which' is already 'sorted'

// prohibit parallel execution to assure identical results on different procs
#pragma omp parallel
  {
#pragma omp master
    {
      int i=0;
      while (i<nsort)
      {
        PHIST_DEB("sort step %d, nsorted=%d [%d]\n",i,nsorted,nsort);
        // sort the next few eigenvalues (up to nselect)
        SortEig(ev+i,nselect-i,idx+i,which,tol,iflag);
        if( *iflag != 0 )
          break;

        // sort next candidate to top.
        for (int j=0;j<i;j++) select[j]=1;
        for (int j=i;j<m;j++) select[j]=0;
        select[std::abs(idx[i])+i]=1; // sort next one to the top
        // note that the index returned by 
        // SortEig misses an offset i because
        // we pass in ev+i
        int nsorted_before=nsorted;
        int ijob=0; // this subroutine can do more, i.e. compute conditioning and
                    // deflating subspaces
#ifdef IS_COMPLEX
        *iflag=PHIST_TG_PREFIX(TGSEN)(SDMAT_FLAG, ijob, wantq, wantz, select, m,
                (blas_cmplx*)S,ldS,(blas_cmplx*)T,ldT,(blas_cmplx*)alpha,(blas_cmplx*)beta,
                (blas_cmplx*)VS,ldVS,(blas_cmplx*)WS,ldWS,&m,&pl,&pr,dif);
#else
        *iflag=PHIST_TG_PREFIX(TGSEN)(SDMAT_FLAG, ijob, wantq, wantz, select, m,
                S,ldS,T,ldT,alpha,alphai,beta,
                VS,ldVS,WS,ldWS,&m,
                &pl,&pr,dif);

        for (int j=0;j<m;j++)
        {
          if (beta[j]!=st::zero())
          {
            ev[j]=std::complex<MT>(alpha[j],alphai[j])/beta[j];
          }
          else
          {
            ev[j]=st::zero();
          }
        }
#endif
        if( *iflag != 0 )
          break;
        i+= std::max(nsorted-nsorted_before,1);
      }//while
#if PHIST_OUTLEV>=PHIST_DEBUG
      //PHIST_OUT(0,"eigenvalues of sorted Schur form:\n");
      //for (int i=0;i<m;i++)
      //{
        //PHIST_OUT(0,"%d\t%16.8g%+16.8gi\n",i,ct::real(ev[i]),ct::imag(ev[i]));
      //}
#endif
    }
  }
}

// reorder multiple eigenvalues in a given (partial) Schur decomposition by the smallest residual norm of the unprojected problem
void SUBR(ReorderPartialSchurDecomp)(_ST_* T, int ldT, _ST_* S, int ldS,
      int m, int nselected, phist_EeigSort which, _MT_ tol, _MT_* resNorm, void* v_ev, int* permutation, int *iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag = 0;
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN_SMALLDETERMINISTIC(ComputeTask)

  CT* ev = (CT*) v_ev;

  // set permutation to identity
  for(int i = 0; i < nselected; i++)
  {
    permutation[i] = i;
  }

#ifndef IS_COMPLEX
  // check for complex-conjugate eigenpairs
  for(int i = 0; i < nselected; i++)
  {
    if( mt::abs(ct::imag(ev[i])) > mt::eps() )
    {
      PHIST_SOUT(PHIST_ERROR,"reordering does not work for complex-conjugate eigenpairs in the real case!");
      PHIST_CHK_IERR(*iflag=PHIST_NOT_IMPLEMENTED,*iflag);
    }
  }
#endif

  // work array for lapack
#ifndef IS_COMPLEX
  _ST_ work[m];
#endif


// prohibit parallel execution to assure identical results on different procs
#pragma omp parallel
  {
#pragma omp master
    {

      // go through all eigenvalues and check if we need to do something
      // using gnome sort
      int pos = 1;
      while( pos < nselected )
      {
/*
#ifndef IS_COMPLEX
        // we cannot reorder conjugate complex eigenpairs currently
        if( ct::imag(ev[pos]) != 0 || ct::imag(ev[pos-1]) != 0 )
        {
          pos++;
          continue;
        }
#endif
*/
        // check if the next eigenvalue has same "magnitude"
        if( which == phist_LM && ct::abs(ev[pos]) < ct::abs(ev[pos-1])-tol )
        {
          pos++;
          continue;
        }
        if( which == phist_SM && ct::abs(ev[pos]) > ct::abs(ev[pos-1])+tol )
        {
          pos++;
          continue;
        }
        if( which == phist_LR && ct::real(ev[pos]) < ct::real(ev[pos-1])-tol )
        {
          pos++;
          continue;
        }
        if( which == phist_SR && ct::real(ev[pos]) > ct::real(ev[pos-1])+tol )
        {
          pos++;
          continue;
        }

        // we have two ev with the same "magnitude", check residuum
        if( which == phist_LM || which == phist_SM )
        {
          if( mt::abs(ct::abs(ev[pos])-ct::abs(ev[pos-1])) <= tol &&
              resNorm[pos] >= resNorm[pos-1] )
          {
            pos++;
            continue;
          }
        }
        if( which == phist_LR || which == phist_SR )
        {
          if( mt::abs(ct::real(ev[pos])-ct::real(ev[pos-1])) <= tol &&
              resNorm[pos] >= resNorm[pos-1] )
          {
            pos++;
            continue;
          }
        }


        // swap elements pos and pos-1
        std::swap(resNorm[pos],resNorm[pos-1]);
        std::swap(permutation[pos],permutation[pos-1]);
        std::swap(ev[pos],ev[pos-1]);

        phist_blas_char compq = 'V';
        int ifst = pos;
        int ilst = pos+1;
        PHIST_SOUT(PHIST_DEBUG,"swapping %d %d in unconverged eigenvalues\n",ifst-1,ilst-1);
#ifdef IS_COMPLEX
        *iflag=PHIST_TG_PREFIX(TREXC) (SDMAT_FLAG, compq, m, (blas_cmplx*) T,
            ldT, (blas_cmplx*) S, ldS, ifst, ilst);
#else
        *iflag=PHIST_TG_PREFIX(TREXC) (SDMAT_FLAG, compq, m, T, ldT, S, ldS,
            &ifst, &ilst);
#endif
        if( *iflag != 0 )
          break;
        PHIST_DEB("ifst = %d,\t ilst = %d\n", ifst-1, ilst-1);

        if( pos > 1 )
          pos--;
      }
    }
  }
PHIST_TASK_END(iflag)
}


 //! reorder multiple eigenvalues in a given (partial) generalized Schur decomposition by the smallest
 //! residual norm of the unprojected problem must be sorted up to nselected to work correctly!
 void SUBR(ReorderPartialGenSchurDecomp)(_ST_* S, int ldS, _ST_* T, int ldT, _ST_* VS, int ldVS, _ST_* WS, int ldWS,
           int m, int nselected, phist_EeigSort which, _MT_ tol, _MT_* resNorm, void* v_ev, int* permutation, int *iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag = 0;

  CT* ev = (CT*) v_ev;

  // set permutation to identity
  for(int i = 0; i < nselected; i++)
  {
    permutation[i] = i;
  }

#ifndef IS_COMPLEX
  // check for complex-conjugate eigenpairs
  for(int i = 0; i < nselected; i++)
  {
    if( mt::abs(ct::imag(ev[i])) > mt::eps() )
    {
      PHIST_SOUT(PHIST_ERROR,"reordering does not work for complex-conjugate eigenpairs in the real case!");
      PHIST_CHK_IERR(*iflag=PHIST_NOT_IMPLEMENTED,*iflag);
    }
  }
#endif

// prohibit parallel execution to assure identical results on different procs
#pragma omp parallel
  {
#pragma omp master
    {

      // go through all eigenvalues and check if we need to do something
      // using gnome sort
      int pos = 1;
      while( pos < nselected )
      {
/*
#ifndef IS_COMPLEX
        // we cannot reorder conjugate complex eigenpairs currently
        if( ct::imag(ev[pos]) != 0 || ct::imag(ev[pos-1]) != 0 )
        {
          pos++;
          continue;
        }
#endif
*/
        // check if the next eigenvalue has same "magnitude"
        if( which == phist_LM && ct::abs(ev[pos]) < ct::abs(ev[pos-1])-tol )
        {
          pos++;
          continue;
        }
        if( which == phist_SM && ct::abs(ev[pos]) > ct::abs(ev[pos-1])+tol )
        {
          pos++;
          continue;
        }
        if( which == phist_LR && ct::real(ev[pos]) < ct::real(ev[pos-1])-tol )
        {
          pos++;
          continue;
        }
        if( which == phist_SR && ct::real(ev[pos]) > ct::real(ev[pos-1])+tol )
        {
          pos++;
          continue;
        }

        // we have two ev with the same "magnitude", check residuum
        if( which == phist_LM || which == phist_SM )
        {
          if( mt::abs(ct::abs(ev[pos])-ct::abs(ev[pos-1])) <= tol &&
              resNorm[pos] >= resNorm[pos-1] )
          {
            pos++;
            continue;
          }
        }
        if( which == phist_LR || which == phist_SR )
        {
          if( mt::abs(ct::real(ev[pos])-ct::real(ev[pos-1])) <= tol &&
              resNorm[pos] >= resNorm[pos-1] )
          {
            pos++;
            continue;
          }
        }


        // swap elements pos and pos-1
        std::swap(resNorm[pos],resNorm[pos-1]);
        std::swap(permutation[pos],permutation[pos-1]);
        std::swap(ev[pos],ev[pos-1]);

        int wantq=1, wantz=1;
        int ifst = pos;
        int ilst = pos+1;
        PHIST_SOUT(PHIST_DEBUG,"swapping %d %d in unconverged eigenvalues\n",ifst-1,ilst-1);
#ifdef IS_COMPLEX
        *iflag=PHIST_TG_PREFIX(TGEXC) (SDMAT_FLAG, wantq, wantz, m,
                (blas_cmplx*) S, ldS, (blas_cmplx*) T, ldT,
                (blas_cmplx*) VS, ldVS, (blas_cmplx*) WS, ldWS,
                ifst, ilst);
#else
        *iflag=PHIST_TG_PREFIX(TGEXC) (SDMAT_FLAG, wantq, wantz, m,
                S, ldS, T, ldT, VS, ldVS, WS, ldWS,
                &ifst, &ilst);
#endif
        if( *iflag != 0 )
          break;
        PHIST_DEB("ifst = %d,\t ilst = %d\n", ifst-1, ilst-1);

        if( pos > 1 )
          pos--;
      }
    }
  }
}


// high-level function to convert an eigenbasis Q into eigenvectors X

// given Q,R: AQ=QR, computes X: AX=X*D with D a diagonal matrix.
// R is overwritten by D, X must be pre-allocated, Q must be orthonormal.
void SUBR(ComputeEigenvectors)(TYPE(const_mvec_ptr) Q, TYPE(sdMat_ptr) R,
                               TYPE(mvec_ptr) X, int* iflag)
{
#include "phist_std_typedefs.hpp"
    PHIST_ENTER_FCN(__FUNCTION__);
    *iflag=0;
    // compute eigenvectors. Given the Schur form R, first compute all its 
    // eigenvectors in S. Then compute the eigenvectors of A as Q*S.
    TYPE(sdMat_ptr) S, S0=NULL; // S0 views the first nX columns of S
    _ST_ *S_raw, *R_raw;
    int n, nX;
    phist_lidx ldR,ldS;
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(Q,&n,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(X,&nX,iflag),*iflag);
    PHIST_CHK_IERR(*iflag=(nX<=n)?0:PHIST_INVALID_INPUT,*iflag);
    phist_const_comm_ptr comm=NULL;
    PHIST_CHK_IERR(SUBR(mvec_get_comm)(Q,&comm,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_create)(&S,n,n,comm,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(S,&S0,0,n-1,0,nX-1,iflag),*iflag);
    phist::SdMatOwner<ST> _S(S),_S0(S0);

    PHIST_CHK_IERR(SUBR(sdMat_extract_view)(S,&S_raw, &ldS,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_extract_view)(R,&R_raw, &ldR,iflag),*iflag);

    PHIST_CHK_IERR(SUBR(sdMat_from_device)(R,iflag),*iflag);
    phist_blas_char side='R';
    phist_blas_char howmny='A';
    int m_out;
#ifdef IS_COMPLEX
    PHIST_CHK_IERR(*iflag=PHIST_TG_PREFIX(TREVC)(SDMAT_FLAG, side, howmny, NULL,
    n, (mt::blas_cmplx_t*)R_raw, ldR, NULL, ldS,
    (mt::blas_cmplx_t*)S_raw, ldS, n, &m_out), *iflag);
#else
    PHIST_CHK_IERR(*iflag=PHIST_TG_PREFIX(TREVC)(SDMAT_FLAG, side, howmny, NULL,
    n,R_raw, ldR, NULL, ldS, S_raw, ldS, n, &m_out),*iflag);
#endif

    PHIST_CHK_IERR(SUBR(sdMat_to_device)(R,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_to_device)(S,iflag),*iflag);
    
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(st::one(),Q,S0,st::zero(),X,iflag),*iflag);

}
