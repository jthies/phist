/*! \file prec_kernels_d.c
 * Implementation of some serial LAPaCK like functions with high precision
 * \author "Melven Roehrig-Zoellner" <Melven.Roehrig-Zoellner@DLR.de>
 * \author "Jonas Thies" <Jonas.Thies@DLR.de>
 *
 * This implementation can be used by any kernel lib that provides
 * PHIST_HIGH_PRECISION_KERNELS and implements sdMat_extract_err.
*/

#include "phist_config.h"
#include <stdlib.h>
#include "prec_helpers.h"
#include "phist_macros.h"
#include "phist_typedefs.h"
#include "mlapack_wrapper.h"

#ifdef PHIST_SDMATS_ROWMAJOR
#error "functionality not implemented for row-major sdMats"
#endif

// threshold at which to call a matrix rank deficient
#ifdef SINGTOL
#undef SINGTOL
#endif
#define SINGTOL 1.0e-25

// calculates a possibly low rank approximation of a lower cholesky factor of an spd matrix
// higher-precision + pivoting + stable low rank approximation
void phist_Dprec_cholesky(double *__restrict__ a, double *__restrict__ aC, phist_lidx n, phist_lidx lda, phist_lidx *perm, int *rank, int* iflag)
{
#if defined(TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif

  // permutation
  int p[n];
  for(int i = 0; i < n; i++)
    p[i] = i;
  // constructed L
  double l[n*n], lC[n*n];
  for(int i = 0; i < n*n; i++)
  {
    l[i] = 0.;
    lC[i] = 0.;
  }
  // diagonal entries
  double d[n], dC[n], diagNorm = 0;
  for(int i = 0; i < n; i++)
  {
    DOUBLE_FAST2SUM(a[i*lda+i],aC[i*lda+i],d[i],dC[i]);
    diagNorm += d[i]*d[i];
  }
  if( diagNorm == 0 )
  {
    printf("Warning zero diagonal in %s\n", __FUNCTION__);
    diagNorm = 1.e-32;
  }

  // keep order of upper left part locked, if requested
  int lock = *rank;
  *rank = 0;
  while(*rank < n)
  {
    // check rank
    double err = 0;
    for(int i = *rank; i < n; i++)
      err = err + d[p[i]];
//printf("step %d, err %e\n", *rank, err);
    if( err < SINGTOL*diagNorm )
      break;

    int m = *rank;
    *rank = *rank + 1;
    // find next pivot
    {
      int i = m;
      if( i >= lock )
      {
        for(int j = m+1; j < n; j++)
          if( d[p[j]] > d[p[i]] )
            i = j;
      }
      // swap p[i] p[m]
      int tmp = p[i];
      p[i] = p[m];
      p[m] = tmp;
//printf("pivot %d, perm", i);
//for(int j = 0; j < n; j++)
//  printf(" %d",p[j]);
    }
//printf("\n");

    // l_m,p[m] = sqrt(d_p[m])
    double div_lmm,div_lCmm;
    DOUBLE_4SQRT_NEWTONRAPHSON_FMA(d[p[m]],dC[p[m]],l[p[m]*n+m],lC[p[m]*n+m],div_lmm,div_lCmm);
//printf("m=%d,p[m]=%d: d_pm = %e, l_m,pm = %e, div_lmm = %e\n", m, p[m], l[p[m]*n+m], div_lmm);

    for(int i = m+1; i < n; i++)
    {
      // l_m,p[i] = 1/l_m,p[m] * ( a_p[m],p[i] - sum_j=0^m-2 l_j,p[m]*l_j,p[i] )
      {
        double s = -a[p[i]*lda+p[m]];
        double t = -aC[p[i]*lda+p[m]];
        for(int j = 0; j < m; j++)
        {
          double lj, ljC;
          DOUBLE_4MULTFMA(l[p[m]*n+j],lC[p[m]*n+j],l[p[i]*n+j],lC[p[i]*n+j],lj,ljC);
//printf("subtract ljm*lji %e\n", lj);
          double oldS = s, oldT = t;
          DOUBLE_4SUM(oldS,oldT,lj,ljC,s,t);
        }
        double oldS = s, oldT = t;
        DOUBLE_4MULTFMA(div_lmm,div_lCmm,oldS,oldT,s,t);
        l[p[i]*n+m] = -s;
        lC[p[i]*n+m] = -t;
//printf("l[m=%d,p[i=%d]] <- %e\n", m,i,-s_);
      }
      // d_p[i] = d_p[i]-l_m,p[i]^2
      {
        double s = -d[p[i]];
        double t = -dC[p[i]];
        double lmi = l[p[i]*n+m];
        double lCmi = lC[p[i]*n+m];
        double lm, lmC;
        DOUBLE_4MULTFMA(lmi,lCmi,lmi,lCmi,lm,lmC);
//printf("lm %e\n", lm);
        double oldS = s, oldT = t;
        DOUBLE_4SUM(oldS,oldT,lm,lmC,s,t);
        d[p[i]] = -s;
        dC[p[i]] = -t;
//printf("d[p[i=%d]] <- %e\n", i,-s);
      }
    }
  }

  // store result in a
  for(int i = 0; i < n; i++)
  {
    perm[i] = p[i];
    for(int j = 0; j < n; j++)
    {
      a[i*lda+j] = l[i*n+j];
      aC[i*lda+j] = lC[i*n+j];
    }
  }
  *iflag=0;
}


// apply backward substitution with permuted upper triangular matrix to k vectors in col-major storage
void phist_Dprec_backwardSubst(const double *__restrict__ r, const double *__restrict__ rC, phist_lidx n, phist_lidx ldr, phist_lidx *p, int rank,
        double *__restrict__ x, double *__restrict__ xC, phist_lidx k, phist_lidx ldx, int* iflag)
{
#if defined(TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif

  for(int l = 0; l < k; l++)
  {
    double newXl[rank], newXCl[rank];
    for(int i = rank-1; i >= 0; i--)
    {
      // for j = i+1..n
      // x_i,l <- x_i,l - r_i,p[j]*x_p[j],l
      double s =  -x[l*ldx+i], t = -xC[l*ldx+i];
//printf("x_p[i=%d],l=%d: %e\n", i, l, -s);
      for(int j = i+1; j < rank; j++)
      {
        double rx, rxC;
        DOUBLE_4MULTFMA(r[p[j]*ldr+i],rC[p[j]*ldr+i],newXl[j],newXCl[j],rx,rxC);
//printf("r: %e, x: %e, r*x: %e\n", r[p[j]*ldr+i],x[l*ldx+p[j]],rx);
        double oldS = s, oldT = t;
        DOUBLE_4SUM(oldS,oldT,rx,rxC,s,t);
      }
//printf("x_p[i=%d],l=%d-sum...: %e\n", i, l, -s);

      // x_p[i] = x_p[i]/r_i,p[i]
      double div_ri, div_riC;
      DOUBLE_4DIV_NEWTONRAPHSON_FMA(r[p[i]*ldr+i],rC[p[i]*ldr+i],div_ri,div_riC);
      double oldS = s, oldT = t;
      DOUBLE_4MULTFMA(oldS,oldT,div_ri,div_riC,s,t);
      newXl[i] = -s;
      newXCl[i] = -t;
//printf("new x_p[i=%d],l=%d: %e\n", i, l, -s);
    }

    for(int i = 0; i < rank; i++)
    {
      x[l*ldx+p[i]] = newXl[i];
      xC[l*ldx+p[i]] = newXCl[i];
    }
    for(int i = rank; i < n; i++)
    {
      x[l*ldx+p[i]] = 0;
      xC[l*ldx+p[i]] = 0;
    }
  }
}

// apply forward substitution with permuted transposed upper triangular matrix
void phist_Dprec_forwardSubst(const double *__restrict__ r, const double *__restrict__ rC, phist_lidx n, phist_lidx ldr, phist_lidx *p, int rank,
        double *__restrict__ x, double *__restrict__ xC, phist_lidx k, phist_lidx ldx, int* iflag)
{
#if defined(TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif

  for(int l = 0; l < k; l++)
  {
    double newXl[rank], newXCl[rank];
    for(int i = 0; i < rank; i++)
    {
      // for j = i+1..n
      // x_p[i],l <- x_p[i],l - r_j,p[i]*x_p[j],l
      double s =  -x[l*ldx+p[i]], t = -xC[l*ldx+p[i]];
//printf("x_p[i=%d],l=%d: %e\n", i, l, -s);
      for(int j = 0; j < i; j++)
      {
        double rx, rxC;
        DOUBLE_4MULTFMA(r[p[i]*ldr+j],rC[p[i]*ldr+j],newXl[j],newXCl[j],rx,rxC);
//printf("r: %e, x: %e, r*x: %e\n", r[p[i]*ldr+j],x[l*ldr+p[j]],rx);
        double oldS = s, oldT = t;
        DOUBLE_4SUM(oldS,oldT,rx,rxC,s,t);
      }
//printf("x_p[i=%d],l=%d-sum...: %e\n", i, l, -s);

      // x_p[i] = x_p[i]/r_i,p[i]
      double div_ri, div_riC;
      DOUBLE_4DIV_NEWTONRAPHSON_FMA(r[p[i]*ldr+i],rC[p[i]*ldr+i],div_ri,div_riC);
      double oldS = s, oldT = t;
      DOUBLE_4MULTFMA(oldS,oldT,div_ri,div_riC,s,t);
      newXl[i] = -s;
      newXCl[i] = -t;
//printf("new x_p[i=%d],l=%d: %e\n", i, l, -s);
    }

    // unpermute result
    for(int i = 0; i < rank; i++)
    {
      x[l*ldx+i] = newXl[i];
      xC[l*ldx+i] = newXCl[i];
    }
    for(int i = rank; i < n; i++)
    {
      x[l*ldx+p[i]] = 0;
      xC[l*ldx+p[i]] = 0;
    }
  }
}

#include "mlapack_wrapper.h"

// given symmetric A=V'V, compute B s.t. Q=V*B is orthonormal. B is computed in-place,
// if A is found to be numerically rank deficient, the last n-*rank columns of B will be zeroed out
// s.t. Q has exactly rank *rank.
void phist_Dprec_qb(double *__restrict__ a, double *__restrict__ aC, 
                    double *__restrict__ bi, double *__restrict__ biC,
                    phist_lidx n, phist_lidx lda, int *rank, int* iflag)
{
#ifdef PHIST_HAVE_MPACK_QD
  // compute sqrt(diag(A)) and its inverse
  double d[n], dC[n], di[n], diC[n];
  for (int i=0; i<n; i++)
  {
    // prevent nans
    if( a[i*lda+i] < SINGTOL )
    {
      d[i] = dC[i] = 0.;
      di[i] = diC[i] = 0.;
    }
    else
    {
      DOUBLE_4SQRT_NEWTONRAPHSON_FMA(a[i*lda+i],aC[i*lda+i],d[i],dC[i],di[i],diC[i]);
    }
  }
  
  // scale input matrix from left and right => diagonal elements 1
  for (int j=0; j<n; j++)
    for (int i=0; i<n; i++)
    {
      double aij,aijC;
      DOUBLE_4MULTFMA(a[j*lda+i],aC[j*lda+i],di[i],diC[i],aij,aijC);
      DOUBLE_4MULTFMA(aij,aijC,di[j],diC[j],a[j*lda+i],aC[j*lda+i]);
    }

  // compute eigenpairs of the scaled matrix A^
  double w[n], wC[n];
  PHIST_CHK_IERR(phist_Drsyev(n,a,aC,lda,w,wC,iflag),*iflag);

  // determine rank of input matrix and 1/sqrt(w)
  double wi[n], wiC[n];
  *rank=(int)n;
  
  // set w=sqrt(w) and wi=1/sqrt(w)
  for(int i=0; i<n; i++)
  {
    PHIST_SOUT(PHIST_DEBUG,"SV[%d]=%25.16e%+8.4e\n",i,w[i],wC[i]);
    if (w[i]<SINGTOL)
    {
      (*rank)--;
      w[i]=0.0; wC[i]=0.0;
      wi[i]=0.0; wiC[i]=0.0;
      PHIST_SOUT(PHIST_DEBUG,"<-0\n");
    }
    else
    {
      double sqrt_w, sqrt_wC;
      DOUBLE_4SQRT_NEWTONRAPHSON_FMA(w[i],wC[i],sqrt_w,sqrt_wC,wi[i],wiC[i]);
      w[i]=sqrt_w;
      wC[i]=sqrt_wC;
    }
  }

  // scale the eigenvector matrix, B^ <- Dinv * B * Einv
  // As B is orthonormal, the inverse of B^ is given by
  // E*B'*D (returned in bi, biC)
  for(int i=0; i<n; i++)
  {
    for(int j=0;j<n;j++)
    {
      double aij=a[j*lda+i], tmp;
      double aijC=aC[j*lda+i], tmpC;
      DOUBLE_4MULTFMA(aij,aijC,di[i],diC[i],tmp,tmpC);
      DOUBLE_4MULTFMA(tmp,tmpC,wi[j],wiC[j],a[j*lda+i],aC[j*lda+i]);
      if (bi!=NULL)
      {
        DOUBLE_4MULTFMA(aij,aijC,d[i],dC[i],tmp,tmpC);
        DOUBLE_4MULTFMA(tmp,tmpC,w[j],wC[j],bi[i*lda+j],biC[i*lda+j]);
      }
    }
  }
  *iflag=0;
#else
  *iflag=-99;
#endif
}
