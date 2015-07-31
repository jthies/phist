/*! \file gemm_kernels_sC_prec.c
 * Fast parallel BLAS-gemm like functions with high precision for different blocksizes for mvec_module
 * \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
 *
*/

#include "phist_config.h"
#include "prec_helpers.h"


// b+eps_b <- alpha*(a+eps_a) + beta*(b+eps_b) more precise
void daxpby_prec(int n, double alpha, const double *restrict a, const double *restrict aC, double beta, double *restrict b, double *restrict bC)
{
#if defined(TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif

  for(int i = 0; i < n; i++)
  {
    // b_ <- beta*b
    double b_, bC_;
    DOUBLE_2MULTFMA(beta,b[i], b_, bC_);
    bC_ = double_fmadd(beta,bC[i],bC_);

    // a_ <- alpha*a
    double a_, aC_;
    DOUBLE_2MULTFMA(alpha,a[i], a_,aC_);
    aC_ = double_fmadd(alpha,aC[i],aC_);

    // newB <- a_ + b_
    DOUBLE_4SUM(a_,aC_,b_,bC_,b[i],bC[i]);
  }
}


// precise calculation of C+Cc <- alpha*(A+Ac)*(B+Bc) + beta*(C+Cc)
void dgemm_prec(int m, int n, int k, double alpha, const double *restrict a, const double *restrict aC,
                                                   const double *restrict b, const double *restrict bC,
                                     double beta,        double *restrict c,       double *restrict cC)
{
#if defined(TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif

  for(int j = 0; j < n; j++)
  {
    for(int i = 0; i < m; i++)
    {
      // c_ <- beta*b[j*m+i]
      double c_, cC_;
      DOUBLE_2MULTFMA(beta,c[j*m+i], c_,cC_);
      cC_ = double_fmadd(beta,cC[j*m+i],cC_);

      for(int l = 0; l < k; l++)
      {
        // a_ <- alpha*a[l*m+i]
        double a_, aC_;
        DOUBLE_2MULTFMA(alpha,a[l*m+i], a_,aC_);
        aC_ = double_fmadd(alpha,aC[l*m+i],aC_);

        // tmp <- a_*b[j*k+l]
        double tmp, tmpC;
        DOUBLE_4MULTFMA(a_,aC_,b[j*k+l],bC[j*k+l],tmp,tmpC);

        // c_ <- c_ + tmp
        double oldC = c_, oldCC = cC_;
        DOUBLE_4SUM(oldC,oldCC,tmp,tmpC,c_,cC_);
      }
      // store result
      c[j*m+i] = c_;
      cC[j*m+i] = cC_;
    }
  }
}


// calculates a possibly low rank approximation of a lower cholesky factor of an spd matrix
// higher-precision + pivoting + stable low rank approximation
void cholesky_prec(int n, double *restrict a, double *restrict aC, int *perm, int *rank)
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
    DOUBLE_FAST2SUM(a[i*n+i],aC[i*n+i],d[i],dC[i]);
    diagNorm += d[i]*d[i];
  }
  diagNorm = sqrt(diagNorm);

  *rank = 0;
  while(*rank < n)
  {
    // check rank
    double err = 0;
    for(int i = *rank; i < n; i++)
      err = err + d[p[i]];
//printf("step %d, err %e\n", *rank, err);
    if( err < 1.e-25*diagNorm )
      break;

    int m = *rank;
    *rank = *rank + 1;
    // find next pivot
    {
      int i = m;
      for(int j = m+1; j < n; j++)
        if( d[p[j]] > d[p[i]] )
          i = j;
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
        double s = -a[p[i]*n+p[m]];
        double t = -aC[p[i]*n+p[m]];
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
      a[i*n+j] = l[i*n+j];
      aC[i*n+j] = lC[i*n+j];
    }
  }
}


// apply backward substitution with permuted upper triangular matrix
void backward_subst_prec(int n, int k, const double *restrict r, const double *restrict rC, int *p, int rank, double *restrict x, double *restrict xC)
{
#if defined(TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif

  for(int l = 0; l < k; l++)
  {
    double newXl[n], newXCl[n];
    for(int i = rank-1; i >= 0; i--)
    {
      // for j = i+1..n
      // x_i,l <- x_i,l - r_i,p[j]*x_p[j],l
      double s =  -x[l*n+i], t = -xC[l*n+i];
//printf("x_p[i=%d],l=%d: %e\n", i, l, -s);
      for(int j = i+1; j < n; j++)
      {
        double rx, rxC;
        DOUBLE_4MULTFMA(r[p[j]*n+i],rC[p[j]*n+i],newXl[j],newXCl[j],rx,rxC);
//printf("r: %e, x: %e, r*x: %e\n", r[p[j]*n+i],x[l*n+p[j]],rx);
        double oldS = s, oldT = t;
        DOUBLE_4SUM(oldS,oldT,rx,rxC,s,t);
      }
//printf("x_p[i=%d],l=%d-sum...: %e\n", i, l, -s);

      // x_p[i] = x_p[i]/r_i,p[i]
      double div_ri, div_riC;
      DOUBLE_4DIV_NEWTONRAPHSON_FMA(r[p[i]*n+i],rC[p[i]*n+i],div_ri,div_riC);
      double oldS = s, oldT = t;
      DOUBLE_4MULTFMA(oldS,oldT,div_ri,div_riC,s,t);
      newXl[i] = -s;
      newXCl[i] = -t;
//printf("new x_p[i=%d],l=%d: %e\n", i, l, -s);
    }

    for(int i = 0; i < n; i++)
    {
      x[l*n+p[i]] = newXl[i];
      xC[l*n+p[i]] = newXCl[i];
    }
  }
}

// apply forward substitution with permuted transposed upper triangular matrix
void forward_subst_prec(int n, int k, const double *restrict r, const double *restrict rC, int *p, int rank, double *restrict x, double *restrict xC)
{
#if defined(TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif

  for(int l = 0; l < k; l++)
  {
    double newXl[n], newXCl[n];
    for(int i = 0; i < rank; i++)
    {
      // for j = i+1..n
      // x_p[i],l <- x_p[i],l - r_j,p[i]*x_p[j],l
      double s =  -x[l*n+p[i]], t = -xC[l*n+p[i]];
//printf("x_p[i=%d],l=%d: %e\n", i, l, -s);
      for(int j = 0; j < i; j++)
      {
        double rx, rxC;
        DOUBLE_4MULTFMA(r[p[i]*n+j],rC[p[i]*n+j],newXl[j],newXCl[j],rx,rxC);
//printf("r: %e, x: %e, r*x: %e\n", r[p[i]*n+j],x[l*n+p[j]],rx);
        double oldS = s, oldT = t;
        DOUBLE_4SUM(oldS,oldT,rx,rxC,s,t);
      }
//printf("x_p[i=%d],l=%d-sum...: %e\n", i, l, -s);

      // x_p[i] = x_p[i]/r_i,p[i]
      double div_ri, div_riC;
      DOUBLE_4DIV_NEWTONRAPHSON_FMA(r[p[i]*n+i],rC[p[i]*n+i],div_ri,div_riC);
      double oldS = s, oldT = t;
      DOUBLE_4MULTFMA(oldS,oldT,div_ri,div_riC,s,t);
      newXl[i] = -s;
      newXCl[i] = -t;
//printf("new x_p[i=%d],l=%d: %e\n", i, l, -s);
    }

    // unpermute result
    for(int i = 0; i < n; i++)
    {
      x[l*n+i] = newXl[i];
      xC[l*n+i] = newXCl[i];
    }
  }
}


// create permuted backward substition factor for faster application later
void extract_backward_subst_matrix(int n, double *restrict r, double *restrict rC, int *p, int rank, double *restrict r_, double *restrict rC_, double *restrict dinv, double *restrict dinvC)
{
#if defined(TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif


  // calculate inverse diagonal entries
  for(int i = 0; i < rank; i++)
  {
    DOUBLE_4DIV_NEWTONRAPHSON_FMA(r[p[i]*n+i],rC[p[i]*n+i],dinv[i],dinvC[i]);
  }
  for(int i = rank; i < n; i++)
  {
    dinv[i] = 0.;
    dinvC[i] = 0.;
  }

  // fill lower left part of r_ with zeros
  for(int i = 0; i < n; i++)
  {
    for(int j = 0; j <= i; j++)
    {
      r_[j*n+i] = 0.;
      rC_[j*n+i] = 0.;
    }
  }

  // fill upper right part of r_ with unpermuted columns of r*dinv
  for(int i = 0; i < n; i++)
  {
    for(int j = i+1; j < n; j++)
    {
      DOUBLE_4MULTFMA(dinv[j],dinvC[j],r[p[j]*n+i],rC[p[j]*n+i],r_[j*n+i],rC_[j*n+i]);
    }
  }
}

