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
    bC_ = (bC_+beta*bC[i]);

    // a_ <- alpha*a
    double a_, aC_;
    DOUBLE_2MULTFMA(alpha,a[i], a_,aC_);

    // newB <- a_ + b_
    double newB, newBC;
    DOUBLE_2SUM(a_, b_, newB, newBC);
    newBC = (newBC+aC_+bC_);

    // round result again
    DOUBLE_FAST2SUM(newB, newBC, b[i], bC[i]);
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
      cC_ = beta*cC[j*m+i]+cC_;

      for(int l = 0; l < k; l++)
      {
        // a_ <- alpha*a[l*m+i]
        double a_, aC_;
        DOUBLE_2MULTFMA(alpha,a[l*m+i], a_,aC_);
        aC_ = alpha*aC[l*m+i]+aC_;

        // tmp <- a_*b[j*k+l]
        double tmp, tmpC;
        DOUBLE_2MULTFMA(a_,b[j*k+l], tmp,tmpC);
        tmpC = a_*bC[j*k+l]+aC_*b[j*k+l]+tmpC;

        // c_ <- c_ + tmp
        double oldC = c_, oldCC_ = cC_;
        DOUBLE_2SUM(oldC,tmp, c_, cC_);
        cC_ = oldCC_ + cC_ + tmpC;
      }

      // round result
      DOUBLE_FAST2SUM(c_,cC_,c[j*m+i],cC[j*m+i]);
    }
  }
}


// calculates a possibly low rank approximation of a lower cholesky factor of an spd matrix
// higher-precision + pivoting + stable low rank approximation
void cholesky_prec(int n, double *restrict a, double *restrict aC, int *perm, int *rank)
{
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
  double d[n], dC[n];
  for(int i = 0; i < n; i++)
  {
    DOUBLE_FAST2SUM(a[i*n+i],aC[i*n+i],d[i],dC[i]);
  }

  *rank = 0;
  while(*rank < n)
  {
    // check rank
    double err = 0;
    for(int i = *rank; i < n; i++)
      err = err + d[p[i]];
//printf("step %d, err %e\n", *rank, err);
    if( err < 1.e-12 )
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
    // uses sqrt(x+y) = sqrt(x) + y/(2*sqrt(x)) + ...
    {
      double s,t,t_;
      DOUBLE_2SQRTFMA(d[p[m]],s,t);
      t = t + dC[p[m]]/(2*s) - dC[p[m]]*dC[p[m]]/(8*s*s*s);
      DOUBLE_FAST2SUM(s,t,l[p[m]*n+m],lC[p[m]*n+m]);
//printf("l[m=%d,p[m=%d]] <- %e\n", m,m,l[p[m]*n+m]);
    }

    for(int i = m+1; i < n; i++)
    {
      // l_m,p[i] = 1/l_m,p[m] * ( a_p[m],p[i] - sum_j=0^m-2 l_j,p[m]*l_j,p[i] )
      {
        double s = -a[p[i]*n+p[m]];
        double t = -aC[p[i]*n+p[m]];
        for(int j = 0; j < m; j++)
        {
          double ljm = l[p[m]*n+j];
          double lCjm = lC[p[m]*n+j];
          double lji = l[p[i]*n+j];
          double lCji = lC[p[i]*n+j];
          double lj, ljC;
          DOUBLE_2MULTFMA(ljm,lji,lj,ljC);
//printf("subtract ljm*lji %e\n", lj);
          ljC = ljC+ljm*lCji+lji*lCjm + lCjm*lCji;
          double oldS = s, t_;
          DOUBLE_2SUM(lj,oldS,s,t_);
          t = t + t_ + ljC;
        }
        double s_, t_;
        DOUBLE_FAST2SUM(s,t,s_,t_);
        // use a/(x+y) = a/x - a*y/x^2 + ... to calculate s_/l_m,p[m]
        DOUBLE_2DIVFMA(s_,l[p[m]*n+m],s,t);
        t = t - s*lC[p[m]*n+m]/l[p[m]*n+m]*(1-lC[p[m]*n+m]/l[p[m]*n+m])+t_/l[p[m]*n+m];
        DOUBLE_FAST2SUM(s,t,s_,t_);
        l[p[i]*n+m] = -s_;
        lC[p[i]*n+m] = -t_;
//printf("l[m=%d,p[i=%d]] <- %e\n", m,i,-s_);
      }
      // d_p[i] = d_p[i]-l_m,p[i]^2
      {
        double s = -d[p[i]];
        double t = -dC[p[i]];
        double lmi = l[p[i]*n+m];
        double lCmi = lC[p[i]*n+m];
        double lm, lmC;
        DOUBLE_2MULTFMA(lmi,lmi,lm,lmC);
//printf("lm %e\n", lm);
        lmC = lmC+2*lmi*lCmi+lCmi*lCmi;
        double oldS = s, t_;
        DOUBLE_2SUM(oldS,lm,s,t_);
        t = t + t_ + lmC;
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

