/*! \file sdmat_kernels_prec.c
   implements potentially slow and serial, but highly accurate, kernels for sdMats
 * \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
 *
*/

#include "phist_config.h"
#include "prec_helpers.h"

#ifdef PHIST_SDMATS_ROWMAJOR
#error "functionality not implemented for row-major sdMats"
#endif

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


