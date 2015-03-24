/*! \file axpy_kernels_prec.c
 * Common routines for more precise blas-like kernels.
 * \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
 *
*/

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#ifdef PHIST_HAVE_OPENMP
#include <omp.h>
#endif
#include "prec_helpers.h"

// b <- alpha*a + beta*b more precise size 1
void daxpby_prec_1(double alpha, const double *restrict a, double beta, double *restrict b, double *restrict bC)
{
  printf("Entering %s\n", __FUNCTION__);
  // we use AVX code here hoping the compiler doesn't optimize it away this way
  __m128d alpha_ = _mm_set_sd(alpha);
  __m128d beta_ = _mm_set_sd(beta);
  __m128d s, c, p, pi;

  __m128d ai = _mm_load_sd(a);
  __m128d bi = _mm_load_sd(b);
  MM128_2MULTFMA(beta_,bi, s,c);
  MM128_2MULTFMA(alpha_,ai, p, pi)
  __m128d oldS = s, sigma;
  MM128_FAST2SUM(oldS,p, s, sigma);
  __m128d tmp = _mm_add_pd(pi,sigma);
  c = _mm_add_pd(c,tmp);

  _mm_store_sd(b,s);
  _mm_store_sd(bC,c);
}

// b <- alpha*a + beta*b more precise size 2
void daxpby_prec_2(double alpha, const double *restrict a, double beta, double *restrict b, double *restrict bC)
{
  printf("Entering %s\n", __FUNCTION__);
  // we use AVX code here hoping the compiler doesn't optimize it away this way
  __m128d alpha_ = _mm_set1_pd(alpha);
  __m128d beta_ = _mm_set1_pd(beta);
  __m128d s, c, p, pi;

  __m128d ai = _mm_loadu_pd(a);
  __m128d bi = _mm_loadu_pd(b);
  MM128_2MULTFMA(beta_,bi, s,c);
  MM128_2MULTFMA(alpha_,ai, p, pi)
  __m128d oldS = s, sigma;
  MM128_FAST2SUM(oldS,p, s, sigma);
  __m128d tmp = _mm_add_pd(pi,sigma);
  c = _mm_add_pd(c,tmp);

  _mm_storeu_pd(b,s);
  _mm_storeu_pd(bC,c);
}

// b <- alpha*a + beta*b more precise size 2
void daxpby_prec_2k(int k, double alpha, const double *restrict a, double beta, double *restrict b, double *restrict bC)
{
  printf("Entering %s\n", __FUNCTION__);
  if( k%2 != 0 )
  {
    printf("wrong dimensions: %i\n",k);
    exit(1);
    return;
  }

  // we use AVX code here hoping the compiler doesn't optimize it away this way
  __m128d alpha_ = _mm_set1_pd(alpha);
  __m128d beta_ = _mm_set1_pd(beta);
  for(int i = 0; i < k; i+=2)
  {
    __m128d s, c, p, pi;

    __m128d ai = _mm_loadu_pd(&a[i]);
    __m128d bi = _mm_loadu_pd(&b[i]);
    MM128_2MULTFMA(beta_,bi, s,c);
    MM128_2MULTFMA(alpha_,ai, p, pi)
    __m128d oldS = s, sigma;
    MM128_FAST2SUM(oldS,p, s, sigma);
    __m128d tmp = _mm_add_pd(pi,sigma);
    c = _mm_add_pd(c,tmp);

    _mm_storeu_pd(&b[i],s);
    _mm_storeu_pd(&bC[i],c);
  }
}

// b <- alpha*a + beta*b more precise size 4
void daxpby_prec_4(double alpha, const double *restrict a, double beta, double *restrict b, double *restrict bC)
{
  printf("Entering %s\n", __FUNCTION__);
  // we use AVX code here hoping the compiler doesn't optimize it away this way
  __m256d alpha_ = _mm256_set1_pd(alpha);
  __m256d beta_ = _mm256_set1_pd(beta);
  __m256d s, c, p, pi;

  __m256d ai = _mm256_loadu_pd(a);
  __m256d bi = _mm256_loadu_pd(b);
  MM256_2MULTFMA(beta_,bi, s,c);
  MM256_2MULTFMA(alpha_,ai, p, pi)
    __m256d oldS = s, sigma;
  MM256_FAST2SUM(oldS,p, s, sigma);
  __m256d tmp = _mm256_add_pd(pi,sigma);
  c = _mm256_add_pd(c,tmp);

  _mm256_storeu_pd(b,s);
  _mm256_storeu_pd(bC,c);
}

// b <- alpha*a + beta*b more precise size k (with k%4 == 0)
void daxpby_prec_4k(int k, double alpha, const double *restrict a, double beta, double *restrict b, double *restrict bC)
{
  printf("Entering %s\n", __FUNCTION__);
  if( k%4 != 0 )
  {
    printf("wrong dimensions: %i\n",k);
    exit(1);
    return;
  }

  // we use AVX code here hoping the compiler doesn't optimize it away this way
  __m256d alpha_ = _mm256_set1_pd(alpha);
  __m256d beta_ = _mm256_set1_pd(beta);
  for(int i = 0; i < k; i+=4)
  {
    __m256d s, c, p, pi;

    __m256d ai = _mm256_loadu_pd(&a[i]);
    __m256d bi = _mm256_loadu_pd(&b[i]);
    MM256_2MULTFMA(beta_,bi, s,c);
    MM256_2MULTFMA(alpha_,ai, p, pi)
      __m256d oldS = s, sigma;
    MM256_FAST2SUM(oldS,p, s, sigma);
    __m256d tmp = _mm256_add_pd(pi,sigma);
    c = _mm256_add_pd(c,tmp);

    _mm256_storeu_pd(&b[i],s);
    _mm256_storeu_pd(&bC[i],c);
  }
}

