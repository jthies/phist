/*! \file prec_helpers.c
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

// precise reduction of gathered MPI results of all processes for block size 1
void prec_reduction_1(int n, const double *restrict s_, const double *restrict c_, double *restrict r, double *restrict rC)
{
#if defined(TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  // we need to sum up s_, c_
  // we use AVX code here hoping the compiler doesn't optimize it away this way
  __m128d s = _mm_load_sd(&s_[0]);
  __m128d c = _mm_load_sd(&c_[0]);
  for(int i = 1; i < n; i++)
  {
    __m128d si = _mm_load_sd(&s_[i]);
    __m128d ci = _mm_load_sd(&c_[i]);
    __m128d sigma, oldS = s;
    MM128_2SUM(oldS, si, s, sigma);
    //MM128_2SUM(oldS, si, s, sigma);
    __m128d tmp = _mm_add_pd(ci,sigma);
    c = _mm_add_pd(c, tmp);
  }
  _mm_store_sd(r,s);
  _mm_store_sd(rC,c);
}


// precise reduction of gathered MPI results of all processes for block size 2
void prec_reduction_2(int n, const double *restrict s_, const double *restrict c_, double *restrict r, double *restrict rC)
{
#if defined(TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  // we need to sum up s_, c_
  // we use AVX code here hoping the compiler doesn't optimize it away this way
  __m128d s = _mm_loadu_pd(s_);
  __m128d c = _mm_loadu_pd(c_);
  for(int i = 1; i < n; i++)
  {
    __m128d si = _mm_loadu_pd(&s_[2*i]);
    __m128d ci = _mm_loadu_pd(&c_[2*i]);
    __m128d sigma, oldS = s;
    MM128_2SUM(oldS, si, s, sigma);
    //MM128_2SUM(oldS, si, s, sigma);
    __m128d tmp = _mm_add_pd(ci,sigma);
    c = _mm_add_pd(c, tmp);
  }
  _mm_storeu_pd(r,s);
  _mm_storeu_pd(rC,c);
}


// precise reduction of gathered MPI results of all processes for block size 4
void prec_reduction_4(int n, const double *restrict s_, const double *restrict c_, double *restrict r, double *restrict rC)
{
#if defined(TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  // we need to sum up s_, c_
  // we use AVX code here hoping the compiler doesn't optimize it away this way
  __m256d s = _mm256_loadu_pd(s_);
  __m256d c = _mm256_loadu_pd(c_);
  for(int i = 1; i < n; i++)
  {
    __m256d si = _mm256_loadu_pd(&s_[4*i]);
    __m256d ci = _mm256_loadu_pd(&c_[4*i]);
    __m256d sigma, oldS = s;
    MM256_2SUM(oldS, si, s, sigma);
    //MM256_2SUM(oldS, si, s, sigma);
    __m256d tmp = _mm256_add_pd(ci,sigma);
    c = _mm256_add_pd(c, tmp);
  }
  _mm256_storeu_pd(r,s);
  _mm256_storeu_pd(rC,c);
}


// precise reduction of gathered MPI results of all processes for larger data
void prec_reduction_4k(int n, int k, const double *restrict s_, const double *restrict c_, double *restrict r, double *restrict rC)
{
#if defined(TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  // we need to sum up s_, c_
  // we use AVX code here hoping the compiler doesn't optimize it away this way
  if( k % 4 != 0 )
  {
    printf("wrong dimensions: %i\n", k);
    exit(1);
    return;
  }

  __m256d s[k/4];
  __m256d c[k/4];
  for(int j = 0; j < k/4; j++)
  {
    s[j] = _mm256_loadu_pd(&s_[4*j]);
    c[j] = _mm256_loadu_pd(&c_[4*j]);
  }
  for(int i = 1; i < n; i++)
  {
    for(int j = 0; j < k/4; j++)
    {
      __m256d si = _mm256_loadu_pd(&s_[k*i+4*j]);
      __m256d ci = _mm256_loadu_pd(&c_[k*i+4*j]);
      __m256d sigma, oldS = s[j];
      MM256_2SUM(oldS, si, s[j], sigma);
      //MM256_2SUM(oldS, si, s, sigma);
      __m256d tmp = _mm256_add_pd(ci,sigma);
      c[j] = _mm256_add_pd(c[j], tmp);
    }
  }
  for(int j = 0; j < k/4; j++)
  {
    _mm256_storeu_pd(&r[4*j],s[j]);
    _mm256_storeu_pd(&rC[4*j],c[j]);
  }
}


// precise reduction of gathered MPI results of all processes for block size 2
void prec_reduction_2k(int n, int k, const double *restrict s_, const double *restrict c_, double *restrict r, double *restrict rC)
{
#if defined(TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  // we need to sum up s_, c_
  // we use AVX code here hoping the compiler doesn't optimize it away this way
  __m128d s[k/2];
  __m128d c[k/2];
  for(int j = 0; j < k/2; j++)
  {
    s[j] = _mm_loadu_pd(&s_[2*j]);
    c[j] = _mm_loadu_pd(&c_[2*j]);
  }
  for(int i = 1; i < n; i++)
  {
    for(int j = 0; j < k/2; j++)
    {
      __m128d si = _mm_loadu_pd(&s_[k*i+2*j]);
      __m128d ci = _mm_loadu_pd(&c_[k*i+2*j]);
      __m128d sigma, oldS = s[j];
      MM128_2SUM(oldS, si, s[j], sigma);
      //MM128_2SUM(oldS, si, s, sigma);
      __m128d tmp = _mm_add_pd(ci,sigma);
      c[j] = _mm_add_pd(c[j], tmp);
    }
  }
  for(int j = 0; j < k/2; j++)
  {
    _mm_storeu_pd(&r[2*j],s[j]);
    _mm_storeu_pd(&rC[2*j],c[j]);
  }
}

// precise reduction of gathered MPI results of all processes for block size k
void prec_reduction_k(int n, int k, const double *restrict s_, const double *restrict c_, double *restrict r, double *restrict rC)
{
#if defined(TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  // we need to sum up s_, c_
  // we use AVX code here hoping the compiler doesn't optimize it away this way
  for(int j = 0; j < k; j++)
  {
    r[j] = s_[j];
    rC[j] = c_[j];
  }
  for(int i = 1; i < n; i++)
  {
    for(int j = 0; j < k/2; j++)
    {
      double si = s_[k*i+j];
      double ci = c_[k*i+j];
      double sigma, oldS = r[j];
      DOUBLE_2SUM(oldS, si, r[j], sigma);
      double tmp = ci+sigma;
      rC[j] += tmp;
    }
  }
}

