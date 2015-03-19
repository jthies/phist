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
  // we need to sum up s_, c_
  // we use AVX code here hoping the compiler doesn't optimize it away this way
  __m128d s = _mm_load_sd(&s_[0]);
  __m128d c = _mm_load_sd(&c_[0]);
  for(int i = 1; i < n; i++)
  {
    __m128d si = _mm_load_sd(&s_[i]);
    __m128d ci = _mm_load_sd(&c_[i]);
    __m128d sigma, oldS = s;
    MM128_FAST2SUM(oldS, si, s, sigma);
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
  // we need to sum up s_, c_
  // we use AVX code here hoping the compiler doesn't optimize it away this way
  __m128d s = _mm_loadu_pd(s_);
  __m128d c = _mm_loadu_pd(c_);
  for(int i = 1; i < n; i++)
  {
    __m128d si = _mm_loadu_pd(&s_[2*i]);
    __m128d ci = _mm_loadu_pd(&c_[2*i]);
    __m128d sigma, oldS = s;
    MM128_FAST2SUM(oldS, si, s, sigma);
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
  // we need to sum up s_, c_
  // we use AVX code here hoping the compiler doesn't optimize it away this way
  __m256d s = _mm256_loadu_pd(s_);
  __m256d c = _mm256_loadu_pd(c_);
  for(int i = 1; i < n; i++)
  {
    __m256d si = _mm256_loadu_pd(&s_[4*i]);
    __m256d ci = _mm256_loadu_pd(&c_[4*i]);
    __m256d sigma, oldS = s;
    MM256_FAST2SUM(oldS, si, s, sigma);
    //MM256_2SUM(oldS, si, s, sigma);
    __m256d tmp = _mm256_add_pd(ci,sigma);
    c = _mm256_add_pd(c, tmp);
  }
  _mm256_storeu_pd(r,s);
  _mm256_storeu_pd(rC,c);
}



