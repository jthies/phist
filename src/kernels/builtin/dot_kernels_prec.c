/*! \file dot_kernels_prec.c
 * Fast parallel BLAS-axpy like functions with high precision for different blocksizes for mvec_module
 * \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
 *
*/

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include <stdint.h>
#include <stdio.h>
#include <emmintrin.h>
#include <immintrin.h>
#include <stdlib.h>
#ifdef PHIST_HAVE_OPENMP
#include <omp.h>
#endif

static inline _Bool is_aligned(const void *restrict pointer, size_t byte_count)
{
  return (uintptr_t)pointer % byte_count == 0;
}


// Fast2Sum (Kahan) without checking |a|>=b
// s = round(a+b)
// s+t = a+b (exact) for |a|>=|b|
#define MM256_FAST2SUM(a,b,s,t) \
{\
  s = _mm256_add_pd(a,b);\
  __m256d z = _mm256_sub_pd(s,a);\
  t = _mm256_sub_pd(b,z);\
}
#define MM128_FAST2SUM(a,b,s,t) \
{\
  s = _mm_add_pd(a,b);\
  __m128d z = _mm_sub_pd(s,a);\
  t = _mm_sub_pd(b,z);\
}

// 2Sum
// s = round(a+b)
// s+t = a+b (exact) for |a|>=|b|
#define MM256_2SUM(a,b,s,t) \
{\
  s = _mm256_add_pd(a,b);\
  __m256d a_ = _mm256_sub_pd(s,b);\
  __m256d b_ = _mm256_sub_pd(s,a);\
  __m256d da = _mm256_sub_pd(a,a_);\
  __m256d db = _mm256_sub_pd(b,b_);\
  t = _mm256_sub_pd(da,db);\
}
#define MM128_2SUM(a,b,s,t) \
{\
  s = _mm_add_pd(a,b);\
  __m128d a_ = _mm_sub_pd(s,b);\
  __m128d b_ = _mm_sub_pd(s,a);\
  __m128d da = _mm_sub_pd(a,a_);\
  __m128d db = _mm_sub_pd(b,b_);\
  t = _mm_sub_pd(da,db);\
}

// Mult2FMA fast AVX FMA implimentation of 2Prod
// s = round(a*b)
// s+t = a*b (exact)
#define MM256_2MULTFMA(a,b,s,t)\
{\
  s = _mm256_mul_pd(a,b);\
  t = _mm256_fmsub_pd(a,b,s);\
}
#define MM128_2MULTFMA(a,b,s,t)\
{\
  s = _mm_mul_pd(a,b);\
  t = _mm_fmsub_pd(a,b,s);\
}



// more accurate dot product x'x AVX2 kernel
void ddot_self_prec_4_(int nrows, const double *restrict x, double *restrict res, double *restrict resC)
{
  if( !is_aligned(x,32) )
  {
    printf("not aligned %lx\n", (uintptr_t)(void*)x);
    exit(1);
    return;
  }

#ifdef PHIST_HAVE_OPENMP
  int nt = omp_get_max_threads();
#else
  int nt = 1;
#endif

  {
    // buffer for omp thread result + padding to prevent false sharing
    __m256d s_[8][nt];
    __m256d c_[8][nt];

#pragma omp parallel shared(s_,c_)
    {
      // initialize sum
      __m256d s = _mm256_set1_pd((double)0.);
      __m256d c = _mm256_set1_pd((double)0.);

#pragma omp for schedule(static)
      for(int i = 0; i < nrows; i++)
      {
        __m256d xi = _mm256_load_pd(&x[4*i]);
        __m256d p, pi;
        MM256_2MULTFMA(xi,xi,p,pi);
        __m256d sigma, oldS = s;
        MM256_FAST2SUM(oldS,p, s,sigma);
        //MM256_2SUM(oldS,p, s,sigma); // more FP ops than Kahan-style FAST2SUM, but exacter
        __m256d tmp = _mm256_add_pd(pi,sigma);
        c = _mm256_add_pd(c,tmp);
      }

      int it = omp_get_thread_num();
      s_[0][it] = s;
      c_[0][it] = c;
    }


    // handcoded omp reduction
    __m256d s = s_[0][0];
    __m256d c = c_[0][0];
    for(int i = 1; i < nt; i++)
    {
      __m256d sigma, oldS = s;
      MM256_FAST2SUM(oldS, s_[0][i], s, sigma);
      //MM256_2SUM(oldS, s_[0][i], s, sigma);
      __m256d tmp = _mm256_add_pd(c_[0][i],sigma);
      c = _mm256_add_pd(c, tmp);
    }

    // return result, needs to be summed up
    _mm256_storeu_pd(res,  s);
    _mm256_storeu_pd(resC, c);
  }

}


// more accurate dot product x'y AVX2 kernel
void ddot_prec_4_(int nrows, const double *restrict x, const double *restrict y, double *restrict res, double *restrict resC)
{
  if( !is_aligned(x,32) )
  {
    printf("not aligned %lx\n", (uintptr_t)(void*)x);
    exit(1);
    return;
  }

  if( !is_aligned(y,32) )
  {
    printf("not aligned %lx\n", (uintptr_t)(void*)y);
    exit(1);
    return;
  }

#ifdef PHIST_HAVE_OPENMP
  int nt = omp_get_max_threads();
#else
  int nt = 1;
#endif

  {
    // buffer for omp thread result + padding to prevent false sharing
    __m256d s_[8][nt];
    __m256d c_[8][nt];

#pragma omp parallel shared(s_,c_)
    {
      // initialize sum
      __m256d s = _mm256_set1_pd((double)0.);
      __m256d c = _mm256_set1_pd((double)0.);

#pragma omp for schedule(static)
      for(int i = 0; i < nrows; i++)
      {
        __m256d xi = _mm256_load_pd(&x[4*i]);
        __m256d yi = _mm256_load_pd(&y[4*i]);
        __m256d p, pi;
        MM256_2MULTFMA(xi,yi,p,pi);
        __m256d sigma, oldS = s;
        MM256_FAST2SUM(oldS,p, s,sigma);
        //MM256_2SUM(oldS,p, s,sigma); // more FP ops than Kahan-style FAST2SUM, but exacter
        __m256d tmp = _mm256_add_pd(pi,sigma);
        c = _mm256_add_pd(c,tmp);
      }

      int it = omp_get_thread_num();
      s_[0][it] = s;
      c_[0][it] = c;
    }


    // handcoded omp reduction
    __m256d s = s_[0][0];
    __m256d c = c_[0][0];
    for(int i = 1; i < nt; i++)
    {
      __m256d sigma, oldS = s;
      MM256_FAST2SUM(oldS, s_[0][i], s, sigma);
      //MM256_2SUM(oldS, s_[0][i], s, sigma);
      __m256d tmp = _mm256_add_pd(c_[0][i],sigma);
      c = _mm256_add_pd(c, tmp);
    }

    // return result, needs to be summed up
    _mm256_storeu_pd(res,  s);
    _mm256_storeu_pd(resC, c);
  }

}


// precise reduction of gathered MPI results of all processes for block size 1
void prec_reduction_1_(int n, const double *restrict s_, const double *restrict c_, double *restrict r, double *restrict rC)
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
void prec_reduction_2_(int n, const double *restrict s_, const double *restrict c_, double *restrict r, double *restrict rC)
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
void prec_reduction_4_(int n, const double *restrict s_, const double *restrict c_, double *restrict r, double *restrict rC)
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


// more accurate dot product x'x for block size 1
void ddot_self_prec_1_(int nrows, const double *restrict x, double *restrict res, double *restrict resC)
{
  double s[4], c[4];

  // assume appropriate padding with zeros!
  ddot_self_prec_4_(nrows/4, x, s, c);

  // we still need to sum up s, c
  prec_reduction_1_(4, s, c, res, resC);
}


// more accurate dot product x'y for block size 1
void ddot_prec_1_(int nrows, const double *restrict x, const double *restrict y, double *restrict res, double *restrict resC)
{
  double s[4], c[4];

  // assume appropriate padding with zeros!
  ddot_prec_4_(nrows/4, x, y, s, c);

  // we still need to sum up s, c
  prec_reduction_1_(4, s, c, res, resC);
}


// more accurate dot product x'x for block size 2
void ddot_self_prec_2_(int nrows, const double *restrict x, double *restrict res, double *restrict resC)
{
  double s[4], c[4];

  // assume appropriate padding with zeros!
  ddot_self_prec_4_(nrows/2, x, s, c);

  // we still need to sum up s, c
  prec_reduction_2_(2, s, c, res, resC);
}


// more accurate dot product x'y for block size 2
void ddot_prec_2_(int nrows, const double *restrict x, const double *restrict y, double *restrict res, double *restrict resC)
{
  double s[4], c[4];

  // assume appropriate padding with zeros!
  ddot_prec_4_(nrows/2, x, y, s, c);

  // we still need to sum up s, c
  prec_reduction_2_(2, s, c, res, resC);
}


