/*! \file prec_helpers.h
 * Some helpful definitions for fast parallel BLAS like functions with high precision for different blocksizes for mvec_module
 * \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
 *
*/

#ifndef PREC_HELPERS_H
#define PREC_HELPERS_H

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include <emmintrin.h>
#include <immintrin.h>
#include <stdint.h>
#include <stdio.h>

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
#define DOUBLE_FAST2SUM(a,b,s,t) \
{\
  __m128d s2_a = _mm_set_sd(a);\
  __m128d s2_b = _mm_set_sd(b);\
  __m128d s2_s = _mm_add_pd(s2_a,s2_b);\
  __m128d s2_z = _mm_sub_pd(s2_s,s2_a);\
  __m128d s2_t = _mm_sub_pd(s2_b,s2_z);\
  _mm_store_sd(&s,s2_s);\
  _mm_store_sd(&t,s2_t);\
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
#define DOUBLE_2SUM(a,b,s,t) \
{\
  __m128d s2_a = _mm_set_sd(a); \
  __m128d s2_b = _mm_set_sd(b); \
  __m128d s2_s = _mm_add_pd(s2_a,s2_b);\
  __m128d s2_a_ = _mm_sub_pd(s2_s,s2_b);\
  __m128d s2_b_ = _mm_sub_pd(s2_s,s2_a);\
  __m128d s2_da = _mm_sub_pd(s2_a,s2_a_);\
  __m128d s2_db = _mm_sub_pd(s2_b,s2_b_);\
  __m128d s2_t = _mm_sub_pd(s2_da,s2_db);\
  _mm_store_sd(&s,s2_s);\
  _mm_store_sd(&t,s2_t);\
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
#define DOUBLE_2MULTFMA(a,b,s,t)\
{\
  __m128d m2_a = _mm_set_sd(a); \
  __m128d m2_b = _mm_set_sd(b); \
  __m128d m2_s = _mm_mul_pd(m2_a,m2_b);\
  __m128d m2_t = _mm_fmsub_pd(m2_a,m2_b,m2_s);\
  _mm_store_sd(&s,m2_s);\
  _mm_store_sd(&t,m2_t);\
}

// Div2FMA accurate division
// s = round(a/b)
// s+t = a/b
// (melven: I hope this is correct)
#define MM256_2DIVFMA(a,b,s,t)\
{\
  s = _mm256_div_pd(a,b); \
  t = _mm256_fnmadd_pd(s,b,a); \
}
#define MM128_2DIVFMA(a,b,s,t)\
{\
  s = _mm_div_pd(a,b); \
  t = _mm_fnmadd_pd(s,b,a); \
}
#define DOUBLE_2DIVFMA(a,b,s,t)\
{\
  __m128d m2_a = _mm_set_sd(a); \
  __m128d m2_b = _mm_set_sd(b); \
  __m128d m2_s = _mm_div_pd(m2_a,m2_b); \
  __m128d m2_t = _mm_fnmadd_pd(m2_s,m2_b,m2_a); \
  _mm_store_sd(&s,m2_s);\
  _mm_store_sd(&t,m2_t);\
}

// Sqrt2FMA accurate square root
// s = round(sqrt(a))
// s+t = sqrt(a)
// (melven: I hope this is correct)
#define MM256_2SQRTFMA(a,s,t)\
{\
  s = _mm256_sqrt_pd(a); \
  t = _mm256_fnmadd_pd(s,s,a); \
}
#define M128_2SQRTFMA(a,s,t)\
{\
  s = _mm_sqrt_pd(a); \
  t = _mm_fnmadd_pd(s,s,a); \
}
#define DOUBLE_2SQRTFMA(a,s,t)\
{\
  __m128d s2_a = _mm_set_sd(a); \
  __m128d s2_s = _mm_sqrt_pd(s2_a); \
  __m128d s2_t = _mm_fnmadd_pd(s2_s,s2_s,s2_a); \
  _mm_store_sd(&s,s2_s);\
  _mm_store_sd(&t,s2_t);\
}


// precise reduction of gathered MPI results of all processes for block size 1
void prec_reduction_1(int n, const double *restrict s_, const double *restrict c_, double *restrict r, double *restrict rC);
// precise reduction of gathered MPI results of all processes for block size 2
void prec_reduction_2(int n, const double *restrict s_, const double *restrict c_, double *restrict r, double *restrict rC);
// precise reduction of gathered MPI results of all processes for block size 4
void prec_reduction_4(int n, const double *restrict s_, const double *restrict c_, double *restrict r, double *restrict rC);
// precise reduction of gathered MPI results of all processes for block size 2
void prec_reduction_2k(int n, int k, const double *restrict s_, const double *restrict c_, double *restrict r, double *restrict rC);
// precise reduction of gathered MPI results of all processes for block size 4
void prec_reduction_4k(int n, int k, const double *restrict s_, const double *restrict c_, double *restrict r, double *restrict rC);
// precise reduction of gathered MPI results of all processes for block size k
void prec_reduction_k(int n, int k, const double *restrict s_, const double *restrict c_, double *restrict r, double *restrict rC);


#endif
