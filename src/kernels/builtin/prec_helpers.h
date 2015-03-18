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


#endif
