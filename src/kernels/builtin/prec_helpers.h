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
#include <math.h>
#include <stdlib.h>

#ifdef __cplusplus
#define restrict
#else
#define bool _Bool
#endif

static inline bool is_aligned(const void *restrict pointer, size_t byte_count)
{
  return (uintptr_t)pointer % byte_count == 0;
}


// Fused multiply add operation for doubles
static inline double double_fmadd(double a, double b, double c)
{
  __m128d a_ = _mm_set_sd(a);
  __m128d b_ = _mm_set_sd(b);
  __m128d c_ = _mm_set_sd(c);
  __m128d d_ = _mm_fmadd_pd(a_,b_,c_);
  double d;
  _mm_store_sd(&d,d_);
  return d;
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
  __m128d s2_a = _mm_set_sd(a); \
  __m128d s2_b = _mm_set_sd(b); \
  __m128d s2_s, s2_t;\
  MM128_FAST2SUM(s2_a,s2_b,s2_s,s2_t);\
  _mm_store_sd(&s,s2_s);\
  _mm_store_sd(&t,s2_t);\
}

// 2Sum
// s = round(a+b)
// s+t = a+b (exact)
#define MM256_2SUM(a,b,s,t) \
{\
  s = _mm256_add_pd(a,b);\
  __m256d s2_a_ = _mm256_sub_pd(s,b);\
  __m256d s2_b_ = _mm256_sub_pd(s,s2_a_);\
  __m256d s2_da = _mm256_sub_pd(a,s2_a_);\
  __m256d s2_db = _mm256_sub_pd(b,s2_b_);\
  t = _mm256_add_pd(s2_da,s2_db);\
}
#define MM128_2SUM(a,b,s,t) \
{\
  s = _mm_add_pd(a,b);\
  __m128d s2_a_ = _mm_sub_pd(s,b);\
  __m128d s2_b_ = _mm_sub_pd(s,s2_a_);\
  __m128d s2_da = _mm_sub_pd(a,s2_a_);\
  __m128d s2_db = _mm_sub_pd(b,s2_b_);\
  t = _mm_add_pd(s2_da,s2_db);\
}
#define DOUBLE_2SUM(a,b,s,t) \
{\
  __m128d s2_a = _mm_set_sd(a); \
  __m128d s2_b = _mm_set_sd(b); \
  __m128d s2_s, s2_t;\
  MM128_2SUM(s2_a,s2_b,s2_s,s2_t);\
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
  __m128d m2_s, m2_t;\
  MM128_2MULTFMA(m2_a,m2_b,m2_s,m2_t);\
  _mm_store_sd(&s,m2_s);\
  _mm_store_sd(&t,m2_t);\
}

// Precise addition of two high-precision numbers (a+aC) + (b+bC) of same exponent (or |a| >= |b|)
#define MM256_FAST4SUM(a,aC,b,bC,s,t)\
{\
  __m256d s4_s, s4_t;\
  MM256_FAST2SUM(a,b,s4_s,s4_t);\
  __m256d s4C_s, s4C_t;\
  MM256_FAST2SUM(aC,bC,s4C_s,s4C_t);\
  __m256d s4_st;\
  MM256_FAST2SUM(s4_s,s4C_s,s,s4_st);\
  __m256d s4_tt = _mm256_add_pd(s4_t,s4C_t);\
  t = _mm256_add_pd(s4_tt,s4_st);\
}
#define MM128_FAST4SUM(a,aC,b,bC,s,t)\
{\
  __m128d s4_s, s4_t;\
  MM128_FAST2SUM(a,b,s4_s,s4_t);\
  __m128d s4C_s, s4C_t;\
  MM128_FAST2SUM(aC,bC,s4C_s,s4C_t);\
  __m128d s4_st;\
  MM128_FAST2SUM(s4_s,s4C_s,s,s4_st);\
  __m128d s4_tt = _mm_add_pd(s4_t,s4C_t);\
  t = _mm_add_pd(s4_tt,s4_st);\
}
#define DOUBLE_FAST4SUM(a,aC,b,bC,s,t)\
{\
  double s4_s, s4_t;\
  DOUBLE_FAST2SUM(a,b,s4_s,s4_t);\
  double s4C_s, s4C_t;\
  DOUBLE_FAST2SUM(aC,bC,s4C_s,s4C_t);\
  double s4_st;\
  DOUBLE_FAST2SUM(s4_s,s4C_s,s,s4_st);\
  double s4_tt = s4_t+s4C_t;\
  t = s4_tt + s4_st;\
}

// Precise addition of two high-precision numbers (a+aC) + (b+bC)
#define MM256_4SUM(a,aC,b,bC,s,t)\
{\
  __m256d s4_s, s4_t;\
  MM256_2SUM(a,b,s4_s,s4_t);\
  __m256d s4C_s, s4C_t;\
  MM256_2SUM(aC,bC,s4C_s,s4C_t);\
  __m256d s4_st;\
  MM256_FAST2SUM(s4_s,s4C_s,s,s4_st);\
  __m256d s4_tt = _mm256_add_pd(s4_t,s4C_t);\
  t = _mm256_add_pd(s4_tt,s4_st);\
}
#define MM128_4SUM(a,aC,b,bC,s,t)\
{\
  __m128d s4_s, s4_t;\
  MM128_2SUM(a,b,s4_s,s4_t);\
  __m128d s4C_s, s4C_t;\
  MM128_2SUM(aC,bC,s4C_s,s4C_t);\
  __m128d s4_st;\
  MM128_FAST2SUM(s4_s,s4C_s,s,s4_st);\
  __m128d s4_tt = _mm_add_pd(s4_t,s4C_t);\
  t = _mm_add_pd(s4_tt,s4_st);\
}
#define DOUBLE_4SUM(a,aC,b,bC,s,t)\
{\
  double s4_s, s4_t;\
  DOUBLE_2SUM(a,b,s4_s,s4_t);\
  double s4C_s, s4C_t;\
  DOUBLE_2SUM(aC,bC,s4C_s,s4C_t);\
  double s4_st;\
  DOUBLE_FAST2SUM(s4_s,s4C_s,s,s4_st);\
  t = s4_t+s4C_t+s4_st;\
}

// Precise multiplication of two high-precision numbers (a+aC) * (b+bC)
#define MM256_4MULTFMA(a,aC,b,bC,s,t)\
{\
  __m256d m4ab_s, m4ab_t;\
  MM256_2MULTFMA(a,b,m4ab_s,m4ab_t);\
  __m256d m4abC_s, m4abC_t;\
  MM256_2MULTFMA(a,bC,m4abC_s,m4abC_t);\
  __m256d m4baC_s, m4baC_t;\
  MM256_2MULTFMA(b,aC,m4baC_s,m4baC_t);\
  __m256d m4C_s, m4C_st;\
  MM256_FAST2SUM(m4abC_s,m4baC_s,m4C_s,m4C_st);\
  __m256d m4_t, m4_tt;\
  MM256_FAST2SUM(m4ab_t,m4C_s,m4_t,m4_tt);\
  __m256d m4_st;\
  MM256_FAST2SUM(m4ab_s,m4_t,s,m4_st);\
  __m256d m4aCbC_t = _mm256_fmadd_pd(aC,bC,m4_st);\
  __m256d m4_ttt = _mm256_add_pd(m4_tt,m4C_st);\
  __m256d m4C_tt = _mm256_add_pd(m4abC_t,m4baC_t);\
  __m256d m4_tttt = _mm256_add_pd(m4aCbC_t,m4C_tt);\
  t =_mm256_add_pd(m4_tt,m4_tttt);\
}
#define MM128_4MULTFMA(a,aC,b,bC,s,t)\
{\
  __m128d m4ab_s, m4ab_t;\
  MM128_2MULTFMA(a,b,m4ab_s,m4ab_t);\
  __m128d m4abC_s, m4abC_t;\
  MM128_2MULTFMA(a,bC,m4abC_s,m4abC_t);\
  __m128d m4baC_s, m4baC_t;\
  MM128_2MULTFMA(b,aC,m4baC_s,m4baC_t);\
  __m128d m4C_s, m4C_st;\
  MM128_FAST2SUM(m4abC_s,m4baC_s,m4C_s,m4C_st);\
  __m128d m4_t, m4_tt;\
  MM128_FAST2SUM(m4ab_t,m4C_s,m4_t,m4_tt);\
  __m128d m4_st;\
  MM128_FAST2SUM(m4ab_s,m4_t,s,m4_st);\
  __m128d m4aCbC_t = _mm_fmadd_pd(aC,bC,m4_st);\
  __m128d m4_ttt = _mm_add_pd(m4_tt,m4C_st);\
  __m128d m4C_tt = _mm_add_pd(m4abC_t,m4baC_t);\
  __m128d m4_tttt = _mm_add_pd(m4aCbC_t,m4C_tt);\
  t =_mm_add_pd(m4_tt,m4_tttt);\
}
#define DOUBLE_4MULTFMA(a,aC,b,bC,s,t)\
{\
  double m4ab_s, m4ab_t;\
  DOUBLE_2MULTFMA(a,b,m4ab_s,m4ab_t);\
  double m4abC_s, m4abC_t;\
  DOUBLE_2MULTFMA(a,bC,m4abC_s,m4abC_t);\
  double m4baC_s, m4baC_t;\
  DOUBLE_2MULTFMA(b,aC,m4baC_s,m4baC_t);\
  double m4C_s, m4C_st;\
  DOUBLE_FAST2SUM(m4abC_s,m4baC_s,m4C_s,m4C_st);\
  double m4_t, m4_tt;\
  DOUBLE_FAST2SUM(m4ab_t,m4C_s,m4_t,m4_tt);\
  double m4_st;\
  DOUBLE_FAST2SUM(m4ab_s,m4_t,s,m4_st);\
  t = m4abC_t+m4baC_t+aC*bC+m4_st+m4_tt+m4C_st;\
}

// Div2FMA (almost) accurate division
// s = round(a/b)
// s+t = a/b
// (s+t)*b = s*b + t*b = rnd(a/b)*b + rnd(rnd(a-sb)/b)*b
#define MM256_2DIVFMA(a,b,s,t)\
{\
  s = _mm256_div_pd(a,b); \
  __m256d d2_bt = _mm256_fnmadd_pd(b,s,a); \
  t = _mm256_div_pd(d2_bt,b); \
}
#define MM128_2DIVFMA(a,b,s,t)\
{\
  s = _mm_div_pd(a,b); \
  __m128d d2_bt = _mm_fnmadd_pd(b,s,a); \
  t = _mm_div_pd(d2_bt,b); \
}
#define DOUBLE_2DIVFMA(a,b,s,t)\
{\
  __m128d d2_a = _mm_set_sd(a); \
  __m128d d2_b = _mm_set_sd(b); \
  __m128d d2_s, d2_t; \
  MM128_2DIVFMA(d2_a,d2_b,d2_s,d2_t); \
  _mm_store_sd(&s,d2_s);\
  _mm_store_sd(&t,d2_t);\
}


// Newton-Raphson iteration for precise calculation of 1/(a+aC)
#define MM256_4DIV_NEWTONRAPHSON_FMA(a,aC,div,divC)\
{\
  /* generate initial guess */ \
  __m256d d4_zero = _mm256_setzero_pd();\
  __m256d d4_one = _mm256_set1_pd(1.);\
  MM256_2DIVFMA(d4_one,a,div,divC);\
  /* use some iterations of Newton-Raphson (quadratic convergence) */ \
  for(int sq4_i = 0; sq4_i < 4; sq4_i++)\
  {\
    __m256d d4_e, d4_eC;\
    /* modified 1 - MM256_4MULTFMA(a,aC,div,divC,d4_adiv,d4_adivC) */ \
    {\
      __m256d m4ab_t = _mm256_fnmadd_pd(a,div,d4_one);\
      __m256d m4abC_s, m4abC_t;\
      MM256_2MULTFMA(a,divC,m4abC_s,m4abC_t);\
      __m256d m4baC_s, m4baC_t;\
      MM256_2MULTFMA(div,aC,m4baC_s,m4baC_t);\
      __m256d m4C_s, m4C_st;\
      MM256_FAST2SUM(m4abC_s,m4baC_s,m4C_s,m4C_st);\
      __m256d m4_t;\
      __m256d m4C_s_ = _mm256_sub_pd(d4_zero,m4C_s);\
      MM256_FAST2SUM(m4ab_t,m4C_s_,d4_e,m4_t);\
      __m256d m4_tt = _mm256_sub_pd(m4_t,m4C_st);\
      __m256d m4C_tt = _mm256_add_pd(m4abC_t,m4baC_t);\
      __m256d m4_ttt = _mm256_sub_pd(m4_tt,m4C_tt);\
      d4_eC = _mm256_fnmadd_pd(aC,divC,m4_ttt);\
    }\
    __m256d d4_dive, d4_diveC;\
    MM256_4MULTFMA(div,divC,d4_e,d4_eC,d4_dive,d4_diveC);\
    __m256d d4_div=div, d4_divC=divC;\
    MM256_FAST4SUM(d4_div,d4_divC,d4_dive,d4_diveC,div,divC);\
  }\
}
#define MM128_4DIV_NEWTONRAPHSON_FMA(a,aC,div,divC)\
{\
  /* generate initial guess */ \
  __m128d d4_zero = _mm_setzero_pd();\
  __m128d d4_one = _mm_set1_pd(1.);\
  MM128_2DIVFMA(d4_one,a,div,divC);\
  /* use some iterations of Newton-Raphson (quadratic convergence) */ \
  for(int sq4_i = 0; sq4_i < 4; sq4_i++)\
  {\
    __m128d d4_e, d4_eC;\
    /* modified 1 - MM128_4MULTFMA(a,aC,div,divC,d4_adiv,d4_adivC) */ \
    {\
      __m128d m4ab_t = _mm_fnmadd_pd(a,div,d4_one);\
      __m128d m4abC_s, m4abC_t;\
      MM128_2MULTFMA(a,divC,m4abC_s,m4abC_t);\
      __m128d m4baC_s, m4baC_t;\
      MM128_2MULTFMA(div,aC,m4baC_s,m4baC_t);\
      __m128d m4C_s, m4C_st;\
      MM128_FAST2SUM(m4abC_s,m4baC_s,m4C_s,m4C_st);\
      __m128d m4_t;\
      __m128d m4C_s_ = _mm_sub_pd(d4_zero,m4C_s);\
      MM128_FAST2SUM(m4ab_t,m4C_s_,d4_e,m4_t);\
      __m128d m4_tt = _mm_sub_pd(m4_t,m4C_st);\
      __m128d m4C_tt = _mm_add_pd(m4abC_t,m4baC_t);\
      __m128d m4_ttt = _mm_sub_pd(m4_tt,m4C_tt);\
      d4_eC = _mm_fnmadd_pd(aC,divC,m4_ttt);\
    }\
    __m128d d4_dive, d4_diveC;\
    MM128_4MULTFMA(div,divC,d4_e,d4_eC,d4_dive,d4_diveC);\
    __m128d d4_div=div, d4_divC=divC;\
    MM128_FAST4SUM(d4_div,d4_divC,d4_dive,d4_diveC,div,divC);\
  }\
}
#define DOUBLE_4DIV_NEWTONRAPHSON_FMA(a,aC,div,divC)\
{\
  __m128d dd4_a = _mm_set_sd(a);\
  __m128d dd4_aC = _mm_set_sd(aC);\
  __m128d dd4_div, dd4_divC;\
  MM128_4DIV_NEWTONRAPHSON_FMA(dd4_a,dd4_aC,dd4_div,dd4_divC);\
  _mm_store_sd(&div,dd4_div);\
  _mm_store_sd(&divC,dd4_divC);\
}


// Sqrt2FMA accurate square root
// s = round(sqrt(a))
// s+t = sqrt(a)
// (s+t)*(s+t) = s^2 + 2ts + t^2 = rnd(sqrt(a))^2 + 2*s*(rnd(a-s^2)/(2s)) + ...
#define MM256_2SQRTFMA(a,s,t)\
{\
  s = _mm256_sqrt_pd(a); \
  __m256d sq2_two = _mm256_set1_pd(2.);\
  __m256d sq2_2ts = _mm256_fnmadd_pd(s,s,a); \
  __m256d sq2_2s = _mm256_mul_pd(sq2_two,s);\
  t = _mm256_div_pd(sq2_2ts,sq2_2s); \
}
#define MM128_2SQRTFMA(a,s,t)\
{\
  s = _mm_sqrt_pd(a); \
  __m128d sq2_two = _mm_set1_pd(2.);\
  __m128d sq2_2ts = _mm_fnmadd_pd(s,s,a); \
  __m128d sq2_2s = _mm_mul_pd(sq2_two,s);\
  t = _mm_div_pd(sq2_2ts,sq2_2s); \
}
#define DOUBLE_2SQRTFMA(a,s,t)\
{\
  __m128d d2_a = _mm_set_sd(a);\
  __m128d d2_s, d2_t;\
  MM128_2SQRTFMA(d2_a,d2_s,d2_t);\
  _mm_store_sd(&s,d2_s);\
  _mm_store_sd(&t,d2_t);\
}


// Newton-Raphson iteration for precise calculation of sqrt(a+aC) and 1/sqrt(a+aC)
// modified variant to keep a+aC in the iteration (-> self-correcting property!)
// g -> sqrt(a)
// h -> 0.5/sqrt(a)
// iteration:
// r_(n+1) = 0.5 - g_n * h_n
// h_(n+1) = h_n + h_n*r_n
// g_(n+1) = 2*a*h_(n+1)        // originally g_n + g_n*r_n, but *a* is ignored then!
#define MM256_4SQRT_NEWTONRAPHSON_FMA(a,aC,sqrt_a,sqrt_aC,divsqrt_a,divsqrt_aC)\
{\
  /* helper variables */ \
  __m256d sq4_zero = _mm256_setzero_pd();\
  __m256d sq4_half = _mm256_set1_pd(0.5);\
  __m256d sq4_one = _mm256_set1_pd(1.);\
  __m256d sq4_two = _mm256_set1_pd(2.);\
  __m256d sq4_2a = _mm256_mul_pd(sq4_two,a);\
  __m256d sq4_2aC = _mm256_mul_pd(sq4_two,aC);\
  /* generate initial guess */ \
  __m256d sq4_g = _mm256_sqrt_pd(a);\
  __m256d sq4_gC = _mm256_setzero_pd();\
  __m256d sq4_h = _mm256_div_pd(sq4_half,sq4_g);\
  __m256d sq4_hC = _mm256_setzero_pd();\
  /* use some iterations of Newton-Raphson (quadratic convergence) */ \
  for(int sq4_i = 0; sq4_i < 4; sq4_i++)\
  {\
    __m256d sq4_r, sq4_rC;\
    /* modified 0.5 - MM256_4MULTFMA(sq4_g,sq4_gC,sq4_h,sq4_hC,...) */ \
    {\
      __m256d m4ab_t = _mm256_fnmadd_pd(sq4_g,sq4_h,sq4_half);\
      __m256d m4abC_s, m4abC_t;\
      MM256_2MULTFMA(sq4_g,sq4_hC,m4abC_s,m4abC_t);\
      __m256d m4baC_s, m4baC_t;\
      MM256_2MULTFMA(sq4_h,sq4_gC,m4baC_s,m4baC_t);\
      __m256d m4C_s, m4C_st;\
      MM256_FAST2SUM(m4abC_s,m4baC_s,m4C_s,m4C_st);\
      __m256d m4_t;\
      __m256d m4C_s_ = _mm256_sub_pd(sq4_zero,m4C_s);\
      MM256_FAST2SUM(m4ab_t,m4C_s_,sq4_r,m4_t);\
      __m256d m4_tt = _mm256_sub_pd(m4_t,m4C_st);\
      __m256d m4C_tt = _mm256_add_pd(m4abC_t,m4baC_t);\
      __m256d m4_ttt = _mm256_sub_pd(m4_tt,m4C_tt);\
      sq4_rC = _mm256_fnmadd_pd(sq4_gC,sq4_hC,m4_ttt);\
    }\
    __m256d sq4_hr, sq4_hrC;\
    MM256_4MULTFMA(sq4_h,sq4_hC,sq4_r,sq4_rC,sq4_hr,sq4_hrC);\
    __m256d sq4_h_ = sq4_h, sq4_hC_ = sq4_hC;\
    MM256_FAST4SUM(sq4_h_,sq4_hC_,sq4_hr,sq4_hrC,sq4_h,sq4_hC);\
    MM256_4MULTFMA(sq4_2a,sq4_2aC,sq4_h,sq4_hC,sq4_g,sq4_gC);\
  }\
  sqrt_a = sq4_g;\
  sqrt_aC = sq4_gC;\
  MM256_4MULTFMA(sq4_two,sq4_zero,sq4_h,sq4_hC,divsqrt_a,divsqrt_aC);\
}
#define MM128_4SQRT_NEWTONRAPHSON_FMA(a,aC,sqrt_a,sqrt_aC,divsqrt_a,divsqrt_aC)\
{\
  /* helper variables */ \
  __m128d sq4_zero = _mm_setzero_pd();\
  __m128d sq4_half = _mm_set1_pd(0.5);\
  __m128d sq4_one = _mm_set1_pd(1.);\
  __m128d sq4_two = _mm_set1_pd(2.);\
  __m128d sq4_2a = _mm_mul_pd(sq4_two,a);\
  __m128d sq4_2aC = _mm_mul_pd(sq4_two,aC);\
  /* generate initial guess */ \
  __m128d sq4_g = _mm_sqrt_pd(a);\
  __m128d sq4_gC = _mm_setzero_pd();\
  __m128d sq4_h = _mm_div_pd(sq4_half,sq4_g);\
  __m128d sq4_hC = _mm_setzero_pd();\
  /* use some iterations of Newton-Raphson (quadratic convergence) */ \
  for(int sq4_i = 0; sq4_i < 4; sq4_i++)\
  {\
    __m128d sq4_r, sq4_rC;\
    /* modified 0.5 - MM128_4MULTFMA(sq4_g,sq4_gC,sq4_h,sq4_hC,...) */ \
    {\
      __m128d m4ab_t = _mm_fnmadd_pd(sq4_g,sq4_h,sq4_half);\
      __m128d m4abC_s, m4abC_t;\
      MM128_2MULTFMA(sq4_g,sq4_hC,m4abC_s,m4abC_t);\
      __m128d m4baC_s, m4baC_t;\
      MM128_2MULTFMA(sq4_h,sq4_gC,m4baC_s,m4baC_t);\
      __m128d m4C_s, m4C_st;\
      MM128_FAST2SUM(m4abC_s,m4baC_s,m4C_s,m4C_st);\
      __m128d m4_t;\
      __m128d m4C_s_ = _mm_sub_pd(sq4_zero,m4C_s);\
      MM128_FAST2SUM(m4ab_t,m4C_s_,sq4_r,m4_t);\
      __m128d m4_tt = _mm_sub_pd(m4_t,m4C_st);\
      __m128d m4C_tt = _mm_add_pd(m4abC_t,m4baC_t);\
      __m128d m4_ttt = _mm_sub_pd(m4_tt,m4C_tt);\
      sq4_rC = _mm_fnmadd_pd(sq4_gC,sq4_hC,m4_ttt);\
    }\
    __m128d sq4_hr, sq4_hrC;\
    MM128_4MULTFMA(sq4_h,sq4_hC,sq4_r,sq4_rC,sq4_hr,sq4_hrC);\
    __m128d sq4_h_ = sq4_h, sq4_hC_ = sq4_hC;\
    MM128_FAST4SUM(sq4_h_,sq4_hC_,sq4_hr,sq4_hrC,sq4_h,sq4_hC);\
    MM128_4MULTFMA(sq4_2a,sq4_2aC,sq4_h,sq4_hC,sq4_g,sq4_gC);\
  }\
  sqrt_a = sq4_g;\
  sqrt_aC = sq4_gC;\
  MM128_4MULTFMA(sq4_two,sq4_zero,sq4_h,sq4_hC,divsqrt_a,divsqrt_aC);\
}
#define DOUBLE_4SQRT_NEWTONRAPHSON_FMA(a,aC,sqrt_a,sqrt_aC,divsqrt_a,divsqrt_aC)\
{\
  __m128d dsq4_a = _mm_set_sd(a);\
  __m128d dsq4_aC = _mm_set_sd(aC);\
  __m128d dsq4_sqrt, dsq4_sqrtC;\
  __m128d dsq4_divsqrt, dsq4_divsqrtC;\
  MM128_4SQRT_NEWTONRAPHSON_FMA(dsq4_a,dsq4_aC,dsq4_sqrt,dsq4_sqrtC,dsq4_divsqrt,dsq4_divsqrtC);\
  _mm_store_sd(&sqrt_a,dsq4_sqrt);\
  _mm_store_sd(&sqrt_aC,dsq4_sqrtC);\
  _mm_store_sd(&divsqrt_a,dsq4_divsqrt);\
  _mm_store_sd(&divsqrt_aC,dsq4_divsqrtC);\
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
