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
/* don't really work, probably optimized away by the compiler!
#define MM256_2SUM(a,b,s,t) \
{\
  s = _mm256_add_pd(a,b);\
  __m256d s2_a_ = _mm256_sub_pd(s,b);\
  __m256d s2_b_ = _mm256_sub_pd(s,s2_a_);\
  __m256d s2_da = _mm256_sub_pd(a,s2_a_);\
  __m256d s2_db = _mm256_sub_pd(b,s2_b_);\
  t = _mm256_sub_pd(s2_da,s2_db);\
}
#define MM128_2SUM(a,b,s,t) \
{\
  s = _mm_add_pd(a,b);\
  __m128d s2_a_ = _mm_sub_pd(s,b);\
  __m128d s2_b_ = _mm_sub_pd(s,s2_a_);\
  __m128d s2_da = _mm_sub_pd(a,s2_a_);\
  __m128d s2_db = _mm_sub_pd(b,s2_b_);\
  t = _mm_sub_pd(s2_da,s2_db);\
}
*/
#define DOUBLE_2SUM(a,b,s,t) \
{\
  if( abs(a) > abs(b) ) \
  {\
    DOUBLE_FAST2SUM(a,b,s,t); \
  }\
  else \
  {\
    DOUBLE_FAST2SUM(b,a,s,t); \
  }\
}
/*
  __m128d s2_a = _mm_set_sd(a); \
  __m128d s2_b = _mm_set_sd(b); \
  __m128d s2_s, s2_t;\
  MM128_2SUM(s2_a,s2_b,s2_s,s2_t);\
  _mm_store_sd(&s,s2_s);\
  _mm_store_sd(&t,s2_t);\
}
*/

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

// Precise addition of two high-precision numbers (a+aC) + (b+bC) of similar size
#define MM256_FAST4SUM(a,aC,b,bC,s,t)\
{\
  __m256d s4_s, s4_t;\
  MM256_FAST2SUM(a,b,s4_s,s4_t);\
  __m256d s4_tt = _mm256_add_pd(aC,bC);\
  __m256d s4_ttt = _mm256_add_pd(s4_tt,s4_t);\
  MM256_FAST2SUM(s4_s,s4_ttt,s,t);\
}
#define MM128_FAST4SUM(a,aC,b,bC,s,t)\
{\
  __m128d s4_s, s4_t;\
  MM128_FAST2SUM(a,b,s4_s,s4_t);\
  __m128d s4_tt = _mm_add_pd(aC,bC);\
  __m128d s4_ttt = _mm_add_pd(s4_tt,s4_t);\
  MM128_FAST2SUM(s4_s,s4_ttt,s,t);\
}
#define DOUBLE_FAST4SUM(a,aC,b,bC,s,t)\
{\
  double s4_s, s4_t;\
  DOUBLE_FAST2SUM(a,b,s4_s,s4_t);\
  double s4_tt = (aC+bC)+s4_t;\
  DOUBLE_FAST2SUM(s4_s,s4_tt,s,t);\
}

// Precise addition of two high-precision numbers (a+aC) + (b+bC)
/*
#define MM256_4SUM(a,aC,b,bC,s,t)\
{\
  __m256d s4_s, s4_t;\
  MM256_2SUM(a,b,s4_s,s4_t);\
  __m256d s4_tt = _mm256_add_pd(aC,bC);\
  __m256d s4_ttt = _mm256_add_pd(s4_tt,s4_t);\
  MM256_FAST2SUM(s4_s,s4_ttt,s,t);\
}
#define MM128_4SUM(a,aC,b,bC,s,t)\
{\
  __m128d s4_s, s4_t;\
  MM128_2SUM(a,b,s4_s,s4_t);\
  __m128d s4_tt = _mm_add_pd(aC,bC);\
  __m128d s4_ttt = _mm_add_pd(s4_tt,s4_t);\
  MM128_FAST2SUM(s4_s,s4_ttt,s,t);\
}
*/
#define DOUBLE_4SUM(a,aC,b,bC,s,t)\
{\
  double s4_s, s4_t;\
  DOUBLE_2SUM(a,b,s4_s,s4_t);\
  double s4_tt = (aC+bC)+s4_t;\
  DOUBLE_FAST2SUM(s4_s,s4_tt,s,t);\
}

// Precise multiplication of two high-precision numbers (a+aC) * (b+bC)
#define MM256_4MULTFMA(a,aC,b,bC,s,t)\
{\
  __m256d m4_s, m4_t;\
  MM256_2MULTFMA(a,b,m4_s,m4_t);\
  __m256d m4_tt = _mm256_fmadd_pd(a,bC,m4_t);\
  __m256d m4_ttt = _mm256_fmadd_pd(b,aC,m4_tt);\
  __m256d m4_tttt = _mm256_fmadd_pd(aC,bC,m4_ttt);\
  MM256_FAST2SUM(m4_s,m4_tttt,s,t);\
}
#define MM128_4MULTFMA(a,aC,b,bC,s,t)\
{\
  __m128d m4_s, m4_t;\
  MM128_2MULTFMA(a,b,m4_s,m4_t);\
  __m128d m4_tt = _mm_fmadd_pd(a,bC,m4_t);\
  __m128d m4_ttt = _mm_fmadd_pd(b,aC,m4_tt);\
  __m128d m4_tttt = _mm_fmadd_pd(aC,bC,m4_ttt);\
  MM128_FAST2SUM(m4_s,m4_tttt,s,t);\
}
#define DOUBLE_4MULTFMA(a,aC,b,bC,s,t)\
{\
  double m4_s, m4_t;\
  DOUBLE_2MULTFMA(a,b,m4_s,m4_t);\
  double m4_tt = m4_t + a*bC + b*aC + aC*bC;\
  DOUBLE_FAST2SUM(m4_s,m4_tt,s,t);\
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


// Newton-Raphson iteration for precise calculation of 1/(a+aC)
#define MM256_4DIV_NEWTONRAPHSON_FMA(a,aC,div,divC)\
{\
  /* generate initial guess */ \
  __m256d d4_zero = _mm256_setzero_pd();\
  __m256d d4_one = _mm256_set1_pd(1.);\
  div = _mm256_div_pd(d4_one,a);\
  divC = _mm256_setzero_pd();\
  /* use some iterations of Newton-Raphson (quadratic convergence) */ \
  for(int sq4_i = 0; sq4_i < 4; sq4_i++)\
  {\
    __m256d d4_adiv, d4_adivC;\
    MM256_4MULTFMA(a,aC,div,divC,d4_adiv,d4_adivC);\
    d4_adiv=-d4_adiv; d4_adivC=-d4_adivC;\
    __m256d d4_e, d4_eC;\
    MM256_FAST4SUM(d4_one,d4_zero,d4_adiv,d4_adivC,d4_e,d4_eC);\
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
  div = _mm_div_pd(d4_one,a);\
  divC = _mm_setzero_pd();\
  /* use some iterations of Newton-Raphson (quadratic convergence) */ \
  for(int sq4_i = 0; sq4_i < 4; sq4_i++)\
  {\
    __m128d d4_adiv, d4_adivC;\
    MM128_4MULTFMA(a,aC,div,divC,d4_adiv,d4_adivC);\
    d4_adiv=-d4_adiv; d4_adivC=-d4_adivC;\
    __m128d d4_e, d4_eC;\
    MM128_FAST4SUM(d4_one,d4_zero,d4_adiv,d4_adivC,d4_e,d4_eC);\
    __m128d d4_dive, d4_diveC;\
    MM128_4MULTFMA(div,divC,d4_e,d4_eC,d4_dive,d4_diveC);\
    __m128d d4_div=div, d4_divC=divC;\
    MM128_FAST4SUM(d4_div,d4_divC,d4_dive,d4_diveC,div,divC);\
  }\
}
#define DOUBLE_4DIV_NEWTONRAPHSON_FMA(a,aC,div,divC)\
{\
  /* generate initial guess */ \
  double d4_zero = 0.;\
  double d4_one = 1.;\
  div = d4_one/a;\
  divC = 0.;\
  /* use some iterations of Newton-Raphson (quadratic convergence) */ \
  for(int sq4_i = 0; sq4_i < 4; sq4_i++)\
  {\
    double d4_adiv, d4_adivC;\
    DOUBLE_4MULTFMA(a,aC,div,divC,d4_adiv,d4_adivC);\
    d4_adiv=-d4_adiv; d4_adivC=-d4_adivC;\
    double d4_e, d4_eC;\
    DOUBLE_FAST4SUM(d4_one,d4_zero,d4_adiv,d4_adivC,d4_e,d4_eC);\
    double d4_dive, d4_diveC;\
    DOUBLE_4MULTFMA(div,divC,d4_e,d4_eC,d4_dive,d4_diveC);\
    double d4_div=div, d4_divC=divC;\
    DOUBLE_FAST4SUM(d4_div,d4_divC,d4_dive,d4_diveC,div,divC);\
  }\
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


// Newton-Raphson iteration for precise calculation of sqrt(a+aC) and 1/sqrt(a+aC)
#define MM256_4SQRT_NEWTONRAPHSON_FMA(a,aC,sqrt_a,sqrt_aC,divsqrt_a,divsqrt_aC)\
{\
  /* generate initial guess */ \
  __m256d sq4_g = _mm256_sqrt_pd(a);\
  __m256d sq4_gC = _mm256_setzero_pd();\
  __m256d sq4_half, sq4_halfC;\
  __m256d sq4_zero = _mm256_setzero_pd();\
  __m256d sq4_one = _mm256_set1_pd(1.);\
  __m256d sq4_two = _mm256_set1_pd(2.);\
  MM256_2DIVFMA(sq4_one,sq4_two,sq4_half,sq4_halfC);\
  __m256d sq4_h = _mm256_div_pd(sq4_half,sq4_g);\
  __m256d sq4_hC = _mm256_setzero_pd();\
  /* use some iterations of Newton-Raphson (quadratic convergence) */ \
  for(int sq4_i = 0; sq4_i < 4; sq4_i++)\
  {\
    __m256d sq4_gh, sq4_ghC;\
    MM256_4MULTFMA(sq4_g,sq4_gC,sq4_h,sq4_hC,sq4_gh,sq4_ghC);\
    sq4_gh=-sq4_gh; sq4_ghC=-sq4_ghC;\
    __m256d sq4_r, sq4_rC;\
    MM256_FAST2SUM(sq4_half,sq4_halfC,sq4_gh,sq4_ghC,sq4_r,sq4_rC);\
    __m256d sq4_gr, sq4_grC;\
    MM256_4MULTFMA(sq4_g,sq4_gC,sq4_r,sq4_rC,sq4_gr,sq4_grC);\
    __m256d sq4_hr, sq4_hrC;\
    MM256_4MULTFMA(sq4_h,sq4_hC,sq4_r,sq4_rC,sq4_hr,sq4_hrC);\
    __m256d sq4_g_ = sq4_g, sq4_gC_ = sq4_gC;\
    MM256_FAST4SUM(sq4_g_,sq4_gC_,sq4_gr,sq4_grC,sq4_g,sq4_gC);\
    __m256d sq4_h_ = sq4_h, sq4_hC_ = sq4_hC;\
    MM256_FAST4SUM(sq4_h_,sq4_hC_,sq4_hr,sq4_hrC,sq4_h,sq4_hC);\
  }\
  sqrt_a = sq4_g; sqrt_aC = sq4_gC;\
  MM256_4MULTFMA(sq4_two,sq4_zero,sq4_h,sq4_hC,divsqrt_a,divsqrt_aC);\
}
#define MM128_4SQRT_NEWTONRAPHSON_FMA(a,aC,sqrt_a,sqrt_aC,divsqrt_a,divsqrt_aC)\
{\
  /* generate initial guess */ \
  __m128d sq4_g = _mm_sqrt_pd(a);\
  __m128d sq4_gC = _mm_setzero_pd();\
  __m128d sq4_half, sq4_halfC;\
  __m128d sq4_zero = _mm_setzero_pd();\
  __m128d sq4_one = _mm_set1_pd(1.);\
  __m128d sq4_two = _mm_set1_pd(2.);\
  MM128_2DIVFMA(sq4_one,sq4_two,sq4_half,sq4_halfC);\
  __m128d sq4_h = _mm_div_pd(sq4_half,sq4_g);\
  __m128d sq4_hC = _mm_setzero_pd();\
  /* use some iterations of Newton-Raphson (quadratic convergence) */ \
  for(int sq4_i = 0; sq4_i < 4; sq4_i++)\
  {\
    __m128d sq4_gh, sq4_ghC;\
    MM128_4MULTFMA(sq4_g,sq4_gC,sq4_h,sq4_hC,sq4_gh,sq4_ghC);\
    sq4_gh=-sq4_gh; sq4_ghC=-sq4_ghC;\
    __m128d sq4_r, sq4_rC;\
    MM128_FAST2SUM(sq4_half,sq4_halfC,sq4_gh,sq4_ghC,sq4_r,sq4_rC);\
    __m128d sq4_gr, sq4_grC;\
    MM128_4MULTFMA(sq4_g,sq4_gC,sq4_r,sq4_rC,sq4_gr,sq4_grC);\
    __m128d sq4_hr, sq4_hrC;\
    MM128_4MULTFMA(sq4_h,sq4_hC,sq4_r,sq4_rC,sq4_hr,sq4_hrC);\
    __m128d sq4_g_ = sq4_g, sq4_gC_ = sq4_gC;\
    MM128_FAST4SUM(sq4_g_,sq4_gC_,sq4_gr,sq4_grC,sq4_g,sq4_gC);\
    __m128d sq4_h_ = sq4_h, sq4_hC_ = sq4_hC;\
    MM128_FAST4SUM(sq4_h_,sq4_hC_,sq4_hr,sq4_hrC,sq4_h,sq4_hC);\
  }\
  sqrt_a = sq4_g; sqrt_aC = sq4_gC;\
  MM128_4MULTFMA(sq4_two,sq4_zero,sq4_h,sq4_hC,divsqrt_a,divsqrt_aC);\
}
#define DOUBLE_4SQRT_NEWTONRAPHSON_FMA(a,aC,sqrt_a,sqrt_aC,divsqrt_a,divsqrt_aC)\
{\
  /* generate initial guess */ \
  double sq4_g = sqrt(a);\
  double sq4_gC = 0.;\
  double sq4_half, sq4_halfC;\
  double sq4_zero = 0.;\
  double sq4_one = 1.;\
  double sq4_two = 2.;\
  DOUBLE_2DIVFMA(sq4_one,sq4_two,sq4_half,sq4_halfC);\
  double sq4_h = sq4_half/sq4_g;\
  double sq4_hC = 0.;\
  /* use some iterations of Newton-Raphson (quadratic convergence) */ \
  for(int sq4_i = 0; sq4_i < 4; sq4_i++)\
  {\
    double sq4_gh, sq4_ghC;\
    DOUBLE_4MULTFMA(sq4_g,sq4_gC,sq4_h,sq4_hC,sq4_gh,sq4_ghC);\
    sq4_gh=-sq4_gh; sq4_ghC=-sq4_ghC;\
    double sq4_r, sq4_rC;\
    DOUBLE_FAST4SUM(sq4_half,sq4_halfC,sq4_gh,sq4_ghC,sq4_r,sq4_rC);\
    double sq4_gr, sq4_grC;\
    DOUBLE_4MULTFMA(sq4_g,sq4_gC,sq4_r,sq4_rC,sq4_gr,sq4_grC);\
    double sq4_hr, sq4_hrC;\
    DOUBLE_4MULTFMA(sq4_h,sq4_hC,sq4_r,sq4_rC,sq4_hr,sq4_hrC);\
    double sq4_g_ = sq4_g, sq4_gC_ = sq4_gC;\
    DOUBLE_FAST4SUM(sq4_g_,sq4_gC_,sq4_gr,sq4_grC,sq4_g,sq4_gC);\
    double sq4_h_ = sq4_h, sq4_hC_ = sq4_hC;\
    DOUBLE_FAST4SUM(sq4_h_,sq4_hC_,sq4_hr,sq4_hrC,sq4_h,sq4_hC);\
  }\
  sqrt_a = sq4_g; sqrt_aC = sq4_gC;\
  DOUBLE_4MULTFMA(sq4_two,sq4_zero,sq4_h,sq4_hC,divsqrt_a,divsqrt_aC);\
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
