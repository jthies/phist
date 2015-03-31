/*! \file gemm_kernels_sC_prec.c
 * Fast parallel BLAS-gemm like functions with high precision for different blocksizes for mvec_module
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


// more accurate gemm product x <- x*m AVX2 kernel for blocksize 4
void dgemm_sb_inplace_prec_4(int nrows, double *restrict x, const double *restrict r, const double *restrict rC)
{
#if defined(TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  if( !is_aligned(x,32) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)x);
    exit(1);
    return;
  }

  // buffer rows of r
  __m256d r_[4], rC_[4];
  for(int i = 0; i < 4; i++)
  {
    r_[i] = _mm256_set_pd(r[i],r[i+4],r[i+8],r[i+12]);
    rC_[i] = _mm256_set_pd(rC[i],rC[i+4],rC[i+8],rC[i+12]);
  }

#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i++)
  {
    __m256d s = _mm256_setzero_pd();
    __m256d t = _mm256_setzero_pd();
    for(int j = 0; j < 4; j++)
    {
      __m256d xi = _mm256_broadcast_sd(&x[4*i+j]);
      __m256d xij, xijC;
      MM256_2MULTFMA(xi,r_[j],xij,xijC);
      __m256d xijC_ = _mm256_fmadd_pd(xi,rC_[j],xijC);
      __m256d oldS = s, t_;
      MM256_FAST2SUM(oldS,xij,s,t_);
      __m256d tmp = _mm256_add_pd(t,t_);
      t = _mm256_add_pd(tmp,xijC_);
    }
    __m256d newX = _mm256_add_pd(s,t);
    _mm256_store_pd(&x[4*i],newX);
  }
}


// more accurate gemm product x <- x*m AVX2 kernel for blocksize 2
void dgemm_sb_inplace_prec_2(int nrows, double *restrict x, const double *restrict r, const double *restrict rC)
{
#if defined(TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  if( !is_aligned(x,32) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)x);
    exit(1);
    return;
  }

  // buffer rows of r
  __m128d r_[2], rC_[2];
  for(int i = 0; i < 2; i++)
  {
    r_[i] = _mm_set_pd(r[i],r[i+2]);
    rC_[i] = _mm_set_pd(rC[i],rC[i+2]);
  }

#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i++)
  {
    __m128d s = _mm_setzero_pd();
    __m128d t = _mm_setzero_pd();
    for(int j = 0; j < 2; j++)
    {
      __m128d xi = _mm_load1_pd(&x[2*i+j]);
      __m128d xij, xijC;
      MM128_2MULTFMA(xi,r_[j],xij,xijC);
      __m128d xijC_ = _mm_fmadd_pd(xi,rC_[j],xijC);
      __m128d oldS = s, t_;
      MM128_FAST2SUM(oldS,xij,s,t_);
      __m128d tmp = _mm_add_pd(t,t_);
      t = _mm_add_pd(tmp,xijC_);
    }
    __m128d newX = _mm_add_pd(s,t);
    _mm_store_pd(&x[2*i],newX);
  }
}


// more accurate gemm product x <- x*m AVX2 kernel for blocksize 1
void dgemm_sb_inplace_prec_1(int nrows, double *restrict x, const double *restrict r, const double *restrict rC)
{
#if defined(TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  if( !is_aligned(x,32) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)x);
    exit(1);
    return;
  }

  // buffer rows of r
  __m256d r_ = _mm256_broadcast_sd(r);
  __m256d rC_ = _mm256_broadcast_sd(rC);

#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows/4; i++)
  {
    __m256d xi = _mm256_load_pd(&x[4*i]);
    __m256d xir, xirC;
    MM256_2MULTFMA(xi,r_,xir,xirC);
    __m256d tmp = _mm256_fmadd_pd(xi,rC_,xirC);
    __m256d newX = _mm256_add_pd(xir,tmp);
    _mm256_store_pd(&x[4*i],newX);
  }
}


