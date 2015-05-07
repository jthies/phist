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
    r_[i] = _mm256_set_pd(r[i+12],r[i+8],r[i+4],r[i+0]);
    rC_[i] = _mm256_set_pd(rC[i+12],rC[i+8],rC[i+4],rC[i+0]);
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
    r_[i] = _mm_set_pd(r[i+2],r[i]);
    rC_[i] = _mm_set_pd(rC[i+2],rC[i]);
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


  int nrows4 = nrows/4;
#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows4; i++)
  {
    __m256d xi = _mm256_load_pd(&x[4*i]);
    __m256d xir, xirC;
    MM256_2MULTFMA(xi,r_,xir,xirC);
    __m256d tmp = _mm256_fmadd_pd(xi,rC_,xirC);
    __m256d newX = _mm256_add_pd(xir,tmp);
    _mm256_store_pd(&x[4*i],newX);
  }
}

// more accurate gemm product y <- alpha*x*m + beta*y AVX2 kernel for y of blocksize 4
void dgemm_sb_prec_k_4(int nrows, int k, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y)
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
  if( !is_aligned(y,32) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
    exit(1);
    return;
  }


  // buffer rows of r and multiply by alpha
  __m256d r_[k], rC_[k];
  __m256d alpha_ = _mm256_set1_pd(alpha);
  __m256d alphaC_ = _mm256_setzero_pd();
  for(int i = 0; i < k; i++)
  {
    __m256d tmp = _mm256_set_pd(r[i+3*k],r[i+2*k],r[i+k],r[i]);
    __m256d tmpC = _mm256_set_pd(rC[i+3*k],rC[i+2*k],rC[i+k],rC[i]);
    MM256_4MULTFMA(alpha_,alphaC_,tmp,tmpC,r_[i],rC_[i]);
  }
  __m256d beta_ = _mm256_set1_pd(beta);

#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i++)
  {
    __m256d oldY = _mm256_load_pd(&y[4*i]);
    __m256d s, t;
    MM256_2MULTFMA(beta_,oldY,s,t);
    for(int j = 0; j < k; j++)
    {
      __m256d xi = _mm256_broadcast_sd(&x[k*i+j]);
      __m256d xij, xijC;
      MM256_2MULTFMA(xi,r_[j],xij,xijC);
      __m256d xijC_ = _mm256_fmadd_pd(xi,rC_[j],xijC);
      __m256d oldS = s, t_;
      MM256_FAST2SUM(oldS,xij,s,t_);
      __m256d tmp = _mm256_add_pd(t,t_);
      t = _mm256_add_pd(tmp,xijC_);
    }
    __m256d newY = _mm256_add_pd(s,t);
    _mm256_store_pd(&y[4*i],newY);
  }
}


// more accurate gemm product y <- alpha*x*m + beta*y AVX2 kernel for y of blocksize 2
void dgemm_sb_prec_k_2(int nrows, int k, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y)
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
  if( !is_aligned(y,32) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
    exit(1);
    return;
  }

  // buffer rows of r and multiply by alpha
  __m128d r_[k], rC_[k];
  __m128d alpha_ = _mm_set1_pd(alpha);
  __m128d alphaC_ = _mm_setzero_pd();
  for(int i = 0; i < k; i++)
  {
    __m128d tmp = _mm_set_pd(r[i+k],r[i]);
    __m128d tmpC = _mm_set_pd(rC[i+k],rC[i]);
    MM128_4MULTFMA(alpha_,alphaC_,tmp,tmpC,r_[i],rC_[i]);
  }
  __m128d beta_ = _mm_set1_pd(beta);

#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i++)
  {
    __m128d oldY = _mm_load_pd(&y[2*i]);
    __m128d s, t;
    MM128_2MULTFMA(beta_,oldY,s,t);
    for(int j = 0; j < k; j++)
    {
      __m128d xi = _mm_load1_pd(&x[k*i+j]);
      __m128d xij, xijC;
      MM128_2MULTFMA(xi,r_[j],xij,xijC);
      __m128d xijC_ = _mm_fmadd_pd(xi,rC_[j],xijC);
      __m128d oldS = s, t_;
      MM128_FAST2SUM(oldS,xij,s,t_);
      __m128d tmp = _mm_add_pd(t,t_);
      t = _mm_add_pd(tmp,xijC_);
    }
    __m128d newY = _mm_add_pd(s,t);
    _mm_store_pd(&y[2*i],newY);
  }
}


// more accurate gemm product y <- alpha*x*m + beta*y AVX2 kernel for y of blocksize 2
void dgemm_sb_prec_k_1(int nrows, int k, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y)
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
  if( !is_aligned(y,32) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
    exit(1);
    return;
  }
  // buffer rows of r and multiply by alpha
  __m256d r_[k], rC_[k];
  __m256d alpha_ = _mm256_set1_pd(alpha);
  __m256d alphaC_ = _mm256_setzero_pd();
  for(int i = 0; i < k; i++)
  {
    __m256d tmp = _mm256_set1_pd(r[i]);
    __m256d tmpC = _mm256_set1_pd(rC[i]);
    MM256_4MULTFMA(alpha_,alphaC_,tmp,tmpC,r_[i],rC_[i]);
  }
  __m256d beta_ = _mm256_set1_pd(beta);


  int nrows4 = nrows/4;
#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows4; i++)
  {
    __m256d oldY = _mm256_load_pd(&y[4*i]);
    __m256d s, t;
    MM256_2MULTFMA(beta_,oldY,s,t);
    for(int j = 0; j < k; j++)
    {
      __m256d xi = _mm256_set_pd(x[k*4*i+3*k+j],x[k*4*i+2*k+j],x[k*4*i+k+j],x[k*4*i+j]);
      __m256d xir, xirC;
      MM256_2MULTFMA(xi,r_[j],xir,xirC);
      __m256d xirC_ = _mm256_fmadd_pd(xi,rC_[j],xirC);
      __m256d oldS = s, t_;
      MM256_FAST2SUM(oldS,xir,s,t_);
      __m256d tmp = _mm256_add_pd(t,t_);
      t = _mm256_add_pd(tmp,xirC_);
    }
    __m256d newY = _mm256_add_pd(s,t);
    _mm256_store_pd(&y[4*i],newY);
  }
}

// more accurate gemm product y <- alpha*x*m AVX2 kernel for y of blocksize 4 with non-temporal stores
void dgemm_sb_prec_k_4_nt(int nrows, int k, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y)
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
  if( !is_aligned(y,32) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
    exit(1);
    return;
  }


  // buffer rows of r and multiply by alpha
  __m256d r_[k], rC_[k];
  __m256d alpha_ = _mm256_set1_pd(alpha);
  __m256d alphaC_ = _mm256_setzero_pd();
  for(int i = 0; i < k; i++)
  {
    __m256d tmp = _mm256_set_pd(r[i+3*k],r[i+2*k],r[i+k],r[i]);
    __m256d tmpC = _mm256_set_pd(rC[i+3*k],rC[i+2*k],rC[i+k],rC[i]);
    MM256_4MULTFMA(alpha_,alphaC_,tmp,tmpC,r_[i],rC_[i]);
  }

#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i++)
  {
    __m256d s = _mm256_setzero_pd();
    __m256d t = _mm256_setzero_pd();
    for(int j = 0; j < k; j++)
    {
      __m256d xi = _mm256_broadcast_sd(&x[k*i+j]);
      __m256d xij, xijC;
      MM256_2MULTFMA(xi,r_[j],xij,xijC);
      __m256d xijC_ = _mm256_fmadd_pd(xi,rC_[j],xijC);
      __m256d oldS = s, t_;
      MM256_FAST2SUM(oldS,xij,s,t_);
      __m256d tmp = _mm256_add_pd(t,t_);
      t = _mm256_add_pd(tmp,xijC_);
    }
    __m256d newY = _mm256_add_pd(s,t);
    _mm256_stream_pd(&y[4*i],newY);
  }
}


// more accurate gemm product y <- alpha*x*m AVX2 kernel for y of blocksize 2 with non-temporal stores
void dgemm_sb_prec_k_2_nt(int nrows, int k, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y)
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
  if( !is_aligned(y,32) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
    exit(1);
    return;
  }

  // buffer rows of r and multiply by alpha
  __m128d r_[k], rC_[k];
  __m128d alpha_ = _mm_set1_pd(alpha);
  __m128d alphaC_ = _mm_setzero_pd();
  for(int i = 0; i < k; i++)
  {
    __m128d tmp = _mm_set_pd(r[i+k],r[i]);
    __m128d tmpC = _mm_set_pd(rC[i+k],rC[i]);
    MM128_4MULTFMA(alpha_,alphaC_,tmp,tmpC,r_[i],rC_[i]);
  }

#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i++)
  {
    __m128d s = _mm_setzero_pd();
    __m128d t = _mm_setzero_pd();
    for(int j = 0; j < k; j++)
    {
      __m128d xi = _mm_load1_pd(&x[k*i+j]);
      __m128d xij, xijC;
      MM128_2MULTFMA(xi,r_[j],xij,xijC);
      __m128d xijC_ = _mm_fmadd_pd(xi,rC_[j],xijC);
      __m128d oldS = s, t_;
      MM128_FAST2SUM(oldS,xij,s,t_);
      __m128d tmp = _mm_add_pd(t,t_);
      t = _mm_add_pd(tmp,xijC_);
    }
    __m128d newY = _mm_add_pd(s,t);
    _mm_stream_pd(&y[2*i],newY);
  }
}


// more accurate gemm product y <- alpha*x*m AVX2 kernel for y of blocksize 2 with non-temporal stores
void dgemm_sb_prec_k_1_nt(int nrows, int k, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y)
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
  if( !is_aligned(y,32) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
    exit(1);
    return;
  }

  // buffer rows of r and multiply by alpha
  __m256d r_[k], rC_[k];
  __m256d alpha_ = _mm256_set1_pd(alpha);
  __m256d alphaC_ = _mm256_setzero_pd();
  for(int i = 0; i < k; i++)
  {
    __m256d tmp = _mm256_set1_pd(r[i]);
    __m256d tmpC = _mm256_set1_pd(rC[i]);
    MM256_4MULTFMA(alpha_,alphaC_,tmp,tmpC,r_[i],rC_[i]);
  }


  int nrows4 = nrows/4;
#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows4; i++)
  {
    __m256d s = _mm256_setzero_pd();
    __m256d t = _mm256_setzero_pd();
    for(int j = 0; j < k; j++)
    {
      __m256d xi = _mm256_set_pd(x[k*4*i+3*k+j],x[k*4*i+2*k+j],x[k*4*i+k+j],x[k*4*i+j]);
      __m256d xir, xirC;
      MM256_2MULTFMA(xi,r_[j],xir,xirC);
      __m256d xirC_ = _mm256_fmadd_pd(xi,rC_[j],xirC);
      __m256d oldS = s, t_;
      MM256_FAST2SUM(oldS,xir,s,t_);
      __m256d tmp = _mm256_add_pd(t,t_);
      t = _mm256_add_pd(tmp,xirC_);
    }
    __m256d newY = _mm256_add_pd(s,t);
    _mm256_stream_pd(&y[4*i],newY);
  }
}

// more accurate gemm product y <- alpha*x*m + beta*y AVX2 kernel for y of blocksize 4
void dgemm_sb_prec_k_strided_4(int nrows, int k, double alpha, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double beta, double *restrict y)
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
  if( !is_aligned(y,32) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
    exit(1);
    return;
  }


  // buffer rows of r and multiply by alpha
  __m256d r_[k], rC_[k];
  __m256d alpha_ = _mm256_set1_pd(alpha);
  __m256d alphaC_ = _mm256_setzero_pd();
  for(int i = 0; i < k; i++)
  {
    __m256d tmp = _mm256_set_pd(r[i+3*k],r[i+2*k],r[i+k],r[i]);
    __m256d tmpC = _mm256_set_pd(rC[i+3*k],rC[i+2*k],rC[i+2],rC[i]);
    MM256_4MULTFMA(alpha_,alphaC_,tmp,tmpC,r_[i],rC_[i]);
  }
  __m256d beta_ = _mm256_set1_pd(beta);

#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i++)
  {
    __m256d oldY = _mm256_load_pd(&y[4*i]);
    __m256d s, t;
    MM256_2MULTFMA(beta_,oldY,s,t);
    for(int j = 0; j < k; j++)
    {
      __m256d xi = _mm256_broadcast_sd(&x[ldx*i+j]);
      __m256d xij, xijC;
      MM256_2MULTFMA(xi,r_[j],xij,xijC);
      __m256d xijC_ = _mm256_fmadd_pd(xi,rC_[j],xijC);
      __m256d oldS = s, t_;
      MM256_FAST2SUM(oldS,xij,s,t_);
      __m256d tmp = _mm256_add_pd(t,t_);
      t = _mm256_add_pd(tmp,xijC_);
    }
    __m256d newY = _mm256_add_pd(s,t);
    _mm256_store_pd(&y[4*i],newY);
  }
}


// more accurate gemm product y <- alpha*x*m + beta*y AVX2 kernel for y of blocksize 2
void dgemm_sb_prec_k_strided_2(int nrows, int k, double alpha, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double beta, double *restrict y)
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
  if( !is_aligned(y,32) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
    exit(1);
    return;
  }

  // buffer rows of r and multiply by alpha
  __m128d r_[k], rC_[k];
  __m128d alpha_ = _mm_set1_pd(alpha);
  __m128d alphaC_ = _mm_setzero_pd();
  for(int i = 0; i < k; i++)
  {
    __m128d tmp = _mm_set_pd(r[i+k],r[i]);
    __m128d tmpC = _mm_set_pd(rC[i+k],rC[i]);
    MM128_4MULTFMA(alpha_,alphaC_,tmp,tmpC,r_[i],rC_[i]);
  }
  __m128d beta_ = _mm_set1_pd(beta);

#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i++)
  {
    __m128d oldY = _mm_load_pd(&y[2*i]);
    __m128d s, t;
    MM128_2MULTFMA(beta_,oldY,s,t);
    for(int j = 0; j < k; j++)
    {
      __m128d xi = _mm_load1_pd(&x[ldx*i+j]);
      __m128d xij, xijC;
      MM128_2MULTFMA(xi,r_[j],xij,xijC);
      __m128d xijC_ = _mm_fmadd_pd(xi,rC_[j],xijC);
      __m128d oldS = s, t_;
      MM128_FAST2SUM(oldS,xij,s,t_);
      __m128d tmp = _mm_add_pd(t,t_);
      t = _mm_add_pd(tmp,xijC_);
    }
    __m128d newY = _mm_add_pd(s,t);
    _mm_store_pd(&y[2*i],newY);
  }
}


// more accurate gemm product y <- alpha*x*m + beta*y AVX2 kernel for y of blocksize 2
void dgemm_sb_prec_k_strided_1(int nrows, int k, double alpha, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double beta, double *restrict y)
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
  if( !is_aligned(y,32) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
    exit(1);
    return;
  }
  // buffer rows of r and multiply by alpha
  __m256d r_[k], rC_[k];
  __m256d alpha_ = _mm256_set1_pd(alpha);
  __m256d alphaC_ = _mm256_setzero_pd();
  for(int i = 0; i < k; i++)
  {
    __m256d tmp = _mm256_set1_pd(r[i]);
    __m256d tmpC = _mm256_set1_pd(rC[i]);
    MM256_4MULTFMA(alpha_,alphaC_,tmp,tmpC,r_[i],rC_[i]);
  }
  __m256d beta_ = _mm256_set1_pd(beta);


  int nrows4 = nrows/4;
#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows4; i++)
  {
    __m256d oldY = _mm256_load_pd(&y[4*i]);
    __m256d s, t;
    MM256_2MULTFMA(beta_,oldY,s,t);
    for(int j = 0; j < k; j++)
    {
      __m256d xi = _mm256_set_pd(x[ldx*4*i+3*ldx+j],x[ldx*4*i+2*ldx+j],x[ldx*4*i+ldx+j],x[ldx*4*i+j]);
      __m256d xir, xirC;
      MM256_2MULTFMA(xi,r_[j],xir,xirC);
      __m256d xirC_ = _mm256_fmadd_pd(xi,rC_[j],xirC);
      __m256d oldS = s, t_;
      MM256_FAST2SUM(oldS,xir,s,t_);
      __m256d tmp = _mm256_add_pd(t,t_);
      t = _mm256_add_pd(tmp,xirC_);
    }
    __m256d newY = _mm256_add_pd(s,t);
    _mm256_store_pd(&y[4*i],newY);
  }
}

// more accurate gemm product y <- alpha*x*m AVX2 kernel for y of blocksize 4 with non-temporal stores
void dgemm_sb_prec_k_strided_4_nt(int nrows, int k, double alpha, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double *restrict y)
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
  if( !is_aligned(y,32) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
    exit(1);
    return;
  }


  // buffer rows of r and multiply by alpha
  __m256d r_[k], rC_[k];
  __m256d alpha_ = _mm256_set1_pd(alpha);
  __m256d alphaC_ = _mm256_setzero_pd();
  for(int i = 0; i < k; i++)
  {
    __m256d tmp = _mm256_set_pd(r[i+3*k],r[i+2*k],r[i+k],r[i]);
    __m256d tmpC = _mm256_set_pd(rC[i+3*k],rC[i+2*k],rC[i+k],rC[i]);
    MM256_4MULTFMA(alpha_,alphaC_,tmp,tmpC,r_[i],rC_[i]);
  }

#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i++)
  {
    __m256d s = _mm256_setzero_pd();
    __m256d t = _mm256_setzero_pd();
    for(int j = 0; j < k; j++)
    {
      __m256d xi = _mm256_broadcast_sd(&x[ldx*i+j]);
      __m256d xij, xijC;
      MM256_2MULTFMA(xi,r_[j],xij,xijC);
      __m256d xijC_ = _mm256_fmadd_pd(xi,rC_[j],xijC);
      __m256d oldS = s, t_;
      MM256_FAST2SUM(oldS,xij,s,t_);
      __m256d tmp = _mm256_add_pd(t,t_);
      t = _mm256_add_pd(tmp,xijC_);
    }
    __m256d newY = _mm256_add_pd(s,t);
    _mm256_stream_pd(&y[4*i],newY);
  }
}


// more accurate gemm product y <- alpha*x*m AVX2 kernel for y of blocksize 2 with non-temporal stores
void dgemm_sb_prec_k_strided_2_nt(int nrows, int k, double alpha, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double *restrict y)
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
  if( !is_aligned(y,32) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
    exit(1);
    return;
  }

  // buffer rows of r and multiply by alpha
  __m128d r_[k], rC_[k];
  __m128d alpha_ = _mm_set1_pd(alpha);
  __m128d alphaC_ = _mm_setzero_pd();
  for(int i = 0; i < k; i++)
  {
    __m128d tmp = _mm_set_pd(r[i+k],r[i]);
    __m128d tmpC = _mm_set_pd(rC[i+k],rC[i]);
    MM128_4MULTFMA(alpha_,alphaC_,tmp,tmpC,r_[i],rC_[i]);
  }

#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i++)
  {
    __m128d s = _mm_setzero_pd();
    __m128d t = _mm_setzero_pd();
    for(int j = 0; j < k; j++)
    {
      __m128d xi = _mm_load1_pd(&x[ldx*i+j]);
      __m128d xij, xijC;
      MM128_2MULTFMA(xi,r_[j],xij,xijC);
      __m128d xijC_ = _mm_fmadd_pd(xi,rC_[j],xijC);
      __m128d oldS = s, t_;
      MM128_FAST2SUM(oldS,xij,s,t_);
      __m128d tmp = _mm_add_pd(t,t_);
      t = _mm_add_pd(tmp,xijC_);
    }
    __m128d newY = _mm_add_pd(s,t);
    _mm_stream_pd(&y[2*i],newY);
  }
}


// more accurate gemm product y <- alpha*x*m AVX2 kernel for y of blocksize 2 with non-temporal stores
void dgemm_sb_prec_k_strided_1_nt(int nrows, int k, double alpha, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double *restrict y)
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
  if( !is_aligned(y,32) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
    exit(1);
    return;
  }
  // buffer rows of r and multiply by alpha
  __m256d r_[k], rC_[k];
  __m256d alpha_ = _mm256_set1_pd(alpha);
  __m256d alphaC_ = _mm256_setzero_pd();
  for(int i = 0; i < k; i++)
  {
    __m256d tmp = _mm256_set1_pd(r[i]);
    __m256d tmpC = _mm256_set1_pd(rC[i]);
    MM256_4MULTFMA(alpha_,alphaC_,tmp,tmpC,r_[i],rC_[i]);
  }


  int nrows4 = nrows/4;
#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows4; i++)
  {
    __m256d s = _mm256_setzero_pd();
    __m256d t = _mm256_setzero_pd();
    for(int j = 0; j < k; j++)
    {
      __m256d xi = _mm256_set_pd(x[ldx*4*i+3*ldx+j],x[ldx*4*i+2*ldx+j],x[ldx*4*i+ldx+j],x[ldx*4*i+j]);
      __m256d xir, xirC;
      MM256_2MULTFMA(xi,r_[j],xir,xirC);
      __m256d xirC_ = _mm256_fmadd_pd(xi,rC_[j],xirC);
      __m256d oldS = s, t_;
      MM256_FAST2SUM(oldS,xir,s,t_);
      __m256d tmp = _mm256_add_pd(t,t_);
      t = _mm256_add_pd(tmp,xirC_);
    }
    __m256d newY = _mm256_add_pd(s,t);
    _mm256_stream_pd(&y[4*i],newY);
  }
}




