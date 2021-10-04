/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
/*! \file gemm_kernels_sB_prec.c
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
#include "phist_defs.h"
#include "phist_prec_helpers.h"


// more accurate gemm product x <- x*m AVX2 kernel for blocksize 4
void dgemm_sb_inplace_prec_4(int nrows, double *restrict x, const double *restrict r, const double *restrict rC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
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
    __m256d x_ = _mm256_load_pd(&x[4*i]);

    __m256d s, t;
    __m256d p, pi;
    // unroll j = 0
    __m256d xi = _mm256_permute4x64_pd(x_,0);
    MM256_2MULTFMA(xi,r_[0],p,pi);
    __m256d pi_ = _mm256_fmadd_pd(xi,rC_[0],pi);
    s = p, t = pi_;

    // j = 1
    xi = _mm256_permute4x64_pd(x_,1+4*1+16*1+64*1);
    MM256_2MULTFMA(xi,r_[1],p,pi);
    pi_ = _mm256_fmadd_pd(xi,rC_[1],pi);
    __m256d oldS = s, sigma;
    MM256_2SUM(oldS,p,s,sigma);
    __m256d tmp = _mm256_add_pd(pi_,sigma);
    t = _mm256_add_pd(t, tmp);

    // j = 2
    xi = _mm256_permute4x64_pd(x_,2+4*2+16*2+64*2);
    MM256_2MULTFMA(xi,r_[2],p,pi);
    pi_ = _mm256_fmadd_pd(xi,rC_[2],pi);
    oldS = s;
    MM256_2SUM(oldS,p,s,sigma);
    tmp = _mm256_add_pd(pi_,sigma);
    t = _mm256_add_pd(t, tmp);

    // j = 3
    xi = _mm256_permute4x64_pd(x_,3+4*3+16*3+64*3);
    MM256_2MULTFMA(xi,r_[3],p,pi);
    pi_ = _mm256_fmadd_pd(xi,rC_[3],pi);
    oldS = s;
    MM256_2SUM(oldS,p,s,sigma);
    tmp = _mm256_add_pd(pi_,sigma);
    t = _mm256_add_pd(t, tmp);

    __m256d newX = _mm256_add_pd(s,t);
    _mm256_store_pd(&x[4*i],newX);
  }
}


// more accurate gemm product x <- x*m AVX2 kernel for blocksize 2
void dgemm_sb_inplace_prec_2(int nrows, double *restrict x, const double *restrict r, const double *restrict rC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  if( !is_aligned(x,32) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)x);
    exit(1);
    return;
  }

  // buffer rows of r
  __m256d r_[2], rC_[2];
  for(int i = 0; i < 2; i++)
  {
    r_[i] = _mm256_set_pd(r[i+2],r[i+0],r[i+2],r[i+0]);
    rC_[i] = _mm256_set_pd(rC[i+2],rC[i+0],rC[i+2],rC[i+0]);
  }

#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i+=2)
  {
    __m256d x_ = _mm256_load_pd(&x[2*i]);

    __m256d s, t;
    __m256d p, pi;
    // unroll j = 0
    __m256d xi = _mm256_permute_pd(x_,0);
    MM256_2MULTFMA(xi,r_[0],p,pi);
    __m256d pi_ = _mm256_fmadd_pd(xi,rC_[0],pi);
    s = p, t = pi_;

    // j = 1
    xi = _mm256_permute_pd(x_,1+2*1+4*1+8*1);
    MM256_2MULTFMA(xi,r_[1],p,pi);
    pi_ = _mm256_fmadd_pd(xi,rC_[1],pi);
    __m256d oldS = s, sigma;
    MM256_2SUM(oldS,p,s,sigma);
    __m256d tmp = _mm256_add_pd(pi_,sigma);
    t = _mm256_add_pd(t, tmp);

    __m256d newX = _mm256_add_pd(s,t);
    _mm256_store_pd(&x[2*i],newX);
  }
}


// more accurate gemm product x <- x*m AVX2 kernel for blocksize 1
void dgemm_sb_inplace_prec_1(int nrows, double *restrict x, const double *restrict r, const double *restrict rC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
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
  for(int i = 0; i < nrows; i+=4)
  {
    __m256d xi = _mm256_load_pd(&x[i]);
    __m256d xir, xirC;
    MM256_2MULTFMA(xi,r_,xir,xirC);
    __m256d tmp = _mm256_fmadd_pd(xi,rC_,xirC);
    __m256d newX = _mm256_add_pd(xir,tmp);
    _mm256_store_pd(&x[i],newX);
  }
}


// more accurate gemm product y <- alpha*x*m + beta*y AVX2 kernel for y of blocksize 4
void dgemm_sb_prec_k_4(int nrows, int k, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  if( !is_aligned(y,32) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
    exit(1);
    return;
  }


  // buffer rows of r and multiply by alpha
  __m256d r_[k], rC_[k];
  __m256d alpha_ = _mm256_set1_pd(alpha);
  for(int i = 0; i < k; i++)
  {
    __m256d tmp = _mm256_set_pd(r[i+3*k],r[i+2*k],r[i+k],r[i]);
    __m256d tmpC = _mm256_set_pd(rC[i+3*k],rC[i+2*k],rC[i+k],rC[i]);
    __m256d pi;
    MM256_2MULTFMA(alpha_,tmp,r_[i],pi);
    rC_[i] = _mm256_fmadd_pd(alpha_,tmpC,pi);
  }

  __m256d beta_ = _mm256_set1_pd(beta);
#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i++)
  {
    __m256d s, t;
    if( beta == 0. )
    {
      s = _mm256_setzero_pd();
      t = _mm256_setzero_pd();
    }
    else
    {
      __m256d oldY = _mm256_load_pd(&y[4*i]);
      MM256_2MULTFMA(beta_,oldY,s,t);
    }
    for(int j = 0; j < k; j++)
    {
      __m256d xi = _mm256_broadcast_sd(&x[k*i+j]);
      __m256d xij, xijC;
      MM256_2MULTFMA(xi,r_[j],xij,xijC);
      __m256d xijC_ = _mm256_fmadd_pd(xi,rC_[j],xijC);
      __m256d oldS = s, t_;
      MM256_2SUM(oldS,xij,s,t_);
      __m256d tmp = _mm256_add_pd(xijC_,t_);
      t = _mm256_add_pd(t,tmp);
    }
    __m256d newY = _mm256_add_pd(s,t);
    if( beta == 0. )
      _mm256_stream_pd(&y[4*i],newY);
    else
      _mm256_store_pd(&y[4*i],newY);
  }
}


// more accurate gemm product y <- alpha*x*m + beta*y AVX2 kernel for y of blocksize 2
void dgemm_sb_prec_k_2(int nrows, int k, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  if( !is_aligned(y,32) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
    exit(1);
    return;
  }


  // buffer rows of r and multiply by alpha
  __m256d r_[k], rC_[k];
  __m256d alpha_ = _mm256_set1_pd(alpha);
  for(int i = 0; i < k; i++)
  {
    __m256d tmp = _mm256_set_pd(r[i+k],r[i],r[i+k],r[i]);
    __m256d tmpC = _mm256_set_pd(rC[i+k],rC[i],rC[i+k],rC[i]);
    __m256d pi;
    MM256_2MULTFMA(alpha_,tmp,r_[i],pi);
    rC_[i] = _mm256_fmadd_pd(alpha_,tmpC,pi);
  }
  __m256d beta_ = _mm256_set1_pd(beta);

#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i+=2)
  {
    __m256d s, t;
    if( beta == 0. )
    {
      s = _mm256_setzero_pd();
      t = _mm256_setzero_pd();
    }
    else
    {
      __m256d oldY = _mm256_load_pd(&y[2*i]);
      MM256_2MULTFMA(beta_,oldY,s,t);
    }
    for(int j = 0; j < k; j++)
    {
      __m128d xil = _mm_load1_pd(&x[k*i+j]);
      __m128d xih = _mm_load1_pd(&x[k*(i+1)+j]);
      __m256d xi  = _mm256_set_m128d(xih,xil);
      __m256d xij, xijC;
      MM256_2MULTFMA(xi,r_[j],xij,xijC);
      __m256d xijC_ = _mm256_fmadd_pd(xi,rC_[j],xijC);
      __m256d oldS = s, t_;
      MM256_2SUM(oldS,xij,s,t_);
      __m256d tmp = _mm256_add_pd(xijC_,t_);
      t = _mm256_add_pd(t,tmp);
    }
    __m256d newY = _mm256_add_pd(s,t);
    if( beta == 0. )
      _mm256_stream_pd(&y[2*i],newY);
    else
      _mm256_store_pd(&y[2*i],newY);
  }
}


// more accurate gemm product y <- alpha*x*m + beta*y AVX2 kernel for y of blocksize 1
void dgemm_sb_prec_k_1(int nrows, int k, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  if( !is_aligned(y,32) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
    exit(1);
    return;
  }
  // buffer rows of r and multiply by alpha
  __m256d r_[k], rC_[k];
  __m256d alpha_ = _mm256_set1_pd(alpha);
  for(int i = 0; i < k; i++)
  {
    __m256d tmp = _mm256_set1_pd(r[i]);
    __m256d tmpC = _mm256_set1_pd(rC[i]);
    __m256d pi;
    MM256_2MULTFMA(alpha_,tmp,r_[i],pi);
    rC_[i] = _mm256_fmadd_pd(alpha_,tmpC,pi);
  }
  __m256d beta_ = _mm256_set1_pd(beta);


#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i+=4)
  {
    __m256d s, t;
    if( beta == 0. )
    {
      s = _mm256_setzero_pd();
      t = _mm256_setzero_pd();
    }
    else
    {
      __m256d oldY = _mm256_load_pd(&y[i]);
      MM256_2MULTFMA(beta_,oldY,s,t);
    }
    for(int j = 0; j < k; j++)
    {
      __m256d xi = _mm256_set_pd(x[k*i+3*k+j],x[k*i+2*k+j],x[k*i+k+j],x[k*i+j]);
      __m256d xij, xijC;
      MM256_2MULTFMA(xi,r_[j],xij,xijC);
      __m256d xijC_ = _mm256_fmadd_pd(xi,rC_[j],xijC);
      __m256d oldS = s, t_;
      MM256_2SUM(oldS,xij,s,t_);
      __m256d tmp = _mm256_add_pd(xijC_,t_);
      t = _mm256_add_pd(t,tmp);
    }
    __m256d newY = _mm256_add_pd(s,t);
    if( beta == 0. )
      _mm256_stream_pd(&y[i],newY);
    else
      _mm256_store_pd(&y[i],newY);
  }
}


// more accurate gemm product y <- alpha*x*m AVX2 kernel for y of blocksize 4 with non-temporal stores
void dgemm_sb_prec_k_4_nt(int nrows, int k, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  if( !is_aligned(y,32) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
    exit(1);
    return;
  }


  // buffer rows of r and multiply by alpha
  __m256d r_[k], rC_[k];
  __m256d alpha_ = _mm256_set1_pd(alpha);
  for(int i = 0; i < k; i++)
  {
    __m256d tmp = _mm256_set_pd(r[i+3*k],r[i+2*k],r[i+k],r[i]);
    __m256d tmpC = _mm256_set_pd(rC[i+3*k],rC[i+2*k],rC[i+k],rC[i]);
    __m256d pi;
    MM256_2MULTFMA(alpha_,tmp,r_[i],pi);
    rC_[i] = _mm256_fmadd_pd(alpha_,tmpC,pi);
  }

#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i++)
  {
    __m256d s, t;
    // j = 0
    __m256d xi = _mm256_broadcast_sd(&x[k*i+0]);
    __m256d pi;
    MM256_2MULTFMA(xi,r_[0],s,pi);
    t = _mm256_fmadd_pd(xi,rC_[0],pi);

    for(int j = 1; j < k; j++)
    {
      __m256d xi = _mm256_broadcast_sd(&x[k*i+j]);
      __m256d xij, xijC;
      MM256_2MULTFMA(xi,r_[j],xij,xijC);
      __m256d xijC_ = _mm256_fmadd_pd(xi,rC_[j],xijC);
      __m256d oldS = s, t_;
      MM256_2SUM(oldS,xij,s,t_);
      __m256d tmp = _mm256_add_pd(xijC_,t_);
      t = _mm256_add_pd(t,tmp);
    }
    __m256d newY = _mm256_add_pd(s,t);
    _mm256_stream_pd(&y[4*i],newY);
  }
}


// more accurate gemm product y <- alpha*x*m AVX2 kernel for y of blocksize 2 with non-temporal stores
void dgemm_sb_prec_k_2_nt(int nrows, int k, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  if( !is_aligned(y,32) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
    exit(1);
    return;
  }


  // buffer rows of r and multiply by alpha
  __m256d r_[k], rC_[k];
  __m256d alpha_ = _mm256_set1_pd(alpha);
  for(int i = 0; i < k; i++)
  {
    __m256d tmp = _mm256_set_pd(r[i+k],r[i],r[i+k],r[i]);
    __m256d tmpC = _mm256_set_pd(rC[i+k],rC[i],rC[i+k],rC[i]);
    __m256d pi;
    MM256_2MULTFMA(alpha_,tmp,r_[i],pi);
    rC_[i] = _mm256_fmadd_pd(alpha_,tmpC,pi);
  }

#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i+=2)
  {
    __m256d s, t;
    // j = 0
    __m128d xil = _mm_load1_pd(&x[k*i+0]);
    __m128d xih = _mm_load1_pd(&x[k*(i+1)+0]);
    __m256d xi  = _mm256_set_m128d(xih,xil);
    __m256d pi;
    MM256_2MULTFMA(xi,r_[0],s,pi);
    t = _mm256_fmadd_pd(xi,rC_[0],pi);

    for(int j = 1; j < k; j++)
    {
      xil = _mm_load1_pd(&x[k*i+j]);
      xih = _mm_load1_pd(&x[k*(i+1)+j]);
      xi  = _mm256_set_m128d(xih,xil);
      __m256d xij, xijC;
      MM256_2MULTFMA(xi,r_[j],xij,xijC);
      __m256d xijC_ = _mm256_fmadd_pd(xi,rC_[j],xijC);
      __m256d oldS = s, t_;
      MM256_2SUM(oldS,xij,s,t_);
      __m256d tmp = _mm256_add_pd(xijC_,t_);
      t = _mm256_add_pd(t,tmp);
    }
    __m256d newY = _mm256_add_pd(s,t);
    _mm256_stream_pd(&y[2*i],newY);
  }
}


// more accurate gemm product y <- alpha*x*m AVX2 kernel for y of blocksize 1 with non-temporal stores
void dgemm_sb_prec_k_1_nt(int nrows, int k, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  if( !is_aligned(y,32) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
    exit(1);
    return;
  }
  // buffer rows of r and multiply by alpha
  __m256d r_[k], rC_[k];
  __m256d alpha_ = _mm256_set1_pd(alpha);
  for(int i = 0; i < k; i++)
  {
    __m256d tmp = _mm256_set1_pd(r[i]);
    __m256d tmpC = _mm256_set1_pd(rC[i]);
    __m256d pi;
    MM256_2MULTFMA(alpha_,tmp,r_[i],pi);
    rC_[i] = _mm256_fmadd_pd(alpha_,tmpC,pi);
  }


#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i+=4)
  {
    __m256d s, t;
    // j = 0
    __m256d xi = _mm256_set_pd(x[k*i+3*k+0],x[k*i+2*k+0],x[k*i+k+0],x[k*i+0]);
    __m256d pi;
    MM256_2MULTFMA(xi,r_[0],s,pi);
    t = _mm256_fmadd_pd(xi,rC_[0],pi);

    for(int j = 1; j < k; j++)
    {
      xi = _mm256_set_pd(x[k*i+3*k+j],x[k*i+2*k+j],x[k*i+k+j],x[k*i+j]);
      __m256d xij, xijC;
      MM256_2MULTFMA(xi,r_[j],xij,xijC);
      __m256d xijC_ = _mm256_fmadd_pd(xi,rC_[j],xijC);
      __m256d oldS = s, t_;
      MM256_2SUM(oldS,xij,s,t_);
      __m256d tmp = _mm256_add_pd(xijC_,t_);
      t = _mm256_add_pd(t,tmp);
    }
    __m256d newY = _mm256_add_pd(s,t);
    _mm256_stream_pd(&y[i],newY);
  }
}



// more accurate gemm product y <- alpha*x*m + beta*y AVX2 kernel for y of blocksize 4
void dgemm_sb_prec_k_strided_4(int nrows, int k, double alpha, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double beta, double *restrict y)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  if( !is_aligned(y,32) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
    exit(1);
    return;
  }


  // buffer rows of r and multiply by alpha
  __m256d r_[k], rC_[k];
  __m256d alpha_ = _mm256_set1_pd(alpha);
  for(int i = 0; i < k; i++)
  {
    __m256d tmp = _mm256_set_pd(r[i+3*k],r[i+2*k],r[i+k],r[i]);
    __m256d tmpC = _mm256_set_pd(rC[i+3*k],rC[i+2*k],rC[i+k],rC[i]);
    __m256d pi;
    MM256_2MULTFMA(alpha_,tmp,r_[i],pi);
    rC_[i] = _mm256_fmadd_pd(alpha_,tmpC,pi);
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
      MM256_2SUM(oldS,xij,s,t_);
      __m256d tmp = _mm256_add_pd(xijC_,t_);
      t = _mm256_add_pd(t,tmp);
    }
    __m256d newY = _mm256_add_pd(s,t);
    _mm256_store_pd(&y[4*i],newY);
  }
}


// more accurate gemm product y <- alpha*x*m + beta*y AVX2 kernel for y of blocksize 2
void dgemm_sb_prec_k_strided_2(int nrows, int k, double alpha, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double beta, double *restrict y)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  if( !is_aligned(y,32) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
    exit(1);
    return;
  }


  // buffer rows of r and multiply by alpha
  __m256d r_[k], rC_[k];
  __m256d alpha_ = _mm256_set1_pd(alpha);
  for(int i = 0; i < k; i++)
  {
    __m256d tmp = _mm256_set_pd(r[i+k],r[i],r[i+k],r[i]);
    __m256d tmpC = _mm256_set_pd(rC[i+k],rC[i],rC[i+k],rC[i]);
    __m256d pi;
    MM256_2MULTFMA(alpha_,tmp,r_[i],pi);
    rC_[i] = _mm256_fmadd_pd(alpha_,tmpC,pi);
  }
  __m256d beta_ = _mm256_set1_pd(beta);

#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i+=2)
  {
    __m256d oldY = _mm256_load_pd(&y[2*i]);
    __m256d s, t;
    MM256_2MULTFMA(beta_,oldY,s,t);
    for(int j = 0; j < k; j++)
    {
      __m128d xil = _mm_load1_pd(&x[ldx*i+j]);
      __m128d xih = _mm_load1_pd(&x[ldx*(i+1)+j]);
      __m256d xi  = _mm256_set_m128d(xih,xil);
      __m256d xij, xijC;
      MM256_2MULTFMA(xi,r_[j],xij,xijC);
      __m256d xijC_ = _mm256_fmadd_pd(xi,rC_[j],xijC);
      __m256d oldS = s, t_;
      MM256_2SUM(oldS,xij,s,t_);
      __m256d tmp = _mm256_add_pd(xijC_,t_);
      t = _mm256_add_pd(t,tmp);
    }
    __m256d newY = _mm256_add_pd(s,t);
    _mm256_store_pd(&y[2*i],newY);
  }
}


// more accurate gemm product y <- alpha*x*m + beta*y AVX2 kernel for y of blocksize 1
void dgemm_sb_prec_k_strided_1(int nrows, int k, double alpha, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double beta, double *restrict y)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  if( !is_aligned(y,32) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
    exit(1);
    return;
  }
  // buffer rows of r and multiply by alpha
  __m256d r_[k], rC_[k];
  __m256d alpha_ = _mm256_set1_pd(alpha);
  for(int i = 0; i < k; i++)
  {
    __m256d tmp = _mm256_set1_pd(r[i]);
    __m256d tmpC = _mm256_set1_pd(rC[i]);
    __m256d pi;
    MM256_2MULTFMA(alpha_,tmp,r_[i],pi);
    rC_[i] = _mm256_fmadd_pd(alpha_,tmpC,pi);
  }
  __m256d beta_ = _mm256_set1_pd(beta);


#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i+=4)
  {
    __m256d oldY = _mm256_load_pd(&y[i]);
    __m256d s, t;
    MM256_2MULTFMA(beta_,oldY,s,t);
    for(int j = 0; j < k; j++)
    {
      __m256d xi = _mm256_set_pd(x[ldx*i+3*ldx+j],x[ldx*i+2*ldx+j],x[ldx*i+ldx+j],x[ldx*i+j]);
      __m256d xij, xijC;
      MM256_2MULTFMA(xi,r_[j],xij,xijC);
      __m256d xijC_ = _mm256_fmadd_pd(xi,rC_[j],xijC);
      __m256d oldS = s, t_;
      MM256_2SUM(oldS,xij,s,t_);
      __m256d tmp = _mm256_add_pd(xijC_,t_);
      t = _mm256_add_pd(t,tmp);
    }
    __m256d newY = _mm256_add_pd(s,t);
    _mm256_store_pd(&y[i],newY);
  }
}


// more accurate gemm product y <- alpha*x*m AVX2 kernel for y of blocksize 4 with non-temporal stores
void dgemm_sb_prec_k_strided_4_nt(int nrows, int k, double alpha, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double *restrict y)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  if( !is_aligned(y,32) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
    exit(1);
    return;
  }


  // buffer rows of r and multiply by alpha
  __m256d r_[k], rC_[k];
  __m256d alpha_ = _mm256_set1_pd(alpha);
  for(int i = 0; i < k; i++)
  {
    __m256d tmp = _mm256_set_pd(r[i+3*k],r[i+2*k],r[i+k],r[i]);
    __m256d tmpC = _mm256_set_pd(rC[i+3*k],rC[i+2*k],rC[i+k],rC[i]);
    __m256d pi;
    MM256_2MULTFMA(alpha_,tmp,r_[i],pi);
    rC_[i] = _mm256_fmadd_pd(alpha_,tmpC,pi);
  }

#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i++)
  {
    __m256d s, t;
    // j = 0
    __m256d xi = _mm256_broadcast_sd(&x[ldx*i+0]);
    __m256d pi;
    MM256_2MULTFMA(xi,r_[0],s,pi);
    t = _mm256_fmadd_pd(xi,rC_[0],pi);

    for(int j = 1; j < k; j++)
    {
      __m256d xi = _mm256_broadcast_sd(&x[ldx*i+j]);
      __m256d xij, xijC;
      MM256_2MULTFMA(xi,r_[j],xij,xijC);
      __m256d xijC_ = _mm256_fmadd_pd(xi,rC_[j],xijC);
      __m256d oldS = s, t_;
      MM256_2SUM(oldS,xij,s,t_);
      __m256d tmp = _mm256_add_pd(xijC_,t_);
      t = _mm256_add_pd(t,tmp);
    }
    __m256d newY = _mm256_add_pd(s,t);
    _mm256_stream_pd(&y[4*i],newY);
  }
}


// more accurate gemm product y <- alpha*x*m AVX2 kernel for y of blocksize 2 with non-temporal stores
void dgemm_sb_prec_k_strided_2_nt(int nrows, int k, double alpha, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double *restrict y)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  if( !is_aligned(y,32) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
    exit(1);
    return;
  }


  // buffer rows of r and multiply by alpha
  __m256d r_[k], rC_[k];
  __m256d alpha_ = _mm256_set1_pd(alpha);
  for(int i = 0; i < k; i++)
  {
    __m256d tmp = _mm256_set_pd(r[i+k],r[i],r[i+k],r[i]);
    __m256d tmpC = _mm256_set_pd(rC[i+k],rC[i],rC[i+k],rC[i]);
    __m256d pi;
    MM256_2MULTFMA(alpha_,tmp,r_[i],pi);
    rC_[i] = _mm256_fmadd_pd(alpha_,tmpC,pi);
  }

#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i+=2)
  {
    __m256d s, t;
    // j = 0
    __m128d xil = _mm_load1_pd(&x[ldx*i+0]);
    __m128d xih = _mm_load1_pd(&x[ldx*(i+1)+0]);
    __m256d xi  = _mm256_set_m128d(xih,xil);
    __m256d pi;
    MM256_2MULTFMA(xi,r_[0],s,pi);
    t = _mm256_fmadd_pd(xi,rC_[0],pi);

    for(int j = 1; j < k; j++)
    {
      xil = _mm_load1_pd(&x[ldx*i+j]);
      xih = _mm_load1_pd(&x[ldx*(i+1)+j]);
      xi  = _mm256_set_m128d(xih,xil);
      __m256d xij, xijC;
      MM256_2MULTFMA(xi,r_[j],xij,xijC);
      __m256d xijC_ = _mm256_fmadd_pd(xi,rC_[j],xijC);
      __m256d oldS = s, t_;
      MM256_2SUM(oldS,xij,s,t_);
      __m256d tmp = _mm256_add_pd(xijC_,t_);
      t = _mm256_add_pd(t,tmp);
    }
    __m256d newY = _mm256_add_pd(s,t);
    _mm256_stream_pd(&y[2*i],newY);
  }
}


// more accurate gemm product y <- alpha*x*m AVX2 kernel for y of blocksize 1 with non-temporal stores
void dgemm_sb_prec_k_strided_1_nt(int nrows, int k, double alpha, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double *restrict y)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  if( !is_aligned(y,32) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
    exit(1);
    return;
  }
  // buffer rows of r and multiply by alpha
  __m256d r_[k], rC_[k];
  __m256d alpha_ = _mm256_set1_pd(alpha);
  for(int i = 0; i < k; i++)
  {
    __m256d tmp = _mm256_set1_pd(r[i]);
    __m256d tmpC = _mm256_set1_pd(rC[i]);
    __m256d pi;
    MM256_2MULTFMA(alpha_,tmp,r_[i],pi);
    rC_[i] = _mm256_fmadd_pd(alpha_,tmpC,pi);
  }


#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i+=4)
  {
    __m256d s, t;
    // j = 0
    __m256d xi = _mm256_set_pd(x[ldx*i+3*ldx+0],x[ldx*i+2*ldx+0],x[ldx*i+ldx+0],x[ldx*i+0]);
    __m256d pi;
    MM256_2MULTFMA(xi,r_[0],s,pi);
    t = _mm256_fmadd_pd(xi,rC_[0],pi);

    for(int j = 1; j < k; j++)
    {
      xi = _mm256_set_pd(x[ldx*i+3*ldx+j],x[ldx*i+2*ldx+j],x[ldx*i+ldx+j],x[ldx*i+j]);
      __m256d xij, xijC;
      MM256_2MULTFMA(xi,r_[j],xij,xijC);
      __m256d xijC_ = _mm256_fmadd_pd(xi,rC_[j],xijC);
      __m256d oldS = s, t_;
      MM256_2SUM(oldS,xij,s,t_);
      __m256d tmp = _mm256_add_pd(xijC_,t_);
      t = _mm256_add_pd(t,tmp);
    }
    __m256d newY = _mm256_add_pd(s,t);
    _mm256_stream_pd(&y[i],newY);
  }
}



