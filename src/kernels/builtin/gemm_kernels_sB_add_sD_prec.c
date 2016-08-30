/*! \file gemm_kernels_sB_add_sD_prec.c
 * Fast parallel BLAS-gemm like functions augmented with the product of the result with high precision for different blocksizes for mvec_module
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
#include "prec_helpers.h"



// more accureate fused gemm products y <- x*r + y*d AVX2 kernel for y of blocksize 4
void dgemm_sb_add_sd_prec_strided_k_4(int nrows, int k, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double *restrict y, const double *restrict d, const double *restrict dC)
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
  if( !is_aligned(y,32) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
    exit(1);
    return;
  }
#ifdef PHIST_HAVE_OPENMP
  int nt = omp_get_max_threads();
#else
  int nt = 1;
#endif

  // buffer rows of r and d
  __m256d r_[k], rC_[k];
  for(int i = 0; i < k; i++)
  {
    r_[i]  = _mm256_set_pd(r [i+3*k],r [i+2*k],r [i+k],r [i]);
    rC_[i] = _mm256_set_pd(rC[i+3*k],rC[i+2*k],rC[i+k],rC[i]);
  }
  __m256d d_[4], dC_[4];
  for(int i = 0; i < 4; i++)
  {
    d_[i]  = _mm256_set_pd(d [i+3*4],d [i+2*4],d [i+4],d [i]);
    dC_[i] = _mm256_set_pd(dC[i+3*4],dC[i+2*4],dC[i+4],dC[i]);
  }

#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i++)
  {
    __m256d yi = _mm256_load_pd(&y[4*i]);
    __m256d s, t;
    // get y*d
    {
      // j = 0
      __m256d yi_ = _mm256_permute4x64_pd(yi,1*0+4*0+16*0+64*0);
      __m256d tmp;
      MM256_2MULTFMA(yi_,d_[0],s,tmp);
      t = _mm256_fmadd_pd(yi_,dC_[0],tmp);

      // j = 1
      yi_ = _mm256_permute4x64_pd(yi,1*1+4*1+16*1+64*1);
      __m256d yij, yijC;
      MM256_2MULTFMA(yi_,d_[1],yij,yijC);
      __m256d yijC_ = _mm256_fmadd_pd(yi_,dC_[1],yijC);
      __m256d oldS = s, t_;
      MM256_2SUM(oldS,yij,s,t_);
      tmp = _mm256_add_pd(yijC_,t_);
      t = _mm256_add_pd(t,tmp);

      // j = 2
      yi_ = _mm256_permute4x64_pd(yi,1*2+4*2+16*2+64*2);
      MM256_2MULTFMA(yi_,d_[2],yij,yijC);
      yijC_ = _mm256_fmadd_pd(yi_,dC_[2],yijC);
      oldS = s;
      MM256_2SUM(oldS,yij,s,t_);
      tmp = _mm256_add_pd(yijC_,t_);
      t = _mm256_add_pd(t,tmp);

      // j = 3
      yi_ = _mm256_permute4x64_pd(yi,1*3+4*3+16*3+64*3);
      MM256_2MULTFMA(yi_,d_[3],yij,yijC);
      yijC_ = _mm256_fmadd_pd(yi_,dC_[3],yijC);
      oldS = s;
      MM256_2SUM(oldS,yij,s,t_);
      tmp = _mm256_add_pd(yijC_,t_);
      t = _mm256_add_pd(t,tmp);
    }

    // add x*r
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
    yi = _mm256_add_pd(s,t);
    _mm256_store_pd(&y[4*i],yi);
  }
}

void dgemm_sb_add_sd_prec_k_4(int nrows, int k, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, const double *restrict d, const double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_add_sd_prec_strided_k_4(nrows, k, x, k, r, rC, y, d, dC);
}

void dgemm_sb_add_sd_prec_4_4(int nrows, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, const double *restrict d, const double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_add_sd_prec_strided_k_4(nrows, 4, x, 4, r, rC, y, d, dC);
}

void dgemm_sb_add_sd_prec_2_4(int nrows, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, const double *restrict d, const double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_add_sd_prec_strided_k_4(nrows, 2, x, 2, r, rC, y, d, dC);
}

void dgemm_sb_add_sd_prec_1_4(int nrows, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, const double *restrict d, const double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_add_sd_prec_strided_k_4(nrows, 1, x, 1, r, rC, y, d, dC);
}


// more accureate fused gemm products y <- x*r + y*d AVX2 kernel for y of blocksize 2
void dgemm_sb_add_sd_prec_strided_k_2(int nrows, int k, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double *restrict y, const double *restrict d, const double *restrict dC)
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
  if( !is_aligned(y,32) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
    exit(1);
    return;
  }
#ifdef PHIST_HAVE_OPENMP
  int nt = omp_get_max_threads();
#else
  int nt = 1;
#endif


  // buffer rows of r and d
  __m256d r_[k], rC_[k];
  for(int i = 0; i < k; i++)
  {
    r_[i]  = _mm256_set_pd(r[i+k],r[i],r[i+k],r[i]);
    rC_[i] = _mm256_set_pd(rC[i+k],rC[i],rC[i+k],rC[i]);
  }
  __m256d d_[2], dC_[2];
  for(int i = 0; i < 2; i++)
  {
    d_[i]  = _mm256_set_pd(d [i+2],d [i],d [i+2],d [i]);
    dC_[i] = _mm256_set_pd(dC[i+2],dC[i],dC[i+2],dC[i]);
  }

#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i+=2)
  {
    __m256d yi = _mm256_load_pd(&y[2*i]);
    __m256d s, t;
    {
      // j = 0
      __m256d yi_ = _mm256_permute_pd(yi,1*0+2*0+4*0+8*0);
      __m256d tmp;
      MM256_2MULTFMA(yi_,d_[0],s,tmp);
      t = _mm256_fmadd_pd(yi_,dC_[0],tmp);

      // j = 1
      yi_ = _mm256_permute_pd(yi,1*1+2*1+4*1+8*1);
      __m256d yij, yijC;
      MM256_2MULTFMA(yi_,d_[1],yij,yijC);
      __m256d yijC_ = _mm256_fmadd_pd(yi_,dC_[1],yijC);
      __m256d oldS = s, t_;
      MM256_2SUM(oldS,yij,s,t_);
      tmp = _mm256_add_pd(yijC_,t_);
      t = _mm256_add_pd(t,tmp);
    }

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
    yi = _mm256_add_pd(s,t);
    _mm256_store_pd(&y[2*i],yi);
  }
}

void dgemm_sb_add_sd_prec_k_2(int nrows, int k, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, const double *restrict d, const double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_add_sd_prec_strided_k_2(nrows, k, x, k, r, rC, y, d, dC);
}

void dgemm_sb_add_sd_prec_4_2(int nrows, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, const double *restrict d, const double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_add_sd_prec_strided_k_2(nrows, 4, x, 4, r, rC, y, d, dC);
}

void dgemm_sb_add_sd_prec_2_2(int nrows, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, const double *restrict d, const double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_add_sd_prec_strided_k_2(nrows, 2, x, 2, r, rC, y, d, dC);
}

void dgemm_sb_add_sd_prec_1_2(int nrows, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, const double *restrict d, const double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_add_sd_prec_strided_k_2(nrows, 1, x, 1, r, rC, y, d, dC);
}


// more accureate fused gemm products y <- x*r + y*d AVX2 kernel for y of blocksize 1
void dgemm_sb_add_sd_prec_strided_k_1(int nrows, int k, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double *restrict y, const double *restrict d, const double *restrict dC)
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
  if( !is_aligned(y,32) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
    exit(1);
    return;
  }
#ifdef PHIST_HAVE_OPENMP
  int nt = omp_get_max_threads();
#else
  int nt = 1;
#endif

  // buffer rows of r and d
  __m256d r_[k], rC_[k];
  for(int i = 0; i < k; i++)
  {
    r_[i]  = _mm256_set1_pd(r[i]);
    rC_[i] = _mm256_set1_pd(rC[i]);
  }
  __m256d d_ = _mm256_set1_pd(*d);
  __m256d dC_ = _mm256_set1_pd(*dC);

#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i+=4)
  {
    __m256d yi = _mm256_load_pd(&y[i]);
    __m256d s, t;
    __m256d tmp;
    MM256_2MULTFMA(yi,d_,s,tmp);
    t = _mm256_fmadd_pd(yi,dC_,tmp);

    for(int j = 0; j < k; j++)
    {
      __m256d xi = _mm256_set_pd(x[ldx*i+3*k+j],x[ldx*i+2*k+j],x[ldx*i+k+j],x[ldx*i+j]);
      __m256d xij, xijC;
      MM256_2MULTFMA(xi,r_[j],xij,xijC);
      __m256d xijC_ = _mm256_fmadd_pd(xi,rC_[j],xijC);
      __m256d oldS = s, t_;
      MM256_2SUM(oldS,xij,s,t_);
      __m256d tmp = _mm256_add_pd(xijC_,t_);
      t = _mm256_add_pd(t,tmp);
    }
    yi = _mm256_add_pd(s,t);
    _mm256_store_pd(&y[i],yi);
  }
}

void dgemm_sb_add_sd_prec_k_1(int nrows, int k, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, const double *restrict d, const double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_add_sd_prec_strided_k_1(nrows, k, x, k, r, rC, y, d, dC);
}

void dgemm_sb_add_sd_prec_4_1(int nrows, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, const double *restrict d, const double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_add_sd_prec_strided_k_1(nrows, 4, x, 4, r, rC, y, d, dC);
}

void dgemm_sb_add_sd_prec_2_1(int nrows, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, const double *restrict d, const double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_add_sd_prec_strided_k_1(nrows, 2, x, 2, r, rC, y, d, dC);
}

void dgemm_sb_add_sd_prec_1_1(int nrows, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, const double *restrict d, const double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_add_sd_prec_strided_k_1(nrows, 1, x, 1, r, rC, y, d, dC);
}

