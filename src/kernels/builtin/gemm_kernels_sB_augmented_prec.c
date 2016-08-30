/*! \file gemm_kernels_sB_augmented_prec.c
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



// more accurate gemm product y <- alpha*x*m + beta*y AVX2 kernel for y of blocksize 4
void dgemm_sb_augmented_prec_strided_k_4(int nrows, int k, double alpha, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
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

  {
    // buffer for omp thread result + padding to prevent false sharing
    __m256d ds_[nt][8];
    __m256d dc_[nt][8];

#pragma omp parallel shared(ds_,dc_)
    {
      // initialize sum
      __m256d ds[4];
      __m256d dc[4];
      for(int j = 0; j < 4; j++)
      {
        ds[j] = _mm256_setzero_pd();
        dc[j] = _mm256_setzero_pd();
      }

#pragma omp for schedule(static)
      for(int i = 0; i < nrows; i++)
      {
        __m256d yi = _mm256_load_pd(&y[4*i]);
        __m256d s, t;
        MM256_2MULTFMA(beta_,yi,s,t);
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

        // multiply row yi with itself
        // j = 0
        __m256d yij = yi;
        MM256_4DOTADD(yi,yij,ds[0],dc[0]);

        // j = 1
        yij = _mm256_permute4x64_pd(yi,1*1+4*2+16*3+64*0);
        MM256_4DOTADD(yi,yij,ds[1],dc[1]);

        // j = 2
        yij = _mm256_permute4x64_pd(yi,1*2+4*3+16*0+64*1);
        MM256_4DOTADD(yi,yij,ds[2],dc[2]);
      }

      int it = omp_get_thread_num();
      for(int j = 0; j < 4; j++)
      {
        ds_[it][j] = ds[j];
        dc_[it][j] = dc[j];
      }
    }

    // handcoded omp reduction
    __m256d ds[4], dc[4];
    for(int j = 0; j < 4; j++)
    {
      ds[j] = _mm256_setzero_pd();
      dc[j] = _mm256_setzero_pd();
      for(int i = 0; i < nt; i++)
      {
        __m256d oldS = ds[j], oldC = dc[j];
        MM256_4SUM(oldS,oldC,ds_[i][j],dc_[i][j],ds[j],dc[j]);
      }
    }

    // obtain sums to construct result
    double d__[3][4], dC__[3][4];
    for(int j = 0; j < 3; j++)
    {
      _mm256_storeu_pd(d__[j],  ds[j]);
      _mm256_storeu_pd(dC__[j], dc[j]);
    }

    // construct and return result, needs to be summed up
    for(int i = 0; i < 3; i++)
    {
      for(int j = 0; j < 4; j++)
      {
        int j_ = (j+i)%4;
        d[j*4+j_] = d[j_*4+j] = d__[i][j];
        dC[j*4+j_] = dC[j_*4+j] = dC__[i][j];
      }
    }
  }
}

void dgemm_sb_augmented_prec_k_4(int nrows, int k, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_augmented_prec_strided_k_4(nrows, k, alpha, x, k, r, rC, beta, y, d, dC);
}

void dgemm_sb_augmented_prec_4_4(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_augmented_prec_strided_k_4(nrows, 4, alpha, x, 4, r, rC, beta, y, d, dC);
}

void dgemm_sb_augmented_prec_2_4(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_augmented_prec_strided_k_4(nrows, 2, alpha, x, 2, r, rC, beta, y, d, dC);
}

void dgemm_sb_augmented_prec_1_4(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_augmented_prec_strided_k_4(nrows, 1, alpha, x, 1, r, rC, beta, y, d, dC);
}


// more accurate gemm product y <- alpha*x*m + beta*y AVX2 kernel for y of blocksize 2
void dgemm_sb_augmented_prec_strided_k_2(int nrows, int k, double alpha, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
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

  {
    // buffer for omp thread result + padding to prevent false sharing
    __m256d ds_[nt][8];
    __m256d dc_[nt][8];

#pragma omp parallel shared(ds_,dc_)
    {
      // initialize sum
      __m256d ds[2];
      __m256d dc[2];
      for(int j = 0; j < 2; j++)
      {
        ds[j] = _mm256_setzero_pd();
        dc[j] = _mm256_setzero_pd();
      }

#pragma omp for schedule(static)
      for(int i = 0; i < nrows; i+=2)
      {
        __m256d yi = _mm256_load_pd(&y[2*i]);
        __m256d s, t;
        MM256_2MULTFMA(beta_,yi,s,t);
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

        // multiply new row yi with itself
        __m256d yij = _mm256_permute_pd(yi,0);
        MM256_4DOTADD(yi,yij,ds[0],dc[0]);

        yij = _mm256_permute_pd(yi,1+2+4+8);
        MM256_4DOTADD(yi,yij,ds[1],dc[1]);
      }
      int it = omp_get_thread_num();
      for(int j = 0; j < 2; j++)
      {
        ds_[it][j] = ds[j];
        dc_[it][j] = dc[j];
      }
    }


    // handcoded omp reduction
    __m256d ds[2], dc[2];
    for(int j = 0; j < 2; j++)
    {
      ds[j] = _mm256_setzero_pd();
      dc[j] = _mm256_setzero_pd();
      for(int i = 0; i < nt; i++)
      {
        __m256d oldS = ds[j], oldC = dc[j];
        MM256_4SUM(oldS,oldC,ds_[i][j],dc_[i][j],ds[j],dc[j]);
      }
    }

    // sum up 4 elements in mm256 to 2 elements in mm128
    // and return result, needs to be summed up again
    for(int j = 0; j < 2; j++)
    {
      __m128d ds2, dc2;
      MM256TO128_4SUM(ds[j],dc[j],ds2,dc2);
      _mm_storeu_pd(&d[j*2],  ds2);
      _mm_storeu_pd(&dC[j*2], dc2);
    }
  }

}

void dgemm_sb_augmented_prec_k_2(int nrows, int k, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_augmented_prec_strided_k_2(nrows, k, alpha, x, k, r, rC, beta, y, d, dC);
}

void dgemm_sb_augmented_prec_4_2(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_augmented_prec_strided_k_2(nrows, 4, alpha, x, 4, r, rC, beta, y, d, dC);
}

void dgemm_sb_augmented_prec_2_2(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_augmented_prec_strided_k_2(nrows, 2, alpha, x, 2, r, rC, beta, y, d, dC);
}

void dgemm_sb_augmented_prec_1_2(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_augmented_prec_strided_k_2(nrows, 1, alpha, x, 1, r, rC, beta, y, d, dC);
}


// more accurate gemm product y <- alpha*x*m + beta*y AVX2 kernel for y of blocksize 1
void dgemm_sb_augmented_prec_strided_k_1(int nrows, int k, double alpha, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
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

  {
    // buffer for omp thread result + padding to prevent false sharing
    __m256d ds_[nt][8];
    __m256d dc_[nt][8];

#pragma omp parallel shared(ds_,dc_)
    {
      // initialize sum
      __m256d ds = _mm256_setzero_pd();
      __m256d dc = _mm256_setzero_pd();


#pragma omp for schedule(static)
      for(int i = 0; i < nrows; i+=4)
      {
        __m256d yi = _mm256_load_pd(&y[i]);
        __m256d s, t;
        MM256_2MULTFMA(beta_,yi,s,t);
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

        // multiply new row yi with itself
        MM256_4DOTADD(yi,yi,ds,dc);
      }
      int it = omp_get_thread_num();
      ds_[it][0] = ds;
      dc_[it][0] = dc;
    }


    // handcoded omp reduction
    __m256d ds = _mm256_setzero_pd();
    __m256d dc = _mm256_setzero_pd();
    for(int i = 0; i < nt; i++)
    {
      __m256d oldS = ds, oldC = dc;
      MM256_4SUM(oldS,oldC,ds_[i][0],dc_[i][0],ds,dc);
    }

    // we still need to sum up s, c
    double s[4], c[4];
    _mm256_storeu_pd(s, ds);
    _mm256_storeu_pd(c, dc);
    prec_reduction_1(4, s, c, d, dC);
  }
}

void dgemm_sb_augmented_prec_k_1(int nrows, int k, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_augmented_prec_strided_k_1(nrows, k, alpha, x, k, r, rC, beta, y, d, dC);
}

void dgemm_sb_augmented_prec_4_1(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_augmented_prec_strided_k_1(nrows, 4, alpha, x, 4, r, rC, beta, y, d, dC);
}

void dgemm_sb_augmented_prec_2_1(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_augmented_prec_strided_k_1(nrows, 2, alpha, x, 2, r, rC, beta, y, d, dC);
}

void dgemm_sb_augmented_prec_1_1(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double beta, double *restrict y, double *restrict d, double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_augmented_prec_strided_k_1(nrows, 1, alpha, x, 1, r, rC, beta, y, d, dC);
}


// more accurate gemm product y <- alpha*x*m AVX2 kernel for y of blocksize 4 with non-temporal stores
void dgemm_sb_augmented_prec_strided_k_4_nt(int nrows, int k, double alpha, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double *restrict y, double *restrict d, double *restrict dC)
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

  {
    // buffer for omp thread result + padding to prevent false sharing
    __m256d ds_[nt][8];
    __m256d dc_[nt][8];

#pragma omp parallel shared(ds_,dc_)
    {
      // initialize sum
      __m256d ds[4];
      __m256d dc[4];
      for(int j = 0; j < 4; j++)
      {
        ds[j] = _mm256_setzero_pd();
        dc[j] = _mm256_setzero_pd();
      }


#pragma omp for schedule(static)
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
          __m256d xi = _mm256_broadcast_sd(&x[k*i+j]);
          __m256d xij, xijC;
          MM256_2MULTFMA(xi,r_[j],xij,xijC);
          __m256d xijC_ = _mm256_fmadd_pd(xi,rC_[j],xijC);
          __m256d oldS = s, t_;
          MM256_2SUM(oldS,xij,s,t_);
          __m256d tmp = _mm256_add_pd(xijC_,t_);
          t = _mm256_add_pd(t,tmp);
        }
        __m256d yi = _mm256_add_pd(s,t);
        _mm256_stream_pd(&y[4*i],yi);

        // multiply row yi with itself
        // j = 0
        __m256d yij = yi;
        MM256_4DOTADD(yi,yij,ds[0],dc[0]);

        // j = 1
        yij = _mm256_permute4x64_pd(yi,1*1+4*2+16*3+64*0);
        MM256_4DOTADD(yi,yij,ds[1],dc[1]);

        // j = 2
        yij = _mm256_permute4x64_pd(yi,1*2+4*3+16*0+64*1);
        MM256_4DOTADD(yi,yij,ds[2],dc[2]);
      }

      int it = omp_get_thread_num();
      for(int j = 0; j < 4; j++)
      {
        ds_[it][j] = ds[j];
        dc_[it][j] = dc[j];
      }
    }

    // handcoded omp reduction
    __m256d ds[4], dc[4];
    for(int j = 0; j < 4; j++)
    {
      ds[j] = _mm256_setzero_pd();
      dc[j] = _mm256_setzero_pd();
      for(int i = 0; i < nt; i++)
      {
        __m256d oldS = ds[j], oldC = dc[j];
        MM256_4SUM(oldS,oldC,ds_[i][j],dc_[i][j],ds[j],dc[j]);
      }
    }

    // obtain sums to construct result
    double d__[3][4], dC__[3][4];
    for(int j = 0; j < 3; j++)
    {
      _mm256_storeu_pd(d__[j],  ds[j]);
      _mm256_storeu_pd(dC__[j], dc[j]);
    }

    // construct and return result, needs to be summed up
    for(int i = 0; i < 3; i++)
    {
      for(int j = 0; j < 4; j++)
      {
        int j_ = (j+i)%4;
        d[j*4+j_] = d[j_*4+j] = d__[i][j];
        dC[j*4+j_] = dC[j_*4+j] = dC__[i][j];
      }
    }
  }

}

void dgemm_sb_augmented_prec_k_4_nt(int nrows, int k, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, double *restrict d, double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_augmented_prec_strided_k_4_nt(nrows, k, alpha, x, k, r, rC, y, d, dC);
}

void dgemm_sb_augmented_prec_4_4_nt(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, double *restrict d, double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_augmented_prec_strided_k_4_nt(nrows, 4, alpha, x, 4, r, rC, y, d, dC);
}

void dgemm_sb_augmented_prec_2_4_nt(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, double *restrict d, double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_augmented_prec_strided_k_4_nt(nrows, 2, alpha, x, 2, r, rC, y, d, dC);
}

void dgemm_sb_augmented_prec_1_4_nt(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, double *restrict d, double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_augmented_prec_strided_k_4_nt(nrows, 1, alpha, x, 1, r, rC, y, d, dC);
}


// more accurate gemm product y <- alpha*x*m AVX2 kernel for y of blocksize 2 with non-temporal stores
void dgemm_sb_augmented_prec_strided_k_2_nt(int nrows, int k, double alpha, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double *restrict y, double *restrict d, double *restrict dC)
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

  {
    // buffer for omp thread result + padding to prevent false sharing
    __m256d ds_[nt][8];
    __m256d dc_[nt][8];

#pragma omp parallel shared(ds_,dc_)
    {
      // initialize sum
      __m256d ds[2];
      __m256d dc[2];
      for(int j = 0; j < 2; j++)
      {
        ds[j] = _mm256_setzero_pd();
        dc[j] = _mm256_setzero_pd();
      }

#pragma omp for schedule(static)
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
        __m256d yi = _mm256_add_pd(s,t);
        _mm256_stream_pd(&y[2*i],yi);

        // multiply new row yi with itself
        __m256d yij = _mm256_permute_pd(yi,0);
        MM256_4DOTADD(yi,yij,ds[0],dc[0]);

        yij = _mm256_permute_pd(yi,1+2+4+8);
        MM256_4DOTADD(yi,yij,ds[1],dc[1]);
      }
      int it = omp_get_thread_num();
      for(int j = 0; j < 2; j++)
      {
        ds_[it][j] = ds[j];
        dc_[it][j] = dc[j];
      }
    }


    // handcoded omp reduction
    __m256d ds[2], dc[2];
    for(int j = 0; j < 2; j++)
    {
      ds[j] = _mm256_setzero_pd();
      dc[j] = _mm256_setzero_pd();
      for(int i = 0; i < nt; i++)
      {
        __m256d oldS = ds[j], oldC = dc[j];
        MM256_4SUM(oldS,oldC,ds_[i][j],dc_[i][j],ds[j],dc[j]);
      }
    }

    // sum up 4 elements in mm256 to 2 elements in mm128
    // and return result, needs to be summed up again
    for(int j = 0; j < 2; j++)
    {
      __m128d ds2, dc2;
      MM256TO128_4SUM(ds_[j][0],dc_[j][0],ds2,dc2);
      _mm_storeu_pd(&d[j*2],  ds2);
      _mm_storeu_pd(&dC[j*2], dc2);
    }
  }

}

void dgemm_sb_augmented_prec_k_2_nt(int nrows, int k, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, double *restrict d, double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_augmented_prec_strided_k_2_nt(nrows, k, alpha, x, k, r, rC, y, d, dC);
}

void dgemm_sb_augmented_prec_4_2_nt(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, double *restrict d, double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_augmented_prec_strided_k_2_nt(nrows, 4, alpha, x, 4, r, rC, y, d, dC);
}

void dgemm_sb_augmented_prec_2_2_nt(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, double *restrict d, double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_augmented_prec_strided_k_2_nt(nrows, 2, alpha, x, 2, r, rC, y, d, dC);
}

void dgemm_sb_augmented_prec_1_2_nt(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, double *restrict d, double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_augmented_prec_strided_k_2_nt(nrows, 1, alpha, x, 1, r, rC, y, d, dC);
}


// more accurate gemm product y <- alpha*x*m AVX2 kernel for y of blocksize 1 with non-temporal stores
void dgemm_sb_augmented_prec_strided_k_1_nt(int nrows, int k, double alpha, const double *restrict x, int ldx, const double *restrict r, const double *restrict rC, double *restrict y, double *restrict d, double *restrict dC)
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

  {
    // buffer for omp thread result + padding to prevent false sharing
    __m256d ds_[nt][8];
    __m256d dc_[nt][8];

#pragma omp parallel shared(ds_,dc_)
    {
      // initialize sum
      __m256d ds = _mm256_setzero_pd();
      __m256d dc = _mm256_setzero_pd();


#pragma omp for schedule(static)
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
          __m256d xi = _mm256_set_pd(x[ldx*i+3*k+j],x[ldx*i+2*k+j],x[ldx*i+k+j],x[ldx*i+j]);
          __m256d xij, xijC;
          MM256_2MULTFMA(xi,r_[j],xij,xijC);
          __m256d xijC_ = _mm256_fmadd_pd(xi,rC_[j],xijC);
          __m256d oldS = s, t_;
          MM256_2SUM(oldS,xij,s,t_);
          __m256d tmp = _mm256_add_pd(xijC_,t_);
          t = _mm256_add_pd(t,tmp);
        }
        __m256d yi = _mm256_add_pd(s,t);
        _mm256_stream_pd(&y[i],yi);

        // multiply new row yi with itself
        MM256_4DOTADD(yi,yi,ds,dc);
      }
      int it = omp_get_thread_num();
      ds_[it][0] = ds;
      dc_[it][0] = dc;
    }


    // handcoded omp reduction
    __m256d ds = _mm256_setzero_pd();
    __m256d dc = _mm256_setzero_pd();
    for(int i = 0; i < nt; i++)
    {
      __m256d oldS = ds, oldC = dc;
      MM256_4SUM(oldS,oldC,ds_[i][0],dc_[i][0],ds,dc);
    }

    // we still need to sum up s, c
    double ds1[4], dc1[4];
    _mm256_storeu_pd(ds1, ds);
    _mm256_storeu_pd(dc1, dc);
    prec_reduction_1(4, ds1, dc1, d, dC);
  }
}

void dgemm_sb_augmented_prec_k_1_nt(int nrows, int k, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, double *restrict d, double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_augmented_prec_strided_k_1_nt(nrows, k, alpha, x, k, r, rC, y, d, dC);
}

void dgemm_sb_augmented_prec_4_1_nt(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, double *restrict d, double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_augmented_prec_strided_k_1_nt(nrows, 4, alpha, x, 4, r, rC, y, d, dC);
}

void dgemm_sb_augmented_prec_2_1_nt(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, double *restrict d, double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_augmented_prec_strided_k_1_nt(nrows, 2, alpha, x, 2, r, rC, y, d, dC);
}

void dgemm_sb_augmented_prec_1_1_nt(int nrows, double alpha, const double *restrict x, const double *restrict r, const double *restrict rC, double *restrict y, double *restrict d, double *restrict dC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sb_augmented_prec_strided_k_1_nt(nrows, 1, alpha, x, 1, r, rC, y, d, dC);
}


