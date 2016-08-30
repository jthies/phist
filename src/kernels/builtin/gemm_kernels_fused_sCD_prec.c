/*! \file gemm_kernels_sC_prec.c
 * Fast parallel fused BLAS-gemm functions with high precision for different blocksizes for mvec_module
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

// define our own datatype for aligned doubles
typedef double aligned_double __attribute__((aligned(64)));


// more accurate gemm product x'x AVX2 kernel
void dgemm_fused_scd_self_prec_4(int nrows, aligned_double *restrict x, double *restrict scal, double *restrict scalC, double *restrict res, double *restrict resC)
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

#ifdef PHIST_HAVE_OPENMP
  int nt = omp_get_max_threads();
#else
  int nt = 1;
#endif

  // buffer rows of scal
  __m256d scal_[4], scalC_[4];
  for(int i = 0; i < 4; i++)
  {
    scal_[i] = _mm256_set_pd(scal[i+12],scal[i+8],scal[i+4],scal[i+0]);
    scalC_[i] = _mm256_set_pd(scalC[i+12],scalC[i+8],scalC[i+4],scalC[i+0]);
  }

  {
    // buffer for omp thread result
    __m256d s_[nt][8];
    __m256d c_[nt][8];

#pragma omp parallel shared(s_,c_)
    {
      __m256d s[4];
      __m256d c[4];
      // initialize sum
      for(int j = 0; j < 4; j++)
      {
        s[j] = _mm256_setzero_pd();
        c[j] = _mm256_setzero_pd();
      }

#pragma omp for schedule(static)
      for(int i = 0; i < nrows; i++)
      {
        // multiplication of x with scal, scalC
        __m256d xi = _mm256_load_pd(&x[4*i]);
        {
          __m256d s, t;
          __m256d p, pi;
          // unroll j = 0
          __m256d xi_ = _mm256_permute4x64_pd(xi,0);
          MM256_2MULTFMA(xi_,scal_[0],p,pi);
          __m256d pi_ = _mm256_fmadd_pd(xi_,scalC_[0],pi);
          s = p, t = pi_;

          // j = 1
          xi_ = _mm256_permute4x64_pd(xi,1+4*1+16*1+64*1);
          MM256_2MULTFMA(xi_,scal_[1],p,pi);
          pi_ = _mm256_fmadd_pd(xi_,scalC_[1],pi);
          __m256d oldS = s, sigma;
          MM256_2SUM(oldS,p,s,sigma);
          __m256d tmp = _mm256_add_pd(pi_,sigma);
          t = _mm256_add_pd(t, tmp);

          // j = 2
          xi_ = _mm256_permute4x64_pd(xi,2+4*2+16*2+64*2);
          MM256_2MULTFMA(xi_,scal_[2],p,pi);
          pi_ = _mm256_fmadd_pd(xi_,scalC_[2],pi);
          oldS = s;
          MM256_2SUM(oldS,p,s,sigma); 
          tmp = _mm256_add_pd(pi_,sigma);
          t = _mm256_add_pd(t, tmp);

          // j = 3
          xi_ = _mm256_permute4x64_pd(xi,3+4*3+16*3+64*3);
          MM256_2MULTFMA(xi_,scal_[3],p,pi);
          pi_ = _mm256_fmadd_pd(xi_,scalC_[3],pi);
          oldS = s;
          MM256_2SUM(oldS,p,s,sigma);
          tmp = _mm256_add_pd(pi_,sigma);
          t = _mm256_add_pd(t,tmp);

          xi = _mm256_add_pd(s,t);
          _mm256_store_pd(&x[4*i],xi);
        }

        // multiply row xi with itself
        // j = 0
        __m256d xij = xi;
        MM256_4DOTADD(xi,xij,s[0],c[0]);

        // j = 1
        xij = _mm256_permute4x64_pd(xi,1*1+4*2+16*3+64*0);
        MM256_4DOTADD(xi,xij,s[1],c[1]);

        // j = 2
        xij = _mm256_permute4x64_pd(xi,1*2+4*3+16*0+64*1);
        MM256_4DOTADD(xi,xij,s[2],c[2]);
      }

      int it = omp_get_thread_num();
      // store result
      for(int j = 0; j < 4; j++)
      {
        s_[it][j] = s[j];
        c_[it][j] = c[j];
      }
    }


    // handcoded omp reduction
    __m256d s[4], c[4];
    for(int j = 0; j < 4; j++)
    {
      s[j] = _mm256_setzero_pd();
      c[j] = _mm256_setzero_pd();
      for(int i = 0; i < nt; i++)
      {
        __m256d oldS = s[j], oldC = c[j];
        MM256_4SUM(oldS,oldC,s_[i][j],c_[i][j],s[j],c[j]);
      }
    }

    // obtain sums to construct result
    double r[3][4], rC[3][4];
    for(int j = 0; j < 3; j++)
    {
      _mm256_storeu_pd(r[j],  s[j]);
      _mm256_storeu_pd(rC[j], c[j]);
    }

    // construct and return result, needs to be summed up
    for(int i = 0; i < 3; i++)
    {
      for(int j = 0; j < 4; j++)
      {
        int j_ = (j+i)%4;
        res[j*4+j_] = res[j_*4+j] = r[i][j];
        resC[j*4+j_] = resC[j_*4+j] = rC[i][j];
      }
    }
  }

}


// more accurate gemm product x'x AVX2 kernel
void dgemm_fused_scd_self_prec_2(int nrows, aligned_double *restrict x, const double *restrict scal, const double *restrict scalC, double *restrict res, double *restrict resC)
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


#ifdef PHIST_HAVE_OPENMP
  int nt = omp_get_max_threads();
#else
  int nt = 1;
#endif

  // buffer rows of scal twice
  __m256d scal_[2], scalC_[2];
  for(int i = 0; i < 2; i++)
  {
    scal_[i] = _mm256_set_pd(scal[i+2],scal[i+0],scal[i+2],scal[i+0]);
    scalC_[i] = _mm256_set_pd(scalC[i+2],scalC[i+0],scalC[i+2],scalC[i+0]);
  }

  {
    // buffer for omp thread result + padding to prevent false sharing
    __m256d s_[nt][8];
    __m256d c_[nt][8];

#pragma omp parallel shared(s_,c_)
    {
      // initialize sum
      __m256d s[2];
      __m256d c[2];
      for(int j = 0; j < 2; j++)
      {
        s[j] = _mm256_setzero_pd();
        c[j] = _mm256_setzero_pd();
      }

      int nrows2 = nrows/2;
#pragma omp for schedule(static)
      for(int i = 0; i < nrows2; i++)
      {
        // multiplication of x with scal (two rows at once)
        __m256d xi = _mm256_load_pd(&x[4*i]);
        {
          __m256d s, t;
          __m256d p, pi;
          // unroll j = 0
          __m256d xi_ = _mm256_permute_pd(xi,0);
          MM256_2MULTFMA(xi_,scal_[0],p,pi);
          __m256d pi_ = _mm256_fmadd_pd(xi_,scalC_[0],pi);
          s = p, t = pi_;

          // j = 1
          xi_ = _mm256_permute_pd(xi,1+2+4+8);
          MM256_2MULTFMA(xi_,scal_[1],p,pi);
          pi_ = _mm256_fmadd_pd(xi_,scalC_[1],pi);
          __m256d oldS = s, sigma;
          MM256_2SUM(oldS,p,s,sigma);
          __m256d tmp = _mm256_add_pd(pi_,sigma);
          t = _mm256_add_pd(t,tmp);

          xi = _mm256_add_pd(s,t);
          _mm256_store_pd(&x[4*i],xi);
        }

        // multiply new row xi with itself
        __m256d xij = _mm256_permute_pd(xi,0);
        MM256_4DOTADD(xi,xij,s[0],c[0]);

        xij = _mm256_permute_pd(xi,1+2+4+8);
        MM256_4DOTADD(xi,xij,s[1],c[1]);
      }

      int it = omp_get_thread_num();
      for(int j = 0; j < 2; j++)
      {
        s_[it][j] = s[j];
        c_[it][j] = c[j];
      }
    }


    // handcoded omp reduction
    __m256d s[2], c[2];
    for(int j = 0; j < 2; j++)
    {
      s[j] = _mm256_setzero_pd();
      c[j] = _mm256_setzero_pd();
      for(int i = 0; i < nt; i++)
      {
        __m256d oldS = s[j], oldC = c[j];
        MM256_4SUM(oldS,oldC,s_[i][j],c_[i][j],s[j],c[j]);
      }
    }

    // sum up 4 elements in mm256 to 2 elements in mm128
    // and return result, needs to be summed up again
    for(int j = 0; j < 2; j++)
    {
      __m128d s2, c2;
      MM256TO128_4SUM(s[j],c[j],s2,c2);
      _mm_storeu_pd(&res[j*2],  s2);
      _mm_storeu_pd(&resC[j*2], c2);
    }
  }

}


// more accurate gemm product x'y AVX2 kernel
void dgemm_fused_scd_prec_k_4(int nrows, int k, const aligned_double *restrict y, aligned_double *restrict x, const double *restrict scal, const double *restrict scalC, double *restrict res, double *restrict resC)
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

  // buffer rows of scal
  __m256d scal_[4], scalC_[4];
  for(int i = 0; i < 4; i++)
  {
    scal_[i] = _mm256_set_pd(scal[i+12],scal[i+8],scal[i+4],scal[i+0]);
    scalC_[i] = _mm256_set_pd(scalC[i+12],scalC[i+8],scalC[i+4],scalC[i+0]);
  }

  {
    // buffer for omp thread result + padding to prevent false sharing
    __m256d s_[nt][(k/8+1)*8];
    __m256d c_[nt][(k/8+1)*8];

#pragma omp parallel shared(s_,c_)
    {
      // initialize sum
      __m256d s[k];
      __m256d c[k];
      for(int j = 0; j < k; j++)
      {
        s[j] = _mm256_setzero_pd();
        c[j] = _mm256_setzero_pd();
      }

#pragma omp for schedule(static)
      for(int i = 0; i < nrows; i++)
      {
        // multiplication of x with scal, scalC
        __m256d xi = _mm256_load_pd(&x[4*i]);
        {
          __m256d s, t;
          __m256d p, pi;
          // unroll j = 0
          __m256d xi_ = _mm256_permute4x64_pd(xi,0);
          MM256_2MULTFMA(xi_,scal_[0],p,pi);
          __m256d pi_ = _mm256_fmadd_pd(xi_,scalC_[0],pi);
          s = p, t = pi_;

          // j = 1
          xi_ = _mm256_permute4x64_pd(xi,1+4*1+16*1+64*1);
          MM256_2MULTFMA(xi_,scal_[1],p,pi);
          pi_ = _mm256_fmadd_pd(xi_,scalC_[1],pi);
          __m256d oldS = s, sigma;
          MM256_2SUM(oldS,p,s,sigma);
          __m256d tmp = _mm256_add_pd(pi_,sigma);
          t = _mm256_add_pd(t, tmp);

          // j = 2
          xi_ = _mm256_permute4x64_pd(xi,2+4*2+16*2+64*2);
          MM256_2MULTFMA(xi_,scal_[2],p,pi);
          pi_ = _mm256_fmadd_pd(xi_,scalC_[2],pi);
          oldS = s;
          MM256_2SUM(oldS,p,s,sigma);
          tmp = _mm256_add_pd(pi_,sigma);
          t = _mm256_add_pd(t, tmp);

          // j = 3
          xi_ = _mm256_permute4x64_pd(xi,3+4*3+16*3+64*3);
          MM256_2MULTFMA(xi_,scal_[3],p,pi);
          pi_ = _mm256_fmadd_pd(xi_,scalC_[3],pi);
          oldS = s;
          MM256_2SUM(oldS,p,s,sigma);
          tmp = _mm256_add_pd(pi_,sigma);
          t = _mm256_add_pd(t,tmp);

          xi = _mm256_add_pd(s,t);
          _mm256_store_pd(&x[4*i],xi);
        }

        // multiply row xi with row yi
        for(int j = 0; j < k; j++)
        {
          __m256d yij = _mm256_broadcast_sd(&y[k*i+j]);
          MM256_4DOTADD(xi,yij,s[j],c[j]);
        }
      }

      int it = omp_get_thread_num();
      for(int j = 0; j < k; j++)
      {
        s_[it][j] = s[j];
        c_[it][j] = c[j];
      }
    }


    // handcoded omp reduction
    __m256d s[k], c[k];
    for(int j = 0; j < k; j++)
    {
      s[j] = _mm256_setzero_pd();
      c[j] = _mm256_setzero_pd();
      for(int i = 0; i < nt; i++)
      {
        __m256d oldS = s[j], oldC = c[j];
        MM256_4SUM(oldS,oldC,s_[i][j],c_[i][j],s[j],c[j]);
      }
    }

    // return result, needs to be summed up
    for(int j = 0; j < k; j++)
    {
      _mm256_storeu_pd(&res[j*4],  s[j]);
      _mm256_storeu_pd(&resC[j*4], c[j]);
    }
  }

}


void dgemm_fused_scd_prec_4_4(int nrows, const aligned_double *restrict y, aligned_double *restrict x, const double *restrict scal, const double *restrict scalC, double *restrict res, double *restrict resC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_fused_scd_prec_k_4(nrows, 4, y, x, scal, scalC, res, resC);
}

void dgemm_fused_scd_prec_2_4(int nrows, const aligned_double *restrict y, aligned_double *restrict x, const double *restrict scal, const double *restrict scalC, double *restrict res, double *restrict resC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_fused_scd_prec_k_4(nrows, 2, y, x, scal, scalC, res, resC);
}

void dgemm_fused_scd_prec_1_4(int nrows, const aligned_double *restrict y, aligned_double *restrict x, const double *restrict scal, const double *restrict scalC, double *restrict res, double *restrict resC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_fused_scd_prec_k_4(nrows, 1, y, x, scal, scalC, res, resC);
}


// more accurate gemm product x'y AVX2 kernel
void dgemm_fused_scd_prec_k_2(int nrows, int k, const aligned_double *restrict y, aligned_double *restrict x, const double *restrict scal, const double *restrict scalC, double *restrict res, double *restrict resC)
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

  // buffer rows of scal
  __m256d scal_[2], scalC_[2];
  for(int i = 0; i < 2; i++)
  {
    scal_[i] = _mm256_set_pd(scal[i+2],scal[i+0],scal[i+2],scal[i+0]);
    scalC_[i] = _mm256_set_pd(scalC[i+2],scalC[i+0],scalC[i+2],scalC[i+0]);
  }

  {
    // buffer for omp thread result + padding to prevent false sharing
    __m256d s_[nt][(k/8+1)*8];
    __m256d c_[nt][(k/8+1)*8];

#pragma omp parallel shared(s_,c_)
    {
      // initialize sum
      __m256d s[k];
      __m256d c[k];
      for(int j = 0; j < k; j++)
      {
        s[j] = _mm256_setzero_pd();
        c[j] = _mm256_setzero_pd();
      }

#pragma omp for schedule(static)
      for(int i = 0; i < nrows; i+=2)
      {
        // multiplication of x with scal (two rows at once)
        __m256d xi = _mm256_load_pd(&x[2*i]);
        {
          __m256d s, t;
          __m256d p, pi;
          // unroll j = 0
          __m256d xi_ = _mm256_permute_pd(xi,0);
          MM256_2MULTFMA(xi_,scal_[0],p,pi);
          __m256d pi_ = _mm256_fmadd_pd(xi_,scalC_[0],pi);
          s = p, t = pi_;

          // j = 1
          xi_ = _mm256_permute_pd(xi,1+2+4+8);
          MM256_2MULTFMA(xi_,scal_[1],p,pi);
          pi_ = _mm256_fmadd_pd(xi_,scalC_[1],pi);
          __m256d oldS = s, sigma;
          MM256_2SUM(oldS,p,s,sigma);
          __m256d tmp = _mm256_add_pd(pi_,sigma);
          t = _mm256_add_pd(t,tmp);

          xi = _mm256_add_pd(s,t);
          _mm256_store_pd(&x[2*i],xi);
        }

        // multiply new row xi with y
        for(int j = 0; j < k; j++)
        {
          __m128d yil = _mm_load1_pd(&y[k*i+j]);
          __m128d yih = _mm_load1_pd(&y[k*(i+1)+j]);
          __m256d yij = _mm256_set_m128d(yih,yil);
          MM256_4DOTADD(xi,yij,s[j],c[j]);
        }
      }

      int it = omp_get_thread_num();
      for(int j = 0; j < k; j++)
      {
        s_[it][j] = s[j];
        c_[it][j] = c[j];
      }
    }


    // handcoded omp reduction
    __m256d s[k], c[k];
    for(int j = 0; j < k; j++)
    {
      s[j] = _mm256_setzero_pd();
      c[j] = _mm256_setzero_pd();
      for(int i = 0; i < nt; i++)
      {
        __m256d oldS = s[j], oldC = c[j];
        MM256_4SUM(oldS,oldC,s_[i][j],c_[i][j],s[j],c[j]);
      }
    }

    // sum up 4 elements in mm256 to 2 elements in mm128
    // and return result, needs to be summed up again
    for(int j = 0; j < k; j++)
    {
      __m128d s2, c2;
      MM256TO128_4SUM(s[j],c[j],s2,c2);
      _mm_storeu_pd(&res[j*2],  s2);
      _mm_storeu_pd(&resC[j*2], c2);
    }
  }

}


void dgemm_fused_scd_prec_4_2(int nrows, const aligned_double *restrict y, aligned_double *restrict x, const double *restrict scal, const double *restrict scalC, double *restrict res, double *restrict resC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_fused_scd_prec_k_2(nrows, 4, y, x, scal, scalC, res, resC);
}

void dgemm_fused_scd_prec_2_2(int nrows, const aligned_double *restrict y, aligned_double *restrict x, const double *restrict scal, const double *restrict scalC, double *restrict res, double *restrict resC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_fused_scd_prec_k_2(nrows, 2, y, x, scal, scalC, res, resC);
}

void dgemm_fused_scd_prec_1_2(int nrows, const aligned_double *restrict y, aligned_double *restrict x, const double *restrict scal, const double *restrict scalC, double *restrict res, double *restrict resC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_fused_scd_prec_k_2(nrows, 1, y, x, scal, scalC, res, resC);
}


// more accurate gemm product x'y AVX2 kernel
void dgemm_fused_scd_prec_k_1(int nrows, int k, const aligned_double *restrict y, aligned_double *restrict x, const double *restrict scal, const double *restrict scalC, double *restrict res, double *restrict resC)
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

  // buffer rows of scal
  __m256d scal_ = _mm256_broadcast_sd(scal);
  __m256d scalC_ = _mm256_broadcast_sd(scalC);
  {
    // buffer for omp thread result + padding to prevent false sharing
    __m256d s_[nt][(k/8+1)*8];
    __m256d c_[nt][(k/8+1)*8];

#pragma omp parallel shared(s_,c_)
    {
      // initialize sum
      __m256d s[k];
      __m256d c[k];
      for(int j = 0; j < k; j++)
      {
        s[j] = _mm256_setzero_pd();
        c[j] = _mm256_setzero_pd();
      }

#pragma omp for schedule(static)
      for(int i = 0; i < nrows; i+=4)
      {
        // multiplication of x with scal (two rows at once)
        __m256d xi = _mm256_load_pd(&x[i]);
        {
          __m256d xir, xirC;
          MM256_2MULTFMA(xi,scal_,xir,xirC);
          __m256d tmp = _mm256_fmadd_pd(xi,scalC_,xirC);
          xi = _mm256_add_pd(xir,tmp);
          _mm256_store_pd(&x[i],xi);
        }

        // multiply new row xi with itself
        for(int j = 0; j < k; j++)
        {
          __m256d yij = _mm256_set_pd(y[k*(i+3)+j],y[k*(i+2)+j],y[k*(i+1)+j],y[k*(i+0)+j]);
          MM256_4DOTADD(xi,yij,s[j],c[j]);
        }
      }

      int it = omp_get_thread_num();
      for(int j = 0; j < k; j++)
      {
        s_[it][j] = s[j];
        c_[it][j] = c[j];
      }
    }



    // handcoded omp reduction
    __m256d s[k], c[k];
    for(int j = 0; j < k; j++)
    {
      s[j] = _mm256_setzero_pd();
      c[j] = _mm256_setzero_pd();
      for(int i = 0; i < nt; i++)
      {
        __m256d oldS = s[j], oldC = c[j];
        MM256_4SUM(oldS,oldC,s_[i][j],c_[i][j],s[j],c[j]);
      }
    }

    // sum up 4 elements in mm256 to 1 double
    // and store result in res,resC
    for(int j = 0; j < k; j++)
    {
      double sj[4], cj[4];
      _mm256_storeu_pd(sj, s[j]);
      _mm256_storeu_pd(cj, c[j]);
      prec_reduction_1(4, sj, cj, &res[j], &resC[j]);
    }
  }

}


void dgemm_fused_scd_prec_4_1(int nrows, const aligned_double *restrict y, aligned_double *restrict x, const double *restrict scal, const double *restrict scalC, double *restrict res, double *restrict resC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_fused_scd_prec_k_1(nrows, 4, y, x, scal, scalC, res, resC);
}

void dgemm_fused_scd_prec_2_1(int nrows, const aligned_double *restrict y, aligned_double *restrict x, const double *restrict scal, const double *restrict scalC, double *restrict res, double *restrict resC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_fused_scd_prec_k_1(nrows, 2, y, x, scal, scalC, res, resC);
}

