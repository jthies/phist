/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
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
#include "phist_defs.h"
#include "phist_prec_helpers.h"

// define our own datatype for aligned doubles
typedef double aligned_double __attribute__((aligned(64)));


// more accurate gemm product x'x AVX2 kernel
void dgemm_sc_self_prec_4(int nrows, const aligned_double *restrict x, double *restrict res, double *restrict resC)
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
        __m256d xi = _mm256_load_pd(&x[4*i]);
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
void dgemm_sc_self_prec_2(int nrows, const aligned_double *restrict x, double *restrict res, double *restrict resC)
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
        __m256d xi = _mm256_load_pd(&x[4*i]);
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
void dgemm_sc_prec_4_k(int nrows, int k, const aligned_double *restrict x, const aligned_double *restrict y, double *restrict res, double *restrict resC)
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
        __m256d xi = _mm256_load_pd(&x[4*i]);
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

// more accurate gemm product x'y AVX2 kernel
void dgemm_sc_prec_4_4(int nrows, const aligned_double *restrict x, const aligned_double *restrict y, double *restrict res, double *restrict resC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sc_prec_4_k(nrows, 4, x, y, res, resC);
}

// more accurate gemm product x'y AVX2 kernel
void dgemm_sc_prec_4_2(int nrows, const aligned_double *restrict x, const aligned_double *restrict y, double *restrict res, double *restrict resC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sc_prec_4_k(nrows, 2, x, y, res, resC);
}

// more accurate gemm product x'y AVX2 kernel
void dgemm_sc_prec_4_1(int nrows, const aligned_double *restrict x, const aligned_double *restrict y, double *restrict res, double *restrict resC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sc_prec_4_k(nrows, 1, x, y, res, resC);
}


// more accurate gemm product x'y AVX2 kernel
void dgemm_sc_prec_2_k(int nrows, int k, const aligned_double *restrict x, const aligned_double *restrict y, double *restrict res, double *restrict resC)
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
        __m256d xi = _mm256_load_pd(&x[2*i]);
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

// more accurate gemm product x'y AVX2 kernel
void dgemm_sc_prec_2_2(int nrows, const aligned_double *restrict x, const aligned_double *restrict y, double *restrict res, double *restrict resC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sc_prec_2_k(nrows, 2, x, y, res, resC);
}

// more accurate gemm product x'y AVX2 kernel
void dgemm_sc_prec_2_1(int nrows, const aligned_double *restrict x, const aligned_double *restrict y, double *restrict res, double *restrict resC)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  dgemm_sc_prec_2_k(nrows, 1, x, y, res, resC);
}


// more accurate gemm product x'y AVX2 kernel
void dgemm_sc_prec_1_k(int nrows, int k, const aligned_double *restrict x, const aligned_double *restrict y, double *restrict res, double *restrict resC)
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
        __m256d xi = _mm256_load_pd(&x[i]);
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

