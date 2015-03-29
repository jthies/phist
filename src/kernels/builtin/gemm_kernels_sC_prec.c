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


// more accurate gemm product x'x AVX2 kernel
void dgemm_sc_self_prec_4(int nrows, const double *restrict x, double *restrict res, double *restrict resC)
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

#ifdef PHIST_HAVE_OPENMP
  int nt = omp_get_max_threads();
#else
  int nt = 1;
#endif

  {
    // buffer for omp thread result + padding to prevent false sharing
    __m256d s_[8][nt];
    __m256d c_[8][nt];

#pragma omp parallel shared(s_,c_)
    {
      // initialize sum
      __m256d s[4];
      __m256d c[4];
      for(int j = 0; j < 4; j++)
      {
        s[j] = _mm256_setzero_pd();
        c[j] = _mm256_setzero_pd();
      }

#pragma omp for schedule(static)
      for(int i = 0; i < nrows; i++)
      {
        __m256d xi = _mm256_load_pd(&x[4*i]);
        for(int j = 0; j < 4; j++)
        {
          __m256d xij = _mm256_broadcast_sd(&x[4*i+j]);
          __m256d p, pi;
          MM256_2MULTFMA(xi,xij,p,pi);
          __m256d sigma, oldS = s[j];
          MM256_FAST2SUM(oldS,p, s[j],sigma);
          //MM256_2SUM(oldS,p, s[j],sigma); // more FP ops than Kahan-style FAST2SUM, but exacter
          __m256d tmp = _mm256_add_pd(pi,sigma);
          c[j] = _mm256_add_pd(c[j],tmp);
        }
      }

      int it = omp_get_thread_num();
      for(int j = 0; j < 4; j++)
      {
        s_[j][it] = s[j];
        c_[j][it] = c[j];
      }
    }


    // handcoded omp reduction
    for(int i = 1; i < nt; i++)
    {
      for(int j = 0; j < 4; j++)
      {
        __m256d sigma, oldS = s_[j][0];
        MM256_FAST2SUM(oldS, s_[j][i], s_[j][0], sigma);
        //MM256_2SUM(oldS, s_[0][i], s_[j][0], sigma);
        __m256d tmp = _mm256_add_pd(c_[j][i],sigma);
        c_[j][0] = _mm256_add_pd(c_[j][0], tmp);
      }
    }

    // return result, needs to be summed up
    for(int j = 0; j < 4; j++)
    {
      _mm256_storeu_pd(&res[j*4],  s_[j][0]);
      _mm256_storeu_pd(&resC[j*4], c_[j][0]);
    }
  }

}


// more accurate gemm product x'x AVX2 kernel
void dgemm_sc_self_prec_2(int nrows, const double *restrict x, double *restrict res, double *restrict resC)
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

#ifdef PHIST_HAVE_OPENMP
  int nt = omp_get_max_threads();
#else
  int nt = 1;
#endif

  {
    // buffer for omp thread result + padding to prevent false sharing
    __m256d s_[8][nt];
    __m256d c_[8][nt];

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

#pragma omp for schedule(static)
      for(int i = 0; i < nrows; i++)
      {
        __m256d xi = _mm256_load_pd(&x[4*i]);
        for(int j = 0; j < 2; j++)
        {
          __m256d xij = _mm256_set_pd(x[4*i+j+2],x[4*i+j+2],x[4*i+j],x[4*i+j]);
          __m256d p, pi;
          MM256_2MULTFMA(xi,xij,p,pi);
          __m256d sigma, oldS = s[j];
          MM256_FAST2SUM(oldS,p, s[j],sigma);
          //MM256_2SUM(oldS,p, s[j],sigma); // more FP ops than Kahan-style FAST2SUM, but exacter
          __m256d tmp = _mm256_add_pd(pi,sigma);
          c[j] = _mm256_add_pd(c[j],tmp);
        }
      }

      int it = omp_get_thread_num();
      for(int j = 0; j < 2; j++)
      {
        s_[j][it] = s[j];
        c_[j][it] = c[j];
      }
    }


    // handcoded omp reduction
    for(int i = 1; i < nt; i++)
    {
      for(int j = 0; j < 2; j++)
      {
        __m256d sigma, oldS = s_[j][0];
        MM256_FAST2SUM(oldS, s_[j][i], s_[j][0], sigma);
        //MM256_2SUM(oldS, s_[0][i], s_[j][0], sigma);
        __m256d tmp = _mm256_add_pd(c_[j][i],sigma);
        c_[j][0] = _mm256_add_pd(c_[j][0], tmp);
      }
    }


    // sum up 4 elements in mm256 to 2 elements in mm128
    // and return result, needs to be summed up again
    for(int j = 0; j < 2; j++)
    {
      __m128d s = _mm256_extractf128_pd(s_[j][0],0);
      __m128d p = _mm256_extractf128_pd(s_[j][0],1);
      __m128d c = _mm256_extractf128_pd(c_[j][0],0);
      __m128d pi = _mm256_extractf128_pd(c_[j][0],1);
      __m128d sigma, oldS = s;
      MM128_FAST2SUM(oldS,p,s,sigma);
      __m128d tmp = _mm_add_pd(pi,sigma);
      c = _mm_add_pd(c,tmp);
      _mm_storeu_pd(&res[j*2],  s);
      _mm_storeu_pd(&resC[j*2], c);
    }
  }

}

// more accurate gemm product x'y AVX2 kernel
void dgemm_sc_prec_4_4(int nrows, const double *restrict x, const double *restrict y, double *restrict res, double *restrict resC)
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
#ifdef PHIST_HAVE_OPENMP
  int nt = omp_get_max_threads();
#else
  int nt = 1;
#endif

  {
    // buffer for omp thread result + padding to prevent false sharing
    __m256d s_[8][nt];
    __m256d c_[8][nt];

#pragma omp parallel shared(s_,c_)
    {
      // initialize sum
      __m256d s[4];
      __m256d c[4];
      for(int j = 0; j < 4; j++)
      {
        s[j] = _mm256_setzero_pd();
        c[j] = _mm256_setzero_pd();
      }

#pragma omp for schedule(static)
      for(int i = 0; i < nrows; i++)
      {
        __m256d xi = _mm256_load_pd(&x[4*i]);
        for(int j = 0; j < 4; j++)
        {
          __m256d yij = _mm256_broadcast_sd(&y[4*i+j]);
          __m256d p, pi;
          MM256_2MULTFMA(xi,yij,p,pi);
          __m256d sigma, oldS = s[j];
          MM256_FAST2SUM(oldS,p, s[j],sigma);
          //MM256_2SUM(oldS,p, s[j],sigma); // more FP ops than Kahan-style FAST2SUM, but exacter
          __m256d tmp = _mm256_add_pd(pi,sigma);
          c[j] = _mm256_add_pd(c[j],tmp);
        }
      }

      int it = omp_get_thread_num();
      for(int j = 0; j < 4; j++)
      {
        s_[j][it] = s[j];
        c_[j][it] = c[j];
      }
    }


    // handcoded omp reduction
    for(int i = 1; i < nt; i++)
    {
      for(int j = 0; j < 4; j++)
      {
        __m256d sigma, oldS = s_[j][0];
        MM256_FAST2SUM(oldS, s_[j][i], s_[j][0], sigma);
        //MM256_2SUM(oldS, s_[0][i], s_[j][0], sigma);
        __m256d tmp = _mm256_add_pd(c_[j][i],sigma);
        c_[j][0] = _mm256_add_pd(c_[j][0], tmp);
      }
    }

    // return result, needs to be summed up
    for(int j = 0; j < 4; j++)
    {
      _mm256_storeu_pd(&res[j*4],  s_[j][0]);
      _mm256_storeu_pd(&resC[j*4], c_[j][0]);
    }
  }

}


// more accurate gemm product x'y AVX2 kernel
void dgemm_sc_prec_2_2(int nrows, const double *restrict x, const double *restrict y, double *restrict res, double *restrict resC)
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


#ifdef PHIST_HAVE_OPENMP
  int nt = omp_get_max_threads();
#else
  int nt = 1;
#endif

  {
    // buffer for omp thread result + padding to prevent false sharing
    __m256d s_[8][nt];
    __m256d c_[8][nt];

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

#pragma omp for schedule(static)
      for(int i = 0; i < nrows/2; i++)
      {
        __m256d xi = _mm256_load_pd(&x[4*i]);
        for(int j = 0; j < 2; j++)
        {
          __m256d yij = _mm256_set_pd(y[4*i+j+2],y[4*i+j+2],y[4*i+j],y[4*i+j]);
          __m256d p, pi;
          MM256_2MULTFMA(xi,yij,p,pi);
          __m256d sigma, oldS = s[j];
          MM256_FAST2SUM(oldS,p, s[j],sigma);
          //MM256_2SUM(oldS,p, s[j],sigma); // more FP ops than Kahan-style FAST2SUM, but exacter
          __m256d tmp = _mm256_add_pd(pi,sigma);
          c[j] = _mm256_add_pd(c[j],tmp);
        }
      }

      int it = omp_get_thread_num();
      for(int j = 0; j < 2; j++)
      {
        s_[j][it] = s[j];
        c_[j][it] = c[j];
      }
    }


    // handcoded omp reduction
    for(int i = 1; i < nt; i++)
    {
      for(int j = 0; j < 2; j++)
      {
        __m256d sigma, oldS = s_[j][0];
        MM256_FAST2SUM(oldS, s_[j][i], s_[j][0], sigma);
        //MM256_2SUM(oldS, s_[0][i], s_[j][0], sigma);
        __m256d tmp = _mm256_add_pd(c_[j][i],sigma);
        c_[j][0] = _mm256_add_pd(c_[j][0], tmp);
      }
    }


    // sum up 4 elements in mm256 to 2 elements in mm128
    // and return result, needs to be summed up again
    for(int j = 0; j < 2; j++)
    {
      __m128d s = _mm256_extractf128_pd(s_[j][0],0);
      __m128d p = _mm256_extractf128_pd(s_[j][0],1);
      __m128d c = _mm256_extractf128_pd(c_[j][0],0);
      __m128d pi = _mm256_extractf128_pd(c_[j][0],1);
      __m128d sigma, oldS = s;
      MM128_FAST2SUM(oldS,p,s,sigma);
      __m128d tmp = _mm_add_pd(pi,sigma);
      c = _mm_add_pd(c,tmp);
      _mm_storeu_pd(&res[j*2],  s);
      _mm_storeu_pd(&resC[j*2], c);
    }
  }

}

// more accurate gemm product x'y AVX2 kernel
void dgemm_sc_prec_2_1(int nrows, const double *restrict x, const double *restrict y, double *restrict res, double *restrict resC)
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


#ifdef PHIST_HAVE_OPENMP
  int nt = omp_get_max_threads();
#else
  int nt = 1;
#endif

  {
    // buffer for omp thread result + padding to prevent false sharing
    __m256d s_[8][nt];
    __m256d c_[8][nt];

#pragma omp parallel shared(s_,c_)
    {
      // initialize sum
      __m256d s = _mm256_setzero_pd();
      __m256d c = _mm256_setzero_pd();

#pragma omp for schedule(static)
      for(int i = 0; i < nrows/2; i++)
      {
        __m256d xi = _mm256_load_pd(&x[4*i]);
        __m256d yij = _mm256_set_pd(y[2*i+1],y[2*i+1],y[2*i],y[2*i]);
        __m256d p, pi;
        MM256_2MULTFMA(xi,yij,p,pi);
        __m256d sigma, oldS = s;
        MM256_FAST2SUM(oldS,p, s,sigma);
        //MM256_2SUM(oldS,p, s,sigma); // more FP ops than Kahan-style FAST2SUM, but exacter
        __m256d tmp = _mm256_add_pd(pi,sigma);
        c = _mm256_add_pd(c,tmp);
      }

      int it = omp_get_thread_num();
      s_[0][it] = s;
      c_[0][it] = c;
    }


    // handcoded omp reduction
    for(int i = 1; i < nt; i++)
    {
      __m256d sigma, oldS = s_[0][0];
      MM256_FAST2SUM(oldS, s_[0][i], s_[0][0], sigma);
      //MM256_2SUM(oldS, s_[0][i], s_[0][0], sigma);
      __m256d tmp = _mm256_add_pd(c_[0][i],sigma);
      c_[0][0] = _mm256_add_pd(c_[0][0], tmp);
    }


    // sum up 4 elements in mm256 to 2 elements in mm128
    // and return result, needs to be summed up again
    __m128d s = _mm256_extractf128_pd(s_[0][0],0);
    __m128d p = _mm256_extractf128_pd(s_[0][0],1);
    __m128d c = _mm256_extractf128_pd(c_[0][0],0);
    __m128d pi = _mm256_extractf128_pd(c_[0][0],1);
    __m128d sigma, oldS = s;
    MM128_FAST2SUM(oldS,p,s,sigma);
    __m128d tmp = _mm_add_pd(pi,sigma);
    c = _mm_add_pd(c,tmp);
    _mm_storeu_pd(res,  s);
    _mm_storeu_pd(resC, c);
  }
}


// more accurate gemm product x'y AVX2 kernel
void dgemm_sc_prec_4_2(int nrows, const double *restrict x, const double *restrict y, double *restrict res, double *restrict resC)
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


#ifdef PHIST_HAVE_OPENMP
  int nt = omp_get_max_threads();
#else
  int nt = 1;
#endif

  {
    // buffer for omp thread result + padding to prevent false sharing
    __m256d s_[8][nt];
    __m256d c_[8][nt];

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

#pragma omp for schedule(static)
      for(int i = 0; i < nrows; i++)
      {
        __m256d xi = _mm256_load_pd(&x[4*i]);
        for(int j = 0; j < 2; j++)
        {
          __m256d yij = _mm256_set_pd(y[2*i+j],y[2*i+j],y[2*i+j],y[2*i+j]);
          __m256d p, pi;
          MM256_2MULTFMA(xi,yij,p,pi);
          __m256d sigma, oldS = s[j];
          MM256_FAST2SUM(oldS,p, s[j],sigma);
          //MM256_2SUM(oldS,p, s[j],sigma); // more FP ops than Kahan-style FAST2SUM, but exacter
          __m256d tmp = _mm256_add_pd(pi,sigma);
          c[j] = _mm256_add_pd(c[j],tmp);
        }
      }

      int it = omp_get_thread_num();
      for(int j = 0; j < 2; j++)
      {
        s_[j][it] = s[j];
        c_[j][it] = c[j];
      }
    }


    // handcoded omp reduction
    for(int i = 1; i < nt; i++)
    {
      for(int j = 0; j < 2; j++)
      {
        __m256d sigma, oldS = s_[j][0];
        MM256_FAST2SUM(oldS, s_[j][i], s_[j][0], sigma);
        //MM256_2SUM(oldS, s_[0][i], s_[j][0], sigma);
        __m256d tmp = _mm256_add_pd(c_[j][i],sigma);
        c_[j][0] = _mm256_add_pd(c_[j][0], tmp);
      }
    }


    // return result, needs to be summed up again
    for(int j = 0; j < 2; j++)
    {
      _mm256_storeu_pd(&res[j*4],  s_[j][0]);
      _mm256_storeu_pd(&resC[j*4], c_[j][0]);
    }
  }

}

// more accurate gemm product x'y AVX2 kernel
void dgemm_sc_prec_4_1(int nrows, const double *restrict x, const double *restrict y, double *restrict res, double *restrict resC)
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


#ifdef PHIST_HAVE_OPENMP
  int nt = omp_get_max_threads();
#else
  int nt = 1;
#endif

  {
    // buffer for omp thread result + padding to prevent false sharing
    __m256d s_[8][nt];
    __m256d c_[8][nt];

#pragma omp parallel shared(s_,c_)
    {
      // initialize sum
      __m256d s = _mm256_setzero_pd();
      __m256d c = _mm256_setzero_pd();

#pragma omp for schedule(static)
      for(int i = 0; i < nrows; i++)
      {
        __m256d xi = _mm256_load_pd(&x[4*i]);
        __m256d yij = _mm256_broadcast_sd(&y[i]);
        __m256d p, pi;
        MM256_2MULTFMA(xi,yij,p,pi);
        __m256d sigma, oldS = s;
        MM256_FAST2SUM(oldS,p, s,sigma);
        //MM256_2SUM(oldS,p, s,sigma); // more FP ops than Kahan-style FAST2SUM, but exacter
        __m256d tmp = _mm256_add_pd(pi,sigma);
        c = _mm256_add_pd(c,tmp);
      }

      int it = omp_get_thread_num();
      s_[0][it] = s;
      c_[0][it] = c;
    }


    // handcoded omp reduction
    for(int i = 1; i < nt; i++)
    {
      __m256d sigma, oldS = s_[0][0];
      MM256_FAST2SUM(oldS, s_[0][i], s_[0][0], sigma);
      //MM256_2SUM(oldS, s_[0][i], s_[0][0], sigma);
      __m256d tmp = _mm256_add_pd(c_[0][i],sigma);
      c_[0][0] = _mm256_add_pd(c_[0][0], tmp);
    }


    // return result, needs to be summed up again
    _mm256_storeu_pd(res,  s_[0][0]);
    _mm256_storeu_pd(resC, c_[0][0]);
  }

}

// more accurate gemm product x'y AVX2 kernel
void dgemm_sc_prec_4_k(int nrows, int k, const double *restrict x, const double *restrict y, double *restrict res, double *restrict resC)
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
#ifdef PHIST_HAVE_OPENMP
  int nt = omp_get_max_threads();
#else
  int nt = 1;
#endif

  {
    // buffer for omp thread result + padding to prevent false sharing
    __m256d s_[(k/8+1)*8][nt];
    __m256d c_[(k/8+1)*8][nt];

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
          __m256d p, pi;
          MM256_2MULTFMA(xi,yij,p,pi);
          __m256d sigma, oldS = s[j];
          MM256_FAST2SUM(oldS,p, s[j],sigma);
          //MM256_2SUM(oldS,p, s[j],sigma); // more FP ops than Kahan-style FAST2SUM, but exacter
          __m256d tmp = _mm256_add_pd(pi,sigma);
          c[j] = _mm256_add_pd(c[j],tmp);
        }
      }

      int it = omp_get_thread_num();
      for(int j = 0; j < k; j++)
      {
        s_[j][it] = s[j];
        c_[j][it] = c[j];
      }
    }


    // handcoded omp reduction
    for(int i = 1; i < nt; i++)
    {
      for(int j = 0; j < k; j++)
      {
        __m256d sigma, oldS = s_[j][0];
        MM256_FAST2SUM(oldS, s_[j][i], s_[j][0], sigma);
        //MM256_2SUM(oldS, s_[0][i], s_[j][0], sigma);
        __m256d tmp = _mm256_add_pd(c_[j][i],sigma);
        c_[j][0] = _mm256_add_pd(c_[j][0], tmp);
      }
    }

    // return result, needs to be summed up
    for(int j = 0; j < k; j++)
    {
      _mm256_storeu_pd(&res[j*4],  s_[j][0]);
      _mm256_storeu_pd(&resC[j*4], c_[j][0]);
    }
  }

}

// more accurate gemm product x'y AVX2 kernel
void dgemm_sc_prec_2_k(int nrows, int k, const double *restrict x, const double *restrict y, double *restrict res, double *restrict resC)
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


#ifdef PHIST_HAVE_OPENMP
  int nt = omp_get_max_threads();
#else
  int nt = 1;
#endif

  {
    // buffer for omp thread result + padding to prevent false sharing
    __m256d s_[(k/8+1)*8][nt];
    __m256d c_[(k/8+1)*8][nt];

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
          __m256d yij = _mm256_set_pd(y[k*(i+1)+j],y[k*(i+1)+j],y[k*i+j],y[k*i+j]);
          __m256d p, pi;
          MM256_2MULTFMA(xi,yij,p,pi);
          __m256d sigma, oldS = s[j];
          MM256_FAST2SUM(oldS,p, s[j],sigma);
          //MM256_2SUM(oldS,p, s[j],sigma); // more FP ops than Kahan-style FAST2SUM, but exacter
          __m256d tmp = _mm256_add_pd(pi,sigma);
          c[j] = _mm256_add_pd(c[j],tmp);
        }
      }

      int it = omp_get_thread_num();
      for(int j = 0; j < k; j++)
      {
        s_[j][it] = s[j];
        c_[j][it] = c[j];
      }
    }


    // handcoded omp reduction
    for(int i = 1; i < nt; i++)
    {
      for(int j = 0; j < k; j++)
      {
        __m256d sigma, oldS = s_[j][0];
        MM256_FAST2SUM(oldS, s_[j][i], s_[j][0], sigma);
        //MM256_2SUM(oldS, s_[0][i], s_[j][0], sigma);
        __m256d tmp = _mm256_add_pd(c_[j][i],sigma);
        c_[j][0] = _mm256_add_pd(c_[j][0], tmp);
      }
    }


    // sum up 4 elements in mm256 to 2 elements in mm128
    // and return result, needs to be summed up again
    for(int j = 0; j < k; j++)
    {
      __m128d s = _mm256_extractf128_pd(s_[j][0],0);
      __m128d p = _mm256_extractf128_pd(s_[j][0],1);
      __m128d c = _mm256_extractf128_pd(c_[j][0],0);
      __m128d pi = _mm256_extractf128_pd(c_[j][0],1);
      __m128d sigma, oldS = s;
      MM128_FAST2SUM(oldS,p,s,sigma);
      __m128d tmp = _mm_add_pd(pi,sigma);
      c = _mm_add_pd(c,tmp);
      _mm_storeu_pd(&res[j*2],  s);
      _mm_storeu_pd(&resC[j*2], c);
    }
  }

}


// more accurate gemm product x'y AVX2 kernel
void dgemm_sc_prec_1_k(int nrows, int k, const double *restrict x, const double *restrict y, double *restrict res, double *restrict resC)
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


#ifdef PHIST_HAVE_OPENMP
  int nt = omp_get_max_threads();
#else
  int nt = 1;
#endif

  {
    // buffer for omp thread result + padding to prevent false sharing
    __m256d s_[(k/8+1)*8][nt];
    __m256d c_[(k/8+1)*8][nt];

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
          __m256d p, pi;
          MM256_2MULTFMA(xi,yij,p,pi);
          __m256d sigma, oldS = s[j];
          MM256_FAST2SUM(oldS,p, s[j],sigma);
          //MM256_2SUM(oldS,p, s[j],sigma); // more FP ops than Kahan-style FAST2SUM, but exacter
          __m256d tmp = _mm256_add_pd(pi,sigma);
          c[j] = _mm256_add_pd(c[j],tmp);
        }
      }

      int it = omp_get_thread_num();
      for(int j = 0; j < k; j++)
      {
        s_[j][it] = s[j];
        c_[j][it] = c[j];
      }
    }


    // handcoded omp reduction
    for(int i = 1; i < nt; i++)
    {
      for(int j = 0; j < k; j++)
      {
        __m256d sigma, oldS = s_[j][0];
        MM256_FAST2SUM(oldS, s_[j][i], s_[j][0], sigma);
        //MM256_2SUM(oldS, s_[0][i], s_[j][0], sigma);
        __m256d tmp = _mm256_add_pd(c_[j][i],sigma);
        c_[j][0] = _mm256_add_pd(c_[j][0], tmp);
      }
    }


    // sum up 4 elements in mm256 to 1 double
    // and store result in res,resC
    for(int j = 0; j < k; j++)
    {
      double sj[4], cj[4];
      _mm256_storeu_pd(sj, s_[j][0]);
      _mm256_storeu_pd(cj, c_[j][0]);
      prec_reduction_1(4, sj, cj, &res[j], &resC[j]);
printf("bla %e %e %e %e res %e\n", sj[0], sj[1], sj[2], sj[3], res[j]);
    }
  }

}


// precise calculation of C+Cc <- alpha*(A+Ac)*(B+Bc) + beta*(C+Cc)
void dgemm_prec(int m, int n, int k, double alpha, const double *restrict a, const double *restrict aC,
                                                   const double *restrict b, const double *restrict bC,
                                     double beta,        double *restrict c,       double *restrict cC)
{
#if defined(TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif

  for(int j = 0; j < n; j++)
  {
    for(int i = 0; i < m; i++)
    {
      // c_ <- beta*b[j*m+i]
      double c_, cC_;
      DOUBLE_2MULTFMA(beta,c[j*m+i], c_,cC_);
      cC_ = beta*cC[j*m+i]+cC_;

      for(int l = 0; l < k; l++)
      {
        // a_ <- alpha*a[l*m+i]
        double a_, aC_;
        DOUBLE_2MULTFMA(alpha,a[l*m+i], a_,aC_);
        aC_ = alpha*aC[l*m+i]+aC_;

        // tmp <- a_*b[j*k+l]
        double tmp, tmpC;
        DOUBLE_2MULTFMA(a_,b[j*k+l], tmp,tmpC);
        tmpC = a_*bC[j*k+l]+aC_*b[j*k+l]+tmpC;

        // c_ <- c_ + tmp
        double oldC = c_, oldCC_ = cC_;
        DOUBLE_2SUM(oldC,tmp, c_, cC_);
        cC_ = oldCC_ + cC_ + tmpC;
      }

      // round result
      DOUBLE_FAST2SUM(c_,cC_,c[j*m+i],cC[j*m+i]);
    }
  }
}


// calculates a possibly low rank approximation of a lower cholesky factor of an spd matrix
// higher-precision + pivoting + stable low rank approximation
void cholesky_prec(int n, double *restrict a, double *restrict aC, int *rank)
{
  // permutation
  int p[n];
  for(int i = 0; i < n; i++)
    p[i] = i;
  // constructed L
  double l[n*n], lC[n*n];
  for(int i = 0; i < n*n; i++)
  {
    l[i] = 0.;
    lC[i] = 0.;
  }
  // diagonal entries
  double d[n], dC[n];
  for(int i = 0; i < n; i++)
  {
    DOUBLE_FAST2SUM(a[i*n+i],aC[i*n+i],d[i],dC[i]);
  }

  *rank = 0;
  while(*rank < n)
  {
    // check rank
    double err = 0;
    for(int i = *rank; i < n; i++)
      err = err + d[p[i]];
//printf("step %d, err %e\n", *rank, err);
    if( err < 1.e-12 )
      break;

    int m = *rank;
    *rank = *rank + 1;
    // find next pivot
    {
      int i = m;
      for(int j = m+1; j < n; j++)
        if( d[p[j]] > d[p[i]] )
          i = j;
      // swap p[i] p[m]
      int tmp = p[i];
      p[i] = p[m];
      p[m] = tmp;
//printf("pivot %d, perm", i);
//for(int j = 0; j < n; j++)
//  printf(" %d",p[j]);
    }
//printf("\n");

    // l_m,p[m] = sqrt(d_p[m])
    // uses sqrt(x+y) = sqrt(x) + y/(2*sqrt(x)) + ...
    {
      double s,t,t_;
      DOUBLE_2SQRTFMA(d[p[m]],s,t);
      double tmp = dC[p[m]]/(2*s);
      DOUBLE_FAST2SUM(s,tmp,l[p[m]*n+m],t_);
      lC[p[m]*n+m] = t+t_;
//printf("l[m=%d,p[m=%d]] <- %e\n", m,m,l[p[m]*n+m]);
    }

    for(int i = m+1; i < n; i++)
    {
      // l_m,p[i] = 1/l_m,p[m] * ( a_p[m],p[i] - sum_j=0^m-2 l_j,p[m]*l_j,p[i] )
      {
        double s = -a[p[i]*n+p[m]];
        double t = -aC[p[i]*n+p[m]];
        for(int j = 0; j < m; j++)
        {
          double ljm = l[p[m]*n+j];
          double lCjm = lC[p[m]*n+j];
          double lji = l[p[i]*n+j];
          double lCji = lC[p[i]*n+j];
          double lj, ljC;
          DOUBLE_2MULTFMA(ljm,lji,lj,ljC);
//printf("subtract ljm*lji %e\n", lj);
          ljC = ljC+ljm*lCji+lji*lCjm;
          double oldS = s, t_;
          DOUBLE_2SUM(lj,oldS,s,t_);
          t = t + t_ + ljC;
        }
        double s_, t_;
        DOUBLE_FAST2SUM(s,t,s_,t_);
        // use a/(x+y) = a/x - a*y/x^2 + ... to calculate s_/l_m,p[m]
        DOUBLE_2DIVFMA(s_,l[p[m]*n+m],s,t);
        t = t - s*lC[p[m]*n+m]/l[p[m]*n+m]+t_/l[p[m]*n+m];
        DOUBLE_FAST2SUM(s,t,s_,t_);
        l[p[i]*n+m] = -s_;
        lC[p[i]*n+m] = -t_;
//printf("l[m=%d,p[i=%d]] <- %e\n", m,i,-s_);
      }
      // d_p[i] = d_p[i]-l_m,p[m]*l_m,p[i]
      {
        double s = -d[p[i]];
        double t = -dC[p[i]];
        double lmm = l[p[m]*n+m];
        double lCmm = lC[p[m]*n+m];
        double lmi = l[p[i]*n+m];
        double lCmi = lC[p[i]*n+m];
        double lm, lmC;
        DOUBLE_2MULTFMA(lmi,lmi,lm,lmC);
//printf("lm %e\n", lm);
        lmC = lmC+lmi*lCmi+lmi*lCmi;
        double oldS = s, t_;
        DOUBLE_2SUM(oldS,lm,s,t_);
        t = t + t_ + lmC;
        d[p[i]] = -s;
        dC[p[i]] = -t;
//printf("d[p[i=%d]] <- %e\n", i,-s);
      }
    }
  }

  // store result in a
  for(int i = 0; i < n; i++)
  {
    for(int j = 0; j < n; j++)
    {
      DOUBLE_FAST2SUM(l[i*n+j],lC[i*n+j],a[i*n+j],aC[i*n+j]);
    }
  }
}

