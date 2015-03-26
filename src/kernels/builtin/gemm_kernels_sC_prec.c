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
  printf("Entering %s\n", __FUNCTION__);
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
  printf("Entering %s\n", __FUNCTION__);
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
  printf("Entering %s\n", __FUNCTION__);
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
  printf("Entering %s\n", __FUNCTION__);
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
  printf("Entering %s\n", __FUNCTION__);
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
  printf("Entering %s\n", __FUNCTION__);
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
  printf("Entering %s\n", __FUNCTION__);
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
  printf("Entering %s\n", __FUNCTION__);
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
  printf("Entering %s\n", __FUNCTION__);
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
  printf("Entering %s\n", __FUNCTION__);
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
