/*! \file dot_kernels_prec.c
 * Fast parallel BLAS-axpy like functions with high precision for different blocksizes for mvec_module
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

// define our own datatype for aligned doubles
typedef double aligned_double __attribute__((aligned(64)));

// more accurate dot product x'x AVX2 kernel
void ddot_self_prec_4(int nrows, const aligned_double *restrict x, aligned_double *restrict res, double *restrict resC)
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
    __m256d s_[nt][8];
    __m256d c_[nt][8];

#pragma omp parallel shared(s_,c_)
    {
      // initialize sum
      __m256d s = _mm256_setzero_pd();
      __m256d c = _mm256_setzero_pd();

#pragma omp for schedule(static)
      for(int i = 0; i < nrows; i++)
      {
        __m256d xi = _mm256_load_pd(&x[4*i]);
        __m256d p, pi;
        MM256_2MULTFMA(xi,xi,p,pi);
        __m256d sigma, oldS = s;
        MM256_FAST2SUM(oldS,p, s,sigma);
        __m256d tmp = _mm256_add_pd(pi,sigma);
        c = _mm256_add_pd(c,tmp);
      }


      int it = omp_get_thread_num();
      s_[it][0] = s;
      c_[it][0] = c;
    }


    // handcoded omp reduction
    __m256d s = _mm256_setzero_pd();
    __m256d t = _mm256_setzero_pd();
    for(int i = 0; i < nt; i++)
    {
      __m256d oldS = s, oldT = t;
      MM256_4SUM(oldS,oldT,s_[i][0],c_[i][0],s,t);
    }

    // return result, needs to be summed up
    _mm256_storeu_pd(res,  s);
    _mm256_storeu_pd(resC, t);
  }
}


// more accurate dot product x'y AVX2 kernel
void ddot_prec_4(int nrows, const double *restrict x, const double *restrict y, double *restrict res, double *restrict resC)
{
#if defined(TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  if( !is_aligned(x,32) )
  {
    printf("%s: x not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)x);
    exit(1);
    return;
  }

  if( !is_aligned(y,32) )
  {
    printf("%s: y, not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
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
      __m256d s = _mm256_setzero_pd();
      __m256d c = _mm256_setzero_pd();

#pragma omp for schedule(static)
      for(int i = 0; i < nrows; i++)
      {
        __m256d xi = _mm256_load_pd(&x[4*i]);
        __m256d yi = _mm256_load_pd(&y[4*i]);
        __m256d p, pi;
        MM256_2MULTFMA(xi,yi,p,pi);
        __m256d sigma, oldS = s;
        MM256_FAST2SUM(oldS,p, s,sigma);
        __m256d tmp = _mm256_add_pd(pi,sigma);
        c = _mm256_add_pd(c,tmp);
      }

      int it = omp_get_thread_num();
      s_[it][0] = s;
      c_[it][0] = c;
    }


    // handcoded omp reduction
    __m256d s = _mm256_setzero_pd();
    __m256d t = _mm256_setzero_pd();
    for(int i = 0; i < nt; i++)
    {
      __m256d oldS = s, oldT = t;
      MM256_4SUM(oldS,oldT,s_[i][0],c_[i][0],s,t);
    }

    // return result, needs to be summed up
    _mm256_storeu_pd(res,  s);
    _mm256_storeu_pd(resC, t);
  }

}


// more accurate dot product x'x for block size 1
void ddot_self_prec_1(int nrows, const double *restrict x, double *restrict res, double *restrict resC)
{
#if defined(TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  double s[4], c[4];

  // assume appropriate padding with zeros!
  ddot_self_prec_4(nrows/4, x, s, c);

  // we still need to sum up s, c
  prec_reduction_1(4, s, c, res, resC);
}


// more accurate dot product x'y for block size 1
void ddot_prec_1(int nrows, const double *restrict x, const double *restrict y, double *restrict res, double *restrict resC)
{
#if defined(TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  double s[4], c[4];

  // assume appropriate padding with zeros!
  ddot_prec_4(nrows/4, x, y, s, c);

  // we still need to sum up s, c
  prec_reduction_1(4, s, c, res, resC);
}


// more accurate dot product x'x for block size 2
void ddot_self_prec_2(int nrows, const double *restrict x, double *restrict res, double *restrict resC)
{
#if defined(TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  double s[4], c[4];

  // assume appropriate padding with zeros!
  ddot_self_prec_4(nrows/2, x, s, c);

  // we still need to sum up s, c
  prec_reduction_2(2, s, c, res, resC);
}


// more accurate dot product x'y for block size 2
void ddot_prec_2(int nrows, const double *restrict x, const double *restrict y, double *restrict res, double *restrict resC)
{
#if defined(TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif
  double s[4], c[4];

  // assume appropriate padding with zeros!
  ddot_prec_4(nrows/2, x, y, s, c);

  // we still need to sum up s, c
  prec_reduction_2(2, s, c, res, resC);
}


