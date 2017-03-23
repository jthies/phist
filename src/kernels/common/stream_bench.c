/*! \file stream_bench.c
 * some simple benchmarking routines to obtain reference performance data
 * \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
 *
*/

#define _XOPEN_SOURCE 600
#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include <immintrin.h>
#include <omp.h>
#include <stdlib.h>

#include "stream_bench.h"

// unroll loops for SIMD usage
#define CHUNK 8

// define our own datatype for aligned doubles
typedef double aligned_double __attribute__((aligned(64)));


//! allocate memory for bench_stream_load_run
void dbench_stream_load_create(double** x, int* ierr)
{
  *ierr = posix_memalign((void**)x, 64, PHIST_BENCH_LARGE_N*sizeof(double));
  // init + NUMA touch
#pragma omp parallel for schedule(static)
  for(int i = 0; i < PHIST_BENCH_LARGE_N; i+=CHUNK)
    for(int j = i; j < i+CHUNK; j++)
      (*x)[j] = j*1./PHIST_BENCH_LARGE_N;
}


//! purely load dominated micro benchmark (for large n), actually calculates a sum and determines the bandwidth
//! \warning trust the compiler for now to do useful things
void dbench_stream_load_run(const aligned_double *restrict x, double *restrict res, double *restrict bw, int *restrict ierr)
{
  // start timing
  double wtime = omp_get_wtime();

#pragma omp parallel 
{
  double sum[CHUNK]={0.0};
#pragma omp for
  for(int i = 0; i < PHIST_BENCH_LARGE_N; i+=CHUNK)
  {
#pragma omp simd
    for (int j=0; j<CHUNK; j++)
    {
      sum[j] += x[i+j];
    }
  }
#pragma omp critical
  {
    // return result, so no smarty pants compiler may optimize away our calculation!
    *res = sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7];
  }
}
  // end timing
  wtime = omp_get_wtime() - wtime;

  // calculate bandwidth
  *bw = sizeof(double)*PHIST_BENCH_LARGE_N / wtime;

  *ierr = 0;
}


//! delete memory for bench_stream_load_run
void dbench_stream_load_destroy(double* x, int* ierr)
{
  free(x);
  *ierr = 0;
}


//! allocate memory for bench_stream_store_run
void dbench_stream_store_create(double** x, int* ierr)
{
  dbench_stream_load_create(x,ierr);
}


//! purely store dominated micro benchmark (for large n), actually writes sum values and determines the bandwidth
void dbench_stream_store_run(aligned_double *restrict x, const double* res, double *restrict bw, int *restrict ierr)
{
  // start timing
  double wtime = omp_get_wtime();

#ifdef PHIST_HAVE_AVX512
  __m512d v0 = _mm512_set_pd(*res,*res,*res,*res,*res,*res,*res,*res);
#elif defined(PHIST_HAVE_AVX)
  __m256d v4 = _mm256_set_pd(*res,*res,*res,*res);
  __m256d v0 = _mm256_set_pd(*res,*res,*res,*res);
#elif defined(PHIST_HAVE_SSE)
  __m128d v6 = _mm_set_pd(*res,*res);
  __m128d v4 = _mm_set_pd(*res,*res);
  __m128d v2 = _mm_set_pd(*res,*res);
  __m128d v0 = _mm_set_pd(*res,*res);
#endif
#pragma omp parallel for schedule(static)
  for(int i = 0; i < PHIST_BENCH_LARGE_N; i+=CHUNK)
  {
#ifdef PHIST_HAVE_AVX512
    _mm512_stream_pd(x+i,v0);
#elif defined(PHIST_HAVE_AVX)
    _mm256_stream_pd(x+i+0,v0);
    _mm256_stream_pd(x+i+4,v4);
#elif defined(PHIST_HAVE_SSE)
    _mm_stream_pd(x+i+0,v0);
    _mm_stream_pd(x+i+2,v2);
    _mm_stream_pd(x+i+4,v4);
    _mm_stream_pd(x+i+6,v6);
#endif
  }

  // end timing
  wtime = omp_get_wtime() - wtime;

  // calculate bandwidth
  *bw = sizeof(double)*PHIST_BENCH_LARGE_N / wtime;

  *ierr = 0;
}


//! delete memory for bench_stream_store_run
void dbench_stream_store_destroy(double* x, int* ierr)
{
  dbench_stream_load_destroy(x,ierr);
}


//! allocate memory for bench_stream_triad_run
void dbench_stream_triad_create(double** x, double** y, double** z, int* ierr)
{
  dbench_stream_load_create(x,ierr); if(*ierr != 0 ) return;
  dbench_stream_load_create(y,ierr); if(*ierr != 0 ) return;
  dbench_stream_load_create(z,ierr); if(*ierr != 0 ) return;
}


//! stream triad micro benchmark, determines the bandwidth
void dbench_stream_triad_run(const aligned_double *restrict x, const aligned_double *restrict y, aligned_double *restrict z, const double *restrict res, double *restrict bw, int *restrict ierr)
{
  // start timing
  double wtime = omp_get_wtime();

#ifdef PHIST_HAVE_AVX512
  __m512d a = _mm512_set_pd(*res,*res,*res,*res,*res,*res,*res,*res);
#elif defined(PHIST_HAVE_AVX)
  __m256d a = _mm256_set_pd(*res,*res,*res,*res);
#elif defined(PHIST_HAVE_SSE)
  __m128d a = _mm_set_pd(*res,*res);
#endif
#pragma omp parallel for schedule(static)
  for(int i = 0; i < PHIST_BENCH_LARGE_N; i+=CHUNK)
  {
#ifdef PHIST_HAVE_AVX512
    __m512d x0 = _mm512_load_pd(x+i+0);
    __m512d y0 = _mm512_load_pd(y+i+0);
    __m512d ay0 = _mm512_mul_pd(a,y0);
    __m512d z0 = _mm512_add_pd(x0,ay0);
    _mm512_stream_pd(z+i+0,z0);
#elif defined(PHIST_HAVE_AVX)
    __m256d x0 = _mm256_load_pd(x+i+0);
    __m256d x4 = _mm256_load_pd(x+i+4);

    __m256d y0 = _mm256_load_pd(y+i+0);
    __m256d y4 = _mm256_load_pd(y+i+4);

    __m256d ay0 = _mm256_mul_pd(a,y0);
    __m256d ay4 = _mm256_mul_pd(a,y4);

    __m256d z0 = _mm256_add_pd(x0,ay0);
    __m256d z4 = _mm256_add_pd(x4,ay4);

    _mm256_stream_pd(z+i+0,z0);
    _mm256_stream_pd(z+i+4,z4);
#elif defined(PHIST_HAVE_SSE)
    __m128d x0 = _mm_load_pd(x+i+0);
    __m128d x2 = _mm_load_pd(x+i+2);
    __m128d x4 = _mm_load_pd(x+i+4);
    __m128d x6 = _mm_load_pd(x+i+6);

    __m128d y0 = _mm_load_pd(y+i+0);
    __m128d y2 = _mm_load_pd(y+i+2);
    __m128d y4 = _mm_load_pd(y+i+4);
    __m128d y6 = _mm_load_pd(y+i+6);

    __m128d ay0 = _mm_mul_pd(a,y0);
    __m128d ay2 = _mm_mul_pd(a,y2);
    __m128d ay4 = _mm_mul_pd(a,y4);
    __m128d ay6 = _mm_mul_pd(a,y6);

    __m128d z0 = _mm_add_pd(x0,ay0);
    __m128d z2 = _mm_add_pd(x2,ay2);
    __m128d z4 = _mm_add_pd(x4,ay4);
    __m128d z6 = _mm_add_pd(x6,ay6);

    _mm_stream_pd(z+i+0,z0);
    _mm_stream_pd(z+i+2,z2);
    _mm_stream_pd(z+i+4,z4);
    _mm_stream_pd(z+i+6,z6);
#endif
  }

  // end timing
  wtime = omp_get_wtime() - wtime;

  // calculate bandwidth
  *bw = 3*sizeof(double)*PHIST_BENCH_LARGE_N / wtime;

  *ierr = 0;
}


//! delete memory for bench_stream_triad_run
void dbench_stream_triad_destroy(double* x, double* y, double *z, int* ierr)
{
  dbench_stream_load_destroy(x,ierr); if(*ierr != 0 ) return;
  dbench_stream_load_destroy(y,ierr); if(*ierr != 0 ) return;
  dbench_stream_load_destroy(z,ierr); if(*ierr != 0 ) return;
}

