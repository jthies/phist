/*! \file bench_kernels.c
 * some simple benchmarking routines to obtain reference performance data
 * \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
 *
*/

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include <emmintrin.h>
#include <omp.h>


// define our own datatype for aligned doubles
typedef double aligned_double __attribute__((aligned(64)));


//! allocate memory for bench_stream_load_run
void dbench_stream_load_create(double** x, int* ierr)
{
  *ierr = posix_memalign((void**)x, 64, PHIST_BENCH_LARGE_N*sizeof(double));
  // init + NUMA touch
#pragma omp parallel for schedule(static)
  for(int i = 0; i < PHIST_BENCH_LARGE_N; i+=4)
    for(int j = i; j < i+4; j++)
      (*x)[j] = j*1./PHIST_BENCH_LARGE_N;
}


//! purely load dominated micro benchmark (for large n), actually calculates a sum and determines the bandwidth
//! \warning trust the compiler for now to do useful things
void dbench_stream_load_run(const aligned_double *restrict x, double *restrict res, double *restrict bw, int *restrict ierr)
{
  // start timing
  double wtime = omp_get_wtime();

  double sum0 = 0, sum1 = 0, sum2 = 0, sum3 = 0.;
#pragma omp parallel for reduction(+:sum0,sum1,sum2,sum3) schedule(static)
  for(int i = 0; i < PHIST_BENCH_LARGE_N; i+=4)
  {
    sum0 += x[i+0];
    sum1 += x[i+1];
    sum2 += x[i+2];
    sum3 += x[i+3];
  }
  // return result, so no smarty pants compiler may optimize away our calculation!
  *res = sum0+sum1+sum2+sum3;

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

  __m128d v2 = _mm_set_pd(*res,*res);
  __m128d v0 = _mm_set_pd(*res,*res);
#pragma omp parallel for schedule(static)
  for(int i = 0; i < PHIST_BENCH_LARGE_N; i+=4)
  {
    _mm_stream_pd(x+i+0,v0);
    _mm_stream_pd(x+i+2,v2);
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

  __m128d a = _mm_set_pd(*res,*res);
#pragma omp parallel for schedule(static)
  for(int i = 0; i < PHIST_BENCH_LARGE_N; i+=4)
  {
    __m128d x0 = _mm_load_pd(x+i+0);
    __m128d x2 = _mm_load_pd(x+i+2);

    __m128d y0 = _mm_load_pd(y+i+0);
    __m128d y2 = _mm_load_pd(y+i+2);

    __m128d ay0 = _mm_mul_pd(a,y0);
    __m128d ay2 = _mm_mul_pd(a,y2);

    __m128d z0 = _mm_add_pd(x0,ay0);
    __m128d z2 = _mm_add_pd(x2,ay2);

    _mm_stream_pd(z+i+0,z0);
    _mm_stream_pd(z+i+2,z2);
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

