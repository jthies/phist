#include "phist_config.h"

#include <cstdio>
#include <cstdlib>

#include <iostream>
#include "phist_macros.h"

#ifdef PHIST_HAVE_LIKWID
#include <likwid.h>
#endif

#include "phist_bench_kernels.h"
#include "stream_bench.h"

extern "C" {
void phist_bench_stream_load(double* max_bw, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  double *data = NULL;
  PHIST_SOUT(PHIST_INFO, "Streaming LOAD benchmark: ");
  PHIST_CHK_IERR(dbench_stream_load_create(&data,iflag),*iflag);
  *max_bw = 0.;
  for(int i = 0; i < 10; i++)
  {
    double bw = 0.;
    double res;
    PHIST_CHK_IERR(dbench_stream_load_run(data,&res,&bw,iflag),*iflag);
    if( bw > *max_bw ) *max_bw = bw;
  }
  PHIST_CHK_IERR(dbench_stream_load_destroy(data,iflag),*iflag);
  PHIST_SOUT(PHIST_INFO, "measured %8.4g Gb/s\n", *max_bw/1.e9);
}

void phist_bench_stream_store(double* max_bw, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  double *data = NULL;
  PHIST_SOUT(PHIST_INFO, "Streaming STORE benchmark: ");
  PHIST_CHK_IERR(dbench_stream_store_create(&data,iflag),*iflag);
  *max_bw = 0.;
  for(int i = 0; i < 10; i++)
  {
    double bw = 0.;
    double res = 77.;
    PHIST_CHK_IERR(dbench_stream_store_run(data,&res,&bw,iflag),*iflag);
    if( bw > *max_bw ) *max_bw = bw;
  }
  PHIST_CHK_IERR(dbench_stream_store_destroy(data,iflag),*iflag);
  PHIST_SOUT(PHIST_INFO, "measured %8.4g Gb/s\n", *max_bw/1.e9);
}

void phist_bench_stream_triad(double* max_bw, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  double *x = NULL;
  double *y = NULL;
  double *z = NULL;
  PHIST_SOUT(PHIST_INFO, "Streaming TRIAD benchmark: ");
  PHIST_CHK_IERR(dbench_stream_triad_create(&x,&y,&z,iflag),*iflag);
  *max_bw = 0.;
  for(int i = 0; i < 10; i++)
  {
    double bw = 0.;
    double res = -53.;
    PHIST_CHK_IERR(dbench_stream_triad_run(x,y,z,&res,&bw,iflag),*iflag);
    if( bw > *max_bw ) *max_bw = bw;
  }
  PHIST_CHK_IERR(dbench_stream_triad_destroy(x,y,z,iflag),*iflag);
  PHIST_SOUT(PHIST_INFO, "measured %8.4g Gb/s\n", *max_bw/1.e9);
}
}


