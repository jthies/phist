/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_config.h"

/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <cstdio>
#include <cstdlib>

#include <iostream>
#include "phist_macros.h"

#ifdef PHIST_HAVE_LIKWID
#include <likwid.h>
#endif

#include <ghost/bench.h>

/* with GHOST we use the builtin stream benchmarks because they will run on any supported hardware, 
   in particular CUDA-capable GPUs. We should keep in mind that the GHOST benchmark is a bit different
   from the one in kernels/common: ghost measures the average over 40 runs, phist the maximum over 10.
 */

extern "C" void phist_bench_stream_load(double* mean_bw, double* max_bw, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_TASK_DECLARE(BenchTask)
  PHIST_TASK_BEGIN(BenchTask)
  PHIST_SOUT(PHIST_VERBOSE, "Streaming LOAD benchmark: ");
#if PHIST_BENCH_LARGE_N<=0
  PHIST_SOUT(PHIST_WARNING, "benchmark skipped because PHIST_BENCH_LARGE_N<=0\n");
  *iflag=1;
  *mean_bw=0.0;
  *max_bw=0.0;
#else
  PHIST_CHK_GERR(ghost_bench_bw(GHOST_BENCH_LOAD, mean_bw, max_bw),*iflag);
  *max_bw*=1.0e9; // GHOST returns GB/s, the PHIST equivalent B/s
  *mean_bw*=1.0e9;
  PHIST_SOUT(PHIST_VERBOSE, "measured %8.4g Gb/s (max) and %8.4g Gb/s (mean)\n", *max_bw/1.e9,*mean_bw/1.e9);
#endif
  PHIST_TASK_END(iflag)
}

extern "C" void phist_bench_stream_store(double* mean_bw, double* max_bw, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_TASK_DECLARE(BenchTask)
  PHIST_TASK_BEGIN(BenchTask)
  PHIST_SOUT(PHIST_VERBOSE, "Streaming STORE benchmark: ");
#if PHIST_BENCH_LARGE_N<=0
  PHIST_SOUT(PHIST_WARNING, "benchmark skipped because PHIST_BENCH_LARGE_N<=0\n");
  *iflag=1;
  *mean_bw=0.0;
  *max_bw=0.0;
#else
  PHIST_CHK_GERR(ghost_bench_bw(GHOST_BENCH_STORE, mean_bw,max_bw),*iflag);
  *max_bw*=1.0e9; // GHOST returns GB/s, the PHIST equivalent B/s
  *mean_bw*=1.0e9;
  PHIST_SOUT(PHIST_VERBOSE, "measured %8.4g Gb/s (max) and %8.4g Gb/s (mean)\n", *max_bw/1.e9,*mean_bw/1.e9);
#endif
  PHIST_TASK_END(iflag)
}

extern "C" void phist_bench_stream_copy(double* mean_bw, double* max_bw, int* iflag)
{
  PHIST_CHK_IERR(*iflag=PHIST_NOT_IMPLEMENTED, *iflag);
#if 0
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  double *x = NULL;
  double *y = NULL;
  PHIST_SOUT(PHIST_VERBOSE, "Streaming COPY benchmark: ");
  *max_bw = 0.;
  *mean_bw= 0.;
  *iflag=0;
#if PHIST_BENCH_LARGE_N<=0
  PHIST_SOUT(PHIST_VERBOSE, "skipped because PHIST_BENCH_LARGE_N<=0\n");
  return;
#endif

int execute=0;
#if defined(PHIST_BENCH_MASTER)&&defined(PHIST_HAVE_MPI)
// execute stream benchmarks only on rank 0 and bcast
MPI_Comm comm = phist_get_default_comm();
MPI_Comm_rank(comm,&execute);
#endif
  if (execute==0)
  {
    PHIST_CHK_IERR(dbench_stream_copy_create(&x,&y,iflag),*iflag);
    for(int i = 0; i < NUM_RUNS; i++)
    {
      double bw = 0.;
      double res = -53.;
      PHIST_CHK_IERR(dbench_stream_copy_run(x,y,&res,&bw,iflag),*iflag);
      *mean_bw+=bw;
      if( bw > *max_bw ) *max_bw = bw;
    }
    *mean_bw/=NUM_RUNS;
    PHIST_CHK_IERR(dbench_stream_copy_destroy(x,y,iflag),*iflag);
    PHIST_SOUT(PHIST_VERBOSE, "measured %8.4g Gb/s (max) and %8.4g Gb/s (mean)\n", *max_bw/1.e9,*mean_bw/1.e9);
  }
#if defined(PHIST_BENCH_MASTER)&&defined(PHIST_HAVE_MPI)
  MPI_Bcast(mean_bw,1,MPI_DOUBLE,0,comm);
  MPI_Bcast(max_bw,1,MPI_DOUBLE,0,comm);
#endif
#endif
}

extern "C" void phist_bench_stream_triad(double* mean_bw, double* max_bw, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_TASK_DECLARE(BenchTask)
  PHIST_TASK_BEGIN(BenchTask)
  PHIST_SOUT(PHIST_VERBOSE, "Streaming TRIAD benchmark: ");
#if PHIST_BENCH_LARGE_N<=0
  PHIST_SOUT(PHIST_WARNING, "benchmark skipped because PHIST_BENCH_LARGE_N<=0\n");
  *iflag=1;
  *mean_bw=0.0;
  *max_bw=0.0;
#else
  PHIST_CHK_GERR(ghost_bench_bw(GHOST_BENCH_STREAM_TRIAD, mean_bw, max_bw),*iflag);
  *max_bw*=1.0e9; // GHOST returns GB/s, the PHIST equivalent B/s
  *mean_bw*=1.0e9;
  PHIST_SOUT(PHIST_VERBOSE, "measured %8.4g Gb/s (max) and %8.4g Gb/s (mean)\n", *max_bw/1.e9,*mean_bw/1.e9);
#endif
  PHIST_TASK_END(iflag)
}



