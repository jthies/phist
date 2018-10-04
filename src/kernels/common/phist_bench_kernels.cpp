/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
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

#include "phist_bench_kernels.h"
#include "stream_bench.h"

#define NUM_RUNS 40

extern "C" void phist_bench_stream_load(double* mean_bw, double* max_bw, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  double *data = NULL;
  *iflag=0;
  *max_bw = 0.;
  *mean_bw = 0.;
  PHIST_SOUT(PHIST_VERBOSE, "Streaming LOAD benchmark: ");
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
    PHIST_CHK_IERR(dbench_stream_load_create(&data,iflag),*iflag);
    for(int i = 0; i < NUM_RUNS; i++)
    {
      double bw = 0.;
      double res;
      PHIST_CHK_IERR(dbench_stream_load_run(data,&res,&bw,iflag),*iflag);
      *mean_bw+=bw;
      if( bw > *max_bw ) *max_bw = bw;
    }
    *mean_bw/=NUM_RUNS;
    PHIST_CHK_IERR(dbench_stream_load_destroy(data,iflag),*iflag);
    PHIST_SOUT(PHIST_VERBOSE, "measured %8.4g Gb/s (max) and %8.4g Gb/s (mean)\n", *max_bw/1.e9,*mean_bw/1.e9);
  }
#if defined(PHIST_BENCH_MASTER)&&defined(PHIST_HAVE_MPI)
  MPI_Bcast(mean_bw,1,MPI_DOUBLE,0,comm);
  MPI_Bcast(max_bw,1,MPI_DOUBLE,0,comm);
#endif
}

extern "C" void phist_bench_stream_store(double *mean_bw, double* max_bw, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  double *data = NULL;
  PHIST_SOUT(PHIST_VERBOSE, "Streaming STORE benchmark: ");
  *iflag=0;
  *max_bw = 0.;
  *mean_bw= 0.;
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
    PHIST_CHK_IERR(dbench_stream_store_create(&data,iflag),*iflag);
    for(int i = 0; i < NUM_RUNS; i++)
    {
      double bw = 0.;
      double res = 77.;
      PHIST_CHK_IERR(dbench_stream_store_run(data,&res,&bw,iflag),*iflag);
      *mean_bw+=bw;
      if( bw > *max_bw ) *max_bw = bw;
    }
    *mean_bw/=NUM_RUNS;
    PHIST_CHK_IERR(dbench_stream_store_destroy(data,iflag),*iflag);
    PHIST_SOUT(PHIST_VERBOSE, "measured %8.4g Gb/s (max) and %8.4g Gb/s (mean)\n", *max_bw/1.e9,*mean_bw/1.e9);
  }
#if defined(PHIST_BENCH_MASTER)&&defined(PHIST_HAVE_MPI)
  MPI_Bcast(mean_bw,1,MPI_DOUBLE,0,comm);
  MPI_Bcast(max_bw,1,MPI_DOUBLE,0,comm);
#endif
}

extern "C" void phist_bench_stream_triad(double* mean_bw, double* max_bw, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  double *x = NULL;
  double *y = NULL;
  double *z = NULL;
  PHIST_SOUT(PHIST_VERBOSE, "Streaming TRIAD benchmark: ");
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
    PHIST_CHK_IERR(dbench_stream_triad_create(&x,&y,&z,iflag),*iflag);
    for(int i = 0; i < NUM_RUNS; i++)
    {
      double bw = 0.;
      double res = -53.;
      PHIST_CHK_IERR(dbench_stream_triad_run(x,y,z,&res,&bw,iflag),*iflag);
      *mean_bw+=bw;
      if( bw > *max_bw ) *max_bw = bw;
    }
    *mean_bw/=NUM_RUNS;
    PHIST_CHK_IERR(dbench_stream_triad_destroy(x,y,z,iflag),*iflag);
    PHIST_SOUT(PHIST_VERBOSE, "measured %8.4g Gb/s (max) and %8.4g Gb/s (mean)\n", *max_bw/1.e9,*mean_bw/1.e9);
  }
#if defined(PHIST_BENCH_MASTER)&&defined(PHIST_HAVE_MPI)
  MPI_Bcast(mean_bw,1,MPI_DOUBLE,0,comm);
  MPI_Bcast(max_bw,1,MPI_DOUBLE,0,comm);
#endif
}


