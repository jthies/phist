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

#include "ghost/bench.h"

/* with GHOST we use the builtin stream benchmarks because they will run on any supported hardware, 
   in particular CUDA-capable GPUs. We should keep in mind that the GHOST benchmark is a bit different
   from the one in kernels/common: ghost measures the average over 40 runs, phist the maximum over 10.
 */

extern "C" void phist_bench_stream_load(double* max_bw, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_SOUT(PHIST_VERBOSE, "Streaming LOAD benchmark: ");
  PHIST_CHK_GERR(ghost_bench_stream(GHOST_BENCH_STREAM_LOAD, max_bw),*iflag);
  *max_bw*=1.0e9; // GHOST returns GB/s, the PHIST equivalent B/s
  PHIST_SOUT(PHIST_VERBOSE, "measured %8.4g Gb/s\n", *max_bw/1.e9);
}

extern "C" void phist_bench_stream_store(double* max_bw, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_SOUT(PHIST_VERBOSE, "Streaming STORE benchmark: ");
  PHIST_CHK_GERR(ghost_bench_stream(GHOST_BENCH_STREAM_STORE, max_bw),*iflag);
  *max_bw*=1.0e9; // GHOST returns GB/s, the PHIST equivalent B/s
  PHIST_SOUT(PHIST_VERBOSE, "measured %8.4g Gb/s\n", *max_bw/1.e9);
}

extern "C" void phist_bench_stream_triad(double* max_bw, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_SOUT(PHIST_VERBOSE, "Streaming TRIAD benchmark: ");
  PHIST_CHK_GERR(ghost_bench_stream(GHOST_BENCH_STREAM_TRIAD, max_bw),*iflag);
  *max_bw*=1.0e9; // GHOST returns GB/s, the PHIST equivalent B/s
  PHIST_SOUT(PHIST_VERBOSE, "measured %8.4g Gb/s\n", *max_bw/1.e9);
}


