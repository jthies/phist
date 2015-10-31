#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "phist_kernels.h"
#include "phist_macros.h"

int main(int argc, char** argv)
{
  int iflag = 0;
  PHIST_ICHK_IERR(phist_kernels_init(&argc,&argv,&iflag),iflag);
PHIST_MAIN_TASK_BEGIN

#ifndef PHIST_KERNEL_LIB_MAGMA
  double bw;
  PHIST_ICHK_IERR(phist_bench_stream_load(&bw,&iflag),iflag);
  PHIST_ICHK_IERR(phist_bench_stream_store(&bw,&iflag),iflag);
  PHIST_ICHK_IERR(phist_bench_stream_triad(&bw,&iflag),iflag);
#endif

PHIST_MAIN_TASK_END
  PHIST_ICHK_IERR(phist_kernels_finalize(&iflag),iflag);
  }
