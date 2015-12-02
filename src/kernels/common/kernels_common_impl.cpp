#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#ifdef PHIST_HAVE_LIKWID
#include <likwid.h>
#endif
#include <stdlib.h>

#include "phist_macros.h"
#include "phist_tasks.h"
#include "phist_fcntrace.hpp"
#include "phist_perfcheck.hpp"
#include "phist_kernels.h"
#include "phist_random.h"

//! this function should not be called by the user but by each kernel lib in kernels_init()
extern "C" void phist_kernels_common_init(int *argc, char*** argv, int* iflag)
{
  phist_random_init();
#ifdef PHIST_HAVE_LIKWID
GHOST_TASK_DECLARE(likwidInitTask)
GHOST_TASK_BEGIN(likwidInitTask)
  LIKWID_MARKER_INIT;
#pragma omp parallel
  {
    LIKWID_MARKER_THREADINIT;
    LIKWID_MARKER_START("phist");
  }
GHOST_TASK_END(likwidInitTask)
#endif
}

extern "C" void phist_kernels_common_finalize(int *iflag)
{
  *iflag=0;
#ifdef PHIST_HAVE_LIKWID
GHOST_TASK_DECLARE(LikwidFinalizeTask)
GHOST_TASK_BEGIN(LikwidFinalizeTask)
#pragma omp parallel
  {
    LIKWID_MARKER_STOP("phist");
  }
  LIKWID_MARKER_CLOSE;
GHOST_TASK_END(LikwidFinalizeTask)
#endif
#if defined(PHIST_TIMEMONITOR) || defined(PHIST_TIMEMONITOR_PERLINE)
PHIST_CXX_TIMER_SUMMARIZE;
#endif

PHIST_PERFCHECK_SUMMARIZE(PHIST_INFO);

#ifdef PHIST_PERFCHECK
  // prevent some strange memory errors during deallocation (due to shared lib context?)
  phist_PerfCheck::benchmarks.clear();
#endif
  *iflag=0;
}

typedef struct {
  int lda;
  int lnrows;
  int lncols;
  gidx_t ilower;
  double* data;
} dwrap;


#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "kernels_common_impl_def.h"

#include "phist_gen_c.h"
#include "kernels_common_impl_def.h"
#endif

#include "phist_gen_d.h"
#include "kernels_common_impl_def.h"

#include "phist_gen_z.h"
#include "kernels_common_impl_def.h"

