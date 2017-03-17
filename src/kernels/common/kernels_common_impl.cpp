#include "phist_config.h"

#include "phist_tools.h"
#include "phist_kernels.h"
#include "phist_perfcheck.hpp"
#include "phist_random.h"

#include <sys/types.h>
#include <unistd.h>

//! this function should not be called by the user but by each kernel lib in kernels_init()
extern "C" void phist_kernels_common_init(int *argc, char*** argv, int* iflag)
{
  *iflag=0;
  phist_random_init();
#ifdef PHIST_HAVE_LIKWID
PHIST_TASK_DECLARE(likwidInitTask)
PHIST_TASK_BEGIN(likwidInitTask)
  LIKWID_MARKER_INIT;
#pragma omp parallel
  {
    LIKWID_MARKER_THREADINIT;
    LIKWID_MARKER_START("phist");
  }
PHIST_TASK_END(iflag)
#endif
  // at this point, check for the environment variable PHIST_RUN_DEBUGGER and if it is set,
  // go into an infinite loop that can only be broken after attaching a debugger
  // this code was copied from the OpenMPI  FAQ: http://www.open-mpi.de/faq/?category=debugging
  if (getenv("PHIST_ATTACH_GDB")!=NULL)
  {
    int i = 0;
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    PHIST_SOUT(PHIST_ERROR,"You have set the environment variable 'PHIST_ATTACH_GDB'.\n");
    PHIST_ORDERED_OUT(PHIST_ERROR,MPI_COMM_WORLD,"PID %d on %s ready for attach.\n", getpid(), hostname);
    PHIST_SOUT(PHIST_ERROR,"To attach gdb to a process, run 'gdb <executable> <PID> on the respective node.\n");
    PHIST_SOUT(PHIST_ERROR,"Press any <ENTER> to continue...\n");
    int rank=0;
#ifdef PHIST_HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
    if (rank==0) fgetc(stdin);
  }
#ifdef PHIST_HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}
 


extern "C" void phist_kernels_common_finalize(int *iflag)
{
  *iflag=0;
#ifdef PHIST_HAVE_LIKWID
PHIST_TASK_DECLARE(LikwidFinalizeTask)
PHIST_TASK_BEGIN(LikwidFinalizeTask)
#pragma omp parallel
  {
    LIKWID_MARKER_STOP("phist");
  }
  LIKWID_MARKER_CLOSE;
PHIST_TASK_END(iflag)
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
  phist_gidx ilower;
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

