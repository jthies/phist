#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "phist_macros.h"
#ifdef PHIST_HAVE_BELOS
#include "phist_trilinos_macros.h"
#endif
#include "../phist_kernels.h"

#include "phist_typedefs.h"
#include "typedefs.hpp"
#include "phist_ScalarTraits.hpp"


#ifdef PHIST_TIMEMONITOR
#include "phist_timemonitor.hpp"
namespace phist_TimeMonitor
{
  Timer::TimeDataMap Timer::_timingResults;
}
#endif

#ifdef PHIST_HAVE_LIKWID
#include <likwid.h>
#endif

#include <malloc.h>
#include <cstring>
#include <omp.h>
#include <sched.h>
#include <iostream>


extern "C" {

#ifndef PHIST_HAVE_MPI
#error "fortran kernels only work with MPI"
#endif


// comment in for glibc/gcc and memory alignment problems...
/*
// force all memory allocations to be 16 byte aligned!
static void*(*old_malloc_hook)(size_t,const void*);
static void * my_malloc_hook (size_t size, const void *caller)
{
   void *result;
   // Restore all old hooks
   __malloc_hook = old_malloc_hook;
   // Call recursively
   result = memalign (16, size);
   // Save underlying hooks
   old_malloc_hook = __malloc_hook;
   // printf might call malloc, so protect it too.
   printf ("malloc (%u) returns %p\n", (unsigned int) size, result);
   // Restore our own hooks
   __malloc_hook = my_malloc_hook;
   return result;
}
static void my_init_hook (void)
{
  old_malloc_hook = __malloc_hook;
  __malloc_hook = my_malloc_hook;
}
// Override initializing hook from the C library.
void (*__malloc_initialize_hook) (void) = my_init_hook;
*/


void init_random_seed(void);

// initialize
void phist_kernels_init(int* argc, char*** argv, int* ierr)
{
  *ierr=0;
  int initialized = 0;

  PHIST_CHK_IERR( *ierr = MPI_Initialized(&initialized), *ierr);
  if (!initialized) {
  	PHIST_CHK_IERR( *ierr = MPI_Init(argc, argv), *ierr);
  }

  int rank;
  PHIST_CHK_IERR( *ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank), *ierr);
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    int coreid = tid;// + 10*(rank%2);
    cpu_set_t *cpu = CPU_ALLOC(40);
    size_t size = CPU_ALLOC_SIZE(40);
    CPU_ZERO_S(size, cpu);
    CPU_SET_S(coreid, size, cpu);
    sched_setaffinity(0, size, cpu);
    CPU_FREE(cpu);
  }

  std::ostringstream oss;
  oss << "rank: " << rank << ", cores(thread):";
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    int coreid = sched_getcpu();
#pragma omp critical
    oss << "\t" << coreid << " (" << tid << ")";
  }
  std::cout << oss.str() << std::endl;

#ifdef PHIST_HAVE_LIKWID
  LIKWID_MARKER_INIT;
#pragma omp parallel
  {
    LIKWID_MARKER_THREADINIT;
    LIKWID_MARKER_START("phist<fortran>");
  }
#endif
  init_random_seed();
}

// finalize fortran
void phist_kernels_finalize(int* ierr)
{
#ifdef PHIST_HAVE_LIKWID
#pragma omp parallel
  {
    LIKWID_MARKER_STOP("phist<fortran>");
  }
  LIKWID_MARKER_CLOSE;
#endif
#ifdef PHIST_TIMEMONITOR
  phist_TimeMonitor::Timer::summarize();
#endif
  PHIST_CHK_IERR( *ierr = MPI_Finalize(), *ierr);
  *ierr=0;
}

void phist_comm_create(comm_ptr_t* vcomm, int* ierr);
void phist_comm_delete(comm_ptr_t vcomm, int* ierr);
void phist_comm_get_rank(const_comm_ptr_t vcomm, int* rank, int* ierr);
void phist_comm_get_size(const_comm_ptr_t vcomm, int* size, int* ierr);
void phist_map_create(map_ptr_t* vmap, const_comm_ptr_t vcomm, gidx_t nglob, int *ierr);
void phist_map_delete(map_ptr_t vmap, int *ierr);
void phist_map_get_comm(const_map_ptr_t vmap, const_comm_ptr_t* vcomm, int* ierr);
void phist_map_get_local_length(const_map_ptr_t vmap, int* nloc, int* ierr);
void phist_map_get_ilower(const_map_ptr_t vmap, gidx_t* ilower, int* ierr);
void phist_map_get_iupper(const_map_ptr_t vmap, gidx_t* iupper, int* ierr);

#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "../kernels_noimpl.c"
#include "../carp_noimpl.c"

#include "phist_gen_c.h"
#include "../kernels_noimpl.c"
#include "../carp_noimpl.c"
#endif

#include "phist_gen_z.h"
#include "../kernels_noimpl.c"
#include "../carp_noimpl.c"

} //extern "C"

#include "phist_gen_d.h"
#include "kernels_def.hpp"
#include "carp_def.hpp"


