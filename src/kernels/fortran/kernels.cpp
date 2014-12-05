#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "phist_macros.h"
#ifdef PHIST_HAVE_TEUCHOS
#include "phist_trilinos_macros.h"
#endif
#include "../phist_kernels.h"

#include "phist_typedefs.h"
#include "typedefs.hpp"
#include "phist_ScalarTraits.hpp"

#ifdef PHIST_HAVE_TEUCHOS
#include "Teuchos_TimeMonitor.hpp"
#endif
#ifdef PHIST_HAVE_LIKWID
#include <likwid.h>
#endif

#include <malloc.h>
#include <cstring>
#ifdef PHIST_HAVE_OPENMP
#include <omp.h>
#else
namespace{
int omp_get_thread_num() {return 0;}
}
#endif
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
#ifdef PHIST_HAVE_OPENMP
  PHIST_SOUT(PHIST_WARNING,"Trying to pin the OpenMP threads to the cores without caring about the real topology of the system.\n"\
                           "Thus, this may fail/crash. Use GHOST kernels!\n"\
                           "Assuming %d NUMA domains per node, with %d cores per NUMA domain and %d SMT threads per core\n",
                           PHIST_FORTRAN_PIN_NUMA_DOMAINS, PHIST_FORTRAN_PIN_CORES_PER_NUMA, PHIST_FORTRAN_PIN_SMTTHREADS_PER_CORE);
#pragma omp parallel
  {
    // WARNING: the pinning may crash/fail depending on the settings/machine
    //          use GHOST kernels
    int tid = omp_get_thread_num();
    int coreid = tid + PHIST_FORTRAN_PIN_CORES_PER_NUMA*( rank % PHIST_FORTRAN_PIN_NUMA_DOMAINS);
    const int totalSmtCores = PHIST_FORTRAN_PIN_CORES_PER_NUMA*PHIST_FORTRAN_PIN_NUMA_DOMAINS*PHIST_FORTRAN_PIN_SMTTHREADS_PER_CORE;
    cpu_set_t *cpu = CPU_ALLOC(totalSmtCores);
    size_t size = CPU_ALLOC_SIZE(totalSmtCores);
    CPU_ZERO_S(size, cpu);
    CPU_SET_S(coreid, size, cpu);
    sched_setaffinity(0, size, cpu);
    CPU_FREE(cpu);
  }
#endif

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
PHIST_CXX_TIMER_SUMMARIZE
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
void phist_map_get_local_length(const_map_ptr_t vmap, lidx_t* nloc, int* ierr);
void phist_map_get_ilower(const_map_ptr_t vmap, gidx_t* ilower, int* ierr);
void phist_map_get_iupper(const_map_ptr_t vmap, gidx_t* iupper, int* ierr);

#ifdef PHIST_TIMEMONITOR
void phist_totalMatVecCount()
{
  ENTER_FCN(__FUNCTION__);
}
#endif

#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "../kernels_noimpl.c"
#include "../carp_noimpl.c"
#include "../kernels_nogpu.c"

#include "phist_gen_c.h"
#include "../kernels_noimpl.c"
#include "../carp_noimpl.c"
#include "../kernels_nogpu.c"
#endif

#include "phist_gen_z.h"
#include "../kernels_noimpl.c"
#include "../carp_noimpl.c"
#include "../kernels_nogpu.c"

} //extern "C"

#include "phist_gen_d.h"
#include "kernels_def.hpp"
#include "carp_def.hpp"


