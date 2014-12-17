#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#else
#error "builtin kernels only work with MPI"
#endif

#include "phist_macros.h"
#ifdef PHIST_HAVE_TEUCHOS
#include "phist_trilinos_macros.h"
#endif
#include "../phist_kernels.h"

#include "phist_typedefs.h"
#include "typedefs.hpp"
#include "phist_ScalarTraits.hpp"

#ifdef PHIST_HAVE_LIKWID
#include <likwid.h>
#endif

#include <cstring>
#include <sys/resource.h>

#ifdef PHIST_HAVE_OPENMP
#include <omp.h>
#else
namespace{
int omp_get_thread_num() {return 0;}
int omp_get_num_threads() {return 1;}
}
#endif
#ifdef PHIST_KERNEL_LIB_BUILTIN_PIN_THREADS
#include <string>
#include <map>
#include <sstream>
#include <sched.h>

namespace {
// thread pinning using simple heuristics
void pinThreads()
{
  int ierr = 0;
  // gather system information
  // threads
  int nThreads = 1;
#pragma omp parallel shared(nThreads)
  {
#pragma omp critical
    nThreads = omp_get_num_threads();
  }
  // MPI rank
  int myRank, nRanks;
  PHIST_CHK_IERR( ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myRank), ierr);
  PHIST_CHK_IERR( ierr = MPI_Comm_size(MPI_COMM_WORLD, &nRanks), ierr);
  // nodes
  int ranksPerNode = 1;
  int myRankInNode = 0;
  {
    char myNodeId[MPI_MAX_PROCESSOR_NAME];
    int nodeId_strlen = 0;
    PHIST_CHK_IERR( ierr = MPI_Get_processor_name(myNodeId, &nodeId_strlen), ierr);
    char allNodeId[nRanks][MPI_MAX_PROCESSOR_NAME];
    PHIST_CHK_IERR( ierr = MPI_Allgather(myNodeId, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, allNodeId, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, MPI_COMM_WORLD), ierr);
    // on the same node the node ids are identical
    std::map<std::string,int> nodeNameSet;
    for(int i = 0; i < nRanks; i++)
    {
      if( i == myRank )
        myRankInNode = nodeNameSet[(std::string(myNodeId))];
      nodeNameSet[(std::string(myNodeId))]++;
    }
    ranksPerNode = nRanks / nodeNameSet.size();
  }
  int myCPU_SETSIZE = std::max(CPU_SETSIZE, ranksPerNode*nThreads);

  PHIST_SOUT(PHIST_WARNING,"Trying to pin the OpenMP threads to the cores based on simple heuristics (may fail, use GHOST instead!).\n"\
                           "Using %d ranks per node with %d threads each\n",
                           ranksPerNode, nThreads);

#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    int targetCoreId = nThreads*myRankInNode + tid;

    cpu_set_t *cpu = CPU_ALLOC(myCPU_SETSIZE);
    size_t size = CPU_ALLOC_SIZE(myCPU_SETSIZE);
    CPU_ZERO_S(size, cpu);
    CPU_SET_S(targetCoreId, size, cpu);
    sched_setaffinity(0, size, cpu);
    CPU_FREE(cpu);
  }

  // print result
  std::ostringstream oss;
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    int coreid = sched_getcpu();
#pragma omp critical
    oss << "\t" << coreid << " (" << tid << ")";
  }
  if( myRank < ranksPerNode )
  {
    PHIST_OUT(PHIST_INFO, "Result of pinning is coreId(threadId): %s\n", oss.str().c_str());
  }
}
}
#endif


// remember if we initialized the MPI
namespace {
int mpiInitializedBefore = false;
}

extern "C" {


// comment in for glibc/gcc and memory alignment problems...
/*
#include <malloc.h>
#include <cstring>
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

  PHIST_CHK_IERR( *ierr = MPI_Initialized(&mpiInitializedBefore), *ierr);
  if (!mpiInitializedBefore)
  {
  	PHIST_CHK_IERR( *ierr = MPI_Init(argc, argv), *ierr);
  }

  // allow unlimited stack
  struct rlimit rlim = { RLIM_INFINITY, RLIM_INFINITY };
  PHIST_CHK_IERR( *ierr = setrlimit(RLIMIT_STACK, &rlim), *ierr);;

#ifdef PHIST_KERNEL_LIB_BUILTIN_PIN_THREADS
  pinThreads();
#endif

#ifdef PHIST_HAVE_LIKWID
  LIKWID_MARKER_INIT;
#pragma omp parallel
  {
    LIKWID_MARKER_THREADINIT;
    LIKWID_MARKER_START("phist<builtin>");
  }
#endif
  init_random_seed();
}

// finalize builtin kernels
void phist_kernels_finalize(int* ierr)
{
#ifdef PHIST_HAVE_LIKWID
#pragma omp parallel
  {
    LIKWID_MARKER_STOP("phist<builtin>");
  }
  LIKWID_MARKER_CLOSE;
#endif
PHIST_CXX_TIMER_SUMMARIZE
  if( !mpiInitializedBefore )
  {
    PHIST_CHK_IERR( *ierr = MPI_Finalize(), *ierr);
  }
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
  PHIST_ENTER_FCN(__FUNCTION__);
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


