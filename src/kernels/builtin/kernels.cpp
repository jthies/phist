/*! \file builtin/kernels.cpp
 * wraps implementation of builtin kernel routines
 * \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
 * \author "Jonas Thies <Jonas.Thies@DLR.de>
 *
*/

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#else
#error "builtin kernels only work with MPI"
#endif

#include <cstdlib>
#include <cstdio>
#include <iostream>

#include "phist_macros.h"
#include "phist_kernel_perfmodels.hpp"
#ifdef PHIST_HAVE_TEUCHOS
#include "phist_trilinos_macros.h"
#endif
#include "../phist_kernels.h"

#include "phist_typedefs.h"
#include "typedefs.hpp"
#include "phist_ScalarTraits.hpp"
extern "C" {
#include "bench_kernels.h"
}

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
void omp_set_num_threads(int nt) {return;}
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
  int iflag = 0;
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
  PHIST_CHK_IERR( iflag = MPI_Comm_rank(MPI_COMM_WORLD, &myRank), iflag);
  PHIST_CHK_IERR( iflag = MPI_Comm_size(MPI_COMM_WORLD, &nRanks), iflag);
  // nodes
  int ranksPerNode = 1;
  int myRankInNode = 0;
  {
    char myNodeId[MPI_MAX_PROCESSOR_NAME];
    int nodeId_strlen = 0;
    PHIST_CHK_IERR( iflag = MPI_Get_processor_name(myNodeId, &nodeId_strlen), iflag);
    char allNodeId[nRanks][MPI_MAX_PROCESSOR_NAME];
    PHIST_CHK_IERR( iflag = MPI_Allgather(myNodeId, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, allNodeId, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, MPI_COMM_WORLD), iflag);
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

  // pin the main thread
  {
    int targetCoreId = nThreads*myRankInNode;

    cpu_set_t *cpu = CPU_ALLOC(myCPU_SETSIZE);
    size_t size = CPU_ALLOC_SIZE(myCPU_SETSIZE);
    CPU_ZERO_S(size, cpu);
    CPU_SET_S(targetCoreId, size, cpu);
    sched_setaffinity(0, size, cpu);
    CPU_FREE(cpu);
  }

  // pin the OpenMP threads
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


/*
// comment in for glibc/gcc and memory alignment problems...

#include <malloc.h>
#include <cstring>
// force all memory allocations to be 32 byte aligned!
static void*(*old_malloc_hook)(size_t,const void*);
static void * my_malloc_hook (size_t size, const void *caller)
{
   void *result;
   // Restore all old hooks
   __malloc_hook = old_malloc_hook;
   // Call recursively
   //result = memalign (32, size);
   posix_memalign(&result, 32, size);
   // Save underlying hooks
   old_malloc_hook = __malloc_hook;
   // printf might call malloc, so protect it too.
   if( (uintptr_t)result % 32 != 0 )
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
void phist_kernels_init(int* argc, char*** argv, int* iflag)
{
  *iflag=0;

  PHIST_CHK_IERR( *iflag = MPI_Initialized(&mpiInitializedBefore), *iflag);
  if (!mpiInitializedBefore)
  {
    PHIST_CHK_IERR( *iflag = MPI_Init(argc, argv), *iflag);
  }

  // allow unlimited stack
  struct rlimit rlim = { RLIM_INFINITY, RLIM_INFINITY };
  PHIST_CHK_IERR( *iflag = setrlimit(RLIMIT_STACK, &rlim), *iflag);;

  // set number of OpenMP threads to $PHIST_NUM_THREADS if it is set
  const char* PHIST_NUM_THREADS=getenv("PHIST_NUM_THREADS");
  int num_threads= (PHIST_NUM_THREADS==NULL)? -1:atoi(PHIST_NUM_THREADS);
  if (num_threads>0)
  {
    omp_set_num_threads(num_threads);
  }
  

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
void phist_kernels_finalize(int* iflag)
{
#ifdef PHIST_HAVE_LIKWID
#pragma omp parallel
  {
    LIKWID_MARKER_STOP("phist<builtin>");
  }
  LIKWID_MARKER_CLOSE;
#endif
#if defined(PHIST_TIMEMONITOR) || defined(PHIST_TIMEMONITOR_PERLINE)
PHIST_CXX_TIMER_SUMMARIZE;
#endif
PHIST_PERFCHECK_SUMMARIZE(PHIST_INFO);
  if( !mpiInitializedBefore )
  {
    PHIST_CHK_IERR( *iflag = MPI_Finalize(), *iflag);
  }
#ifdef PHIST_PERFCHECK
  // prevent some strange memory errors during deallocation (due to shared lib context?)
  phist_PerfCheck::benchmarks.clear();
#endif
  *iflag=0;
}

void phist_comm_create(comm_ptr_t* vcomm, int* iflag);
void phist_comm_delete(comm_ptr_t vcomm, int* iflag);
void phist_comm_get_rank(const_comm_ptr_t vcomm, int* rank, int* iflag);
void phist_comm_get_size(const_comm_ptr_t vcomm, int* size, int* iflag);
#ifdef PHIST_HAVE_MPI
void phist_comm_get_mpi_comm(const_comm_ptr_t comm, MPI_Comm* mpiComm, int* iflag)
{
  *iflag=0;
  MPI_Fint fcomm = *((MPI_Fint*)comm);
  *mpiComm = MPI_Comm_f2c(fcomm);
}
#endif
void phist_map_create(map_ptr_t* vmap, const_comm_ptr_t vcomm, gidx_t nglob, int *iflag);
void phist_map_delete(map_ptr_t vmap, int *iflag);
void phist_map_get_comm(const_map_ptr_t vmap, const_comm_ptr_t* vcomm, int* iflag);
void phist_map_get_local_length(const_map_ptr_t vmap, lidx_t* nloc, int* iflag);
void phist_map_get_global_length(const_map_ptr_t vmap, gidx_t* nglob, int* iflag);
void phist_map_get_ilower(const_map_ptr_t vmap, gidx_t* ilower, int* iflag);
void phist_map_get_iupper(const_map_ptr_t vmap, gidx_t* iupper, int* iflag);

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
  double total_bw = *max_bw;
#ifdef PHIST_HAVE_MPI
  PHIST_CHK_IERR(*iflag = MPI_Reduce(max_bw,&total_bw,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD), *iflag);
#endif
  PHIST_SOUT(PHIST_INFO, "measured %8.4g Gb/s per process, %8.4g Gb/s total\n", *max_bw/1.e9, total_bw/1.e9);
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
  double total_bw = *max_bw;
#ifdef PHIST_HAVE_MPI
  PHIST_CHK_IERR(*iflag = MPI_Reduce(max_bw,&total_bw,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD), *iflag);
#endif
  PHIST_SOUT(PHIST_INFO, "measured %8.4g Gb/s per process, %8.4g Gb/s total\n", *max_bw/1.e9, total_bw/1.e9);
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
  double total_bw = *max_bw;
#ifdef PHIST_HAVE_MPI
  PHIST_CHK_IERR(*iflag = MPI_Reduce(max_bw,&total_bw,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD), *iflag);
#endif
  PHIST_SOUT(PHIST_INFO, "measured %8.4g Gb/s per process, %8.4g Gb/s total\n", *max_bw/1.e9, total_bw/1.e9);
}



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

