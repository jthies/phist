/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
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
#include <string>
#include <map>

#include "phist_tools.h"
#include "phist_kernel_perfmodels.hpp"
#include "../phist_kernels.h"

#include "phist_typedefs.h"
#include "typedefs.hpp"
extern "C" {
#include "phist_random.h"
}

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
#ifdef PHIST_TRY_TO_PIN_THREADS
#include <string>
#include <map>
#include <sstream>
#include <sched.h>

namespace {
// thread pinning using simple heuristics
void pinThreads(bool quiet)
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
  PHIST_CHK_IERR( iflag = MPI_Comm_rank(phist_get_default_comm(), &myRank), iflag);
  PHIST_CHK_IERR( iflag = MPI_Comm_size(phist_get_default_comm(), &nRanks), iflag);
  // nodes
  int ranksPerNode = 1;
  int myRankInNode = 0;
  {
    char *myNodeId=new char[MPI_MAX_PROCESSOR_NAME];
    char *allNodeId=new char[nRanks*MPI_MAX_PROCESSOR_NAME];
    int nodeId_strlen = 0;
    PHIST_CHK_IERR( iflag = MPI_Get_processor_name(myNodeId, &nodeId_strlen), iflag);
    PHIST_CHK_IERR( iflag = MPI_Allgather(myNodeId, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, allNodeId, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, MPI_COMM_WORLD), iflag);
    // on the same node the node ids are identical
    std::map<std::string,int> nodeNameSet;
    for(int i = 0; i < nRanks; i++)
    {
      std::string curNodeIdStr(&(allNodeId[i*MPI_MAX_PROCESSOR_NAME]));
      if( i == myRank ) myRankInNode = nodeNameSet[curNodeIdStr];
      nodeNameSet[curNodeIdStr]++;
    }
    ranksPerNode = nRanks / nodeNameSet.size();
    delete [] myNodeId;
    delete [] allNodeId;
  }
  int myCPU_SETSIZE = std::max(CPU_SETSIZE, ranksPerNode*nThreads);
  if (!quiet)
  {
    PHIST_SOUT(PHIST_WARNING,"Pinning the OpenMP threads to the cores based on simple heuristics (may fail).\n"
                           "Using %d ranks per node with %d threads each\n",
                           ranksPerNode, nThreads);
  }
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
  if( myRank < ranksPerNode && !quiet)
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


// initialize
void phist_kernels_init(int* argc, char*** argv, int* iflag)
{
  bool quiet=(*iflag&PHIST_KERNELS_QUIET);
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
  

#ifdef PHIST_TRY_TO_PIN_THREADS
  pinThreads(quiet);
#endif

#ifdef PHIST_HAVE_LIKWID
  LIKWID_MARKER_INIT;
#pragma omp parallel
  {
    LIKWID_MARKER_THREADINIT;
    LIKWID_MARKER_START("phist<builtin>");
  }
#endif

  PHIST_CHK_IERR(phist_kernels_common_init(argc,argv,iflag),*iflag);
}

// finalize builtin kernels
void phist_kernels_finalize(int* iflag)
{
  PHIST_CHK_IERR(phist_kernels_common_finalize(iflag),*iflag);

  if( !mpiInitializedBefore )
  {
    PHIST_CHK_IERR( *iflag = MPI_Finalize(), *iflag);
  }
  *iflag=0;
}

void phist_comm_create(phist_comm_ptr* vcomm, int* iflag);
void phist_comm_delete(phist_comm_ptr vcomm, int* iflag);
void phist_comm_get_rank(phist_const_comm_ptr vcomm, int* rank, int* iflag);
void phist_comm_get_size(phist_const_comm_ptr vcomm, int* size, int* iflag);
#ifdef PHIST_HAVE_MPI
void phist_comm_get_mpi_comm(phist_const_comm_ptr comm, MPI_Comm* mpiComm, int* iflag)
{
  *iflag=0;
  MPI_Fint fcomm = *((MPI_Fint*)comm);
  *mpiComm = MPI_Comm_f2c(fcomm);
}
#endif
void phist_map_create(phist_map_ptr* vmap, phist_const_comm_ptr vcomm, phist_gidx nglob, int *iflag);
void phist_map_delete(phist_map_ptr vmap, int *iflag);
void phist_map_get_comm(phist_const_map_ptr vmap, phist_const_comm_ptr* vcomm, int* iflag);
void phist_map_get_local_length(phist_const_map_ptr vmap, phist_lidx* nloc, int* iflag);
void phist_map_get_global_length(phist_const_map_ptr vmap, phist_gidx* nglob, int* iflag);
void phist_map_get_ilower(phist_const_map_ptr vmap, phist_gidx* ilower, int* iflag);
void phist_map_get_iupper(phist_const_map_ptr vmap, phist_gidx* iupper, int* iflag);
void phist_maps_compatible(phist_const_map_ptr map1, phist_const_map_ptr map2, int* iflag);

#include "../common/default_context.h"
#include "../common/default_context.cpp"

#include "../common/phist_bench_kernels.cpp"

#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "../common/kernels_no_impl.cpp"

# ifdef PHIST_HAVE_CMPLX
# include "phist_gen_c.h"
# include "../common/kernels_no_impl.cpp"
# endif
#endif

#ifdef PHIST_HAVE_CMPLX
#include "phist_gen_z.h"
#include "../common/kernels_no_impl.cpp"
#endif

} //extern "C"

#include "phist_gen_d.h"
#include "../common/default_context_def.hpp"
#include "kernels_def.hpp"
#include "../common/default_mvec_get_data_def.hpp"
#include "carp_def.hpp"
#include "../common/kernels_no_io.cpp"
#include "../common/kernels_no_gpu.cpp"
#include "../common/kernels_no_fused_spmv.cpp"

