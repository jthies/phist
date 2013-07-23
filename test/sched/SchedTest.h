#ifndef PHIST_SCHED_TEST_H
#define PHIST_SCHED_TEST_H

#include <iostream>
#include <iomanip>
#include "gtest/gtest.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#else
typedef int MPI_Comm;
#endif

#include "phist_typedefs.h"
#include "phist_kernels.h"
#include "ghost/ghost_taskq.h"

/** Base for tests involving dynamic scheduling.
    Right now this class only initializes the ghost task-queue.
 */
class SchedTest: public testing::Test {
public:

	/** Set up method.
	calls ghost_taskq_init()
	 */
	virtual void SetUp() 
	{
	// if ghost is the kernel lib, this is done at the beginning of main
	// in phist_kernels_init()
#ifndef PHIST_KERNEL_LIB_GHOST
	num_ghost_queues_=1;
        ghost_cpuid_init();
        ghost_thpool_init(ghost_cpuid_topology.numHWThreads/ghost_cpuid_topology.numThreadsPerCore);
	ghost_taskq_init(num_ghost_queues_)
#else
        num_ghost_queues_=ghost_cpuid_topology.numSockets;
#endif
	}

virtual void TearDown()
  {
  ghost_taskq_finish();
  }
  

int num_ghost_queues_;
bool myMpiSession_;
};

#endif
