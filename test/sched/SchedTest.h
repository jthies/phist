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
#ifdef PHIST_HAVE_GHOST
extern "C" {
#include "ghost/likwid/cpuid.h"
#include "ghost/ghost_taskq.h"
}
#endif
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
#ifdef PHIST_HAVE_GHOST
	num_ghost_queues_=ghost_cpuid_topology.numSockets;
#endif
	}

virtual void TearDown()
  {
  }
  

int num_ghost_queues_;
};

#endif
