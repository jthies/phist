/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/phist_gtest.h"

#include "phist_kernels.h"
#include "phist_macros.h"

#ifdef PHIST_HAVE_GHOST
#include "ghost/taskq.h"
#include "ghost/task.h"
#endif


class MainTaskTest: public testing::Test
{
  public:

    virtual void SetUp()
    {
      curTask_ = NULL;
#ifdef PHIST_KERNEL_LIB_GHOST
      ghost_task_cur(&curTask_);
#endif
    }

    virtual void TearDown()
    {
      curTask_ = NULL;
    }

  protected:
#ifdef PHIST_KERNEL_LIB_GHOST
    ghost_task* curTask_;
#else
    void* curTask_;
#endif
};


// check that we are inside a task
#ifdef PHIST_KERNEL_LIB_GHOST
TEST_F(MainTaskTest, inside_task)
#else
TEST_F(MainTaskTest, DISABLED_inside_task)
#endif
{
  ASSERT_TRUE(curTask_ != NULL);
}


// check that the task has no parent
#ifdef PHIST_KERNEL_LIB_GHOST
TEST_F(MainTaskTest, root_task)
#else
TEST_F(MainTaskTest, DISABLED_root_task)
#endif
{
  ASSERT_TRUE(curTask_ != NULL);
#ifdef PHIST_KERNEL_LIB_GHOST
  ASSERT_TRUE(curTask_->parent == NULL);
#else
  FAIL();
#endif
}


// check that the task has requires no resources
#ifdef PHIST_KERNEL_LIB_GHOST
TEST_F(MainTaskTest, task_has_no_threads)
#else
TEST_F(MainTaskTest, DISABLED_task_has_no_threads)
#endif
{
  ASSERT_TRUE(curTask_ != NULL);
#ifdef PHIST_KERNEL_LIB_GHOST
  PHIST_OUT(PHIST_DEBUG, " has nThreads = %d\n", curTask_->nThreads);
  ASSERT_EQ(curTask_->nThreads, 0);
#else
  FAIL();
#endif
}


// check that the task currently runs
#ifdef PHIST_KERNEL_LIB_GHOST
TEST_F(MainTaskTest, task_is_running)
#else
TEST_F(MainTaskTest, DISABLED_task_is_running)
#endif
{
  ASSERT_TRUE(curTask_ != NULL);
#ifdef PHIST_KERNEL_LIB_GHOST
  ASSERT_EQ(GHOST_TASK_RUNNING, curTask_->state);
#else
  FAIL();
#endif
}


