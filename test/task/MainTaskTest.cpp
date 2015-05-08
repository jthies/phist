#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"

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
    ghost_task_t* curTask_;
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


// check that the task has at least some resources
#ifdef PHIST_KERNEL_LIB_GHOST
TEST_F(MainTaskTest, task_has_threads)
#else
TEST_F(MainTaskTest, DISABLED_task_has_threads)
#endif
{
  ASSERT_TRUE(curTask_ != NULL);
#ifdef PHIST_KERNEL_LIB_GHOST
  ASSERT_GT(curTask_->nThreads, 0);
#else
  FAIL();
#endif
}


// check that the task currently runs
#ifdef PHIST_KERNEL_LIB_GHOST
TEST_F(MainTaskTest, task_is_running)
#else
TEST_F(MainTaskTest, DISABLED_task_has_threads)
#endif
{
  ASSERT_TRUE(curTask_ != NULL);
#ifdef PHIST_KERNEL_LIB_GHOST
  ASSERT_EQ(GHOST_TASK_RUNNING, curTask_->state);
#else
  FAIL();
#endif
}


