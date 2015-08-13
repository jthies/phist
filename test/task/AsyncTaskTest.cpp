#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"

#include "phist_kernels.h"
#include "phist_macros.h"

#include <vector>
#include <unistd.h>


// small helper function to observe timing behavior
void doSomething()
{
  usleep(100000);
}

class AsyncTaskTest: public testing::Test
{
  public:

    virtual void SetUp()
    {
      iflag_ = 0;
      iflag_t1_ = 0;
      iflag_t2_ = 0;
      iflag_t3_ = 0;
    }

    virtual void TearDown()
    {
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(0,iflag_t1_);
      ASSERT_EQ(0,iflag_t2_);
      ASSERT_EQ(0,iflag_t3_);
    }

  protected:
    int iflag_;
    int iflag_t1_;
    int iflag_t2_;
    int iflag_t3_;
};


// just make sure, tasks are executed
TEST_F(AsyncTaskTest, execute_some_tasks)
{
PHIST_TASK_DECLARE(Task1)
PHIST_TASK_DECLARE(Task2)
PHIST_TASK_DECLARE(Task3)

  bool task1Run = false;
  bool task2Run = false;
  bool task3Run = false;

PHIST_TASK_BEGIN(Task1)
  // to verify we do not simply observe timing behavior!
  doSomething();
  task1Run = true;
PHIST_TASK_END_NOWAIT(&iflag_)

PHIST_TASK_BEGIN(Task2)
  task2Run = true;
PHIST_TASK_END_NOWAIT(&iflag_)

PHIST_TASK_BEGIN(Task3)
  task3Run = true;
PHIST_TASK_END_NOWAIT(&iflag_)

PHIST_TASK_WAIT(Task3,&iflag_)
PHIST_TASK_WAIT(Task2,&iflag_)
PHIST_TASK_WAIT(Task1,&iflag_)

  ASSERT_TRUE(task1Run);
  ASSERT_TRUE(task2Run);
  ASSERT_TRUE(task3Run);
}


// check that tasks are executed asynchronously
#ifdef PHIST_KERNEL_LIB_GHOST
TEST_F(AsyncTaskTest, is_async)
#else
TEST_F(AsyncTaskTest, DISABLED_is_async)
#endif
{
PHIST_TASK_DECLARE(Task1)

  std::vector<int> order(2,-1);
  int step = 0;

PHIST_TASK_BEGIN(Task1)
  // just use a timer here so the creating thread comes first - this *is* a race condition!
  doSomething();
#pragma omp atomic capture
  order[1] = step++;
PHIST_TASK_END_NOWAIT(&iflag_)

#pragma omp atomic capture
  order[0] = step++;

PHIST_TASK_WAIT(Task1,&iflag_)

  std::vector<int> expectedOrder;
  expectedOrder.push_back(0);
  expectedOrder.push_back(1);
  ASSERT_EQ(expectedOrder, order);
}


// check that async tasks can start compute tasks (e.g. non-async tasks)
TEST_F(AsyncTaskTest, nested_blocking_task)
{
PHIST_TASK_DECLARE(Task1)
PHIST_TASK_DECLARE(Task2)

  bool task1Started = false;
  bool task2Run = false;
  bool task1Finished = false;

PHIST_TASK_BEGIN(Task1)
  task1Started = true;
PHIST_TASK_BEGIN(Task2)
  task2Run = true;
PHIST_TASK_END(&iflag_t1_)
  task1Finished = true;
PHIST_TASK_END_NOWAIT(&iflag_)

PHIST_TASK_WAIT(Task1,&iflag_)

  ASSERT_TRUE(task1Started);
  ASSERT_TRUE(task2Run);
  ASSERT_TRUE(task1Finished);
}


// check nesting of async task
TEST_F(AsyncTaskTest, nested_async_task)
{
PHIST_TASK_DECLARE(Task1)
PHIST_TASK_DECLARE(Task2)

  bool task1Started = false;
  bool task2Run = false;
  bool task1Finished = false;

PHIST_TASK_BEGIN(Task1)
  task1Started = true;
PHIST_TASK_BEGIN(Task2)
  task2Run = true;
PHIST_TASK_END_NOWAIT(&iflag_t1_)
  task1Finished = true;
PHIST_TASK_END_NOWAIT(&iflag_)

PHIST_TASK_WAIT(Task1,&iflag_)
PHIST_TASK_WAIT(Task2,&iflag_)

  ASSERT_TRUE(task1Started);
  ASSERT_TRUE(task2Run);
  ASSERT_TRUE(task1Finished);
}


// check if the task semaphore works
#ifdef PHIST_KERNEL_LIB_GHOST
TEST_F(AsyncTaskTest, task_post_wait_step)
#else
TEST_F(AsyncTaskTest, DISABLED_task_post_wait_step)
#endif
{
PHIST_TASK_DECLARE(Task1)

  std::vector<int> order(4,-1);
  int step = 0;

// start an async task
PHIST_TASK_BEGIN(Task1)

  doSomething();
#pragma omp atomic capture
  order[0] = step++;
PHIST_TASK_POST_STEP(&iflag_t1)


  doSomething();
#pragma omp atomic capture
  order[2] = step++;
PHIST_TASK_END_NOWAIT(&iflag_)


PHIST_TASK_WAIT_STEP(Task1,&iflag_)
#pragma omp atomic capture
  order[1] = step++;

PHIST_TASK_WAIT(Task1,&iflag_)

#pragma omp atomic capture
  order[3] = step++;

  std::vector<int> expectedOrder;
  expectedOrder.push_back(0);
  expectedOrder.push_back(1);
  expectedOrder.push_back(2);
  expectedOrder.push_back(3);
  ASSERT_EQ(expectedOrder, order);
}


// check ordering of async and compute tasks with the task semaphore
#ifdef PHIST_KERNEL_LIB_GHOST
TEST_F(AsyncTaskTest, simulate_communication_hiding)
#else
TEST_F(AsyncTaskTest, DISABLED_simulate_communication_hiding)
#endif
{
PHIST_TASK_DECLARE(Task1)
PHIST_TASK_DECLARE(Task3)

  std::vector<int> order(4,-1);
  int step = 0;

// start an async task
PHIST_TASK_BEGIN(Task1)

PHIST_TASK_DECLARE(Task2)
PHIST_TASK_BEGIN(Task2)
  doSomething(); // local calculation
#pragma omp atomic capture
  order[0] = step++;
PHIST_TASK_END(&iflag_t1_)

PHIST_TASK_POST_STEP(&iflag_t1_)

  doSomething(); // communication
#pragma omp atomic capture
  order[2] = step++;
PHIST_TASK_END_NOWAIT(&iflag_)


// wait till calculation is finished
PHIST_TASK_WAIT_STEP(Task1,&iflag_)
PHIST_TASK_BEGIN(Task3)
#pragma omp atomic capture
  order[1] = step++;
  doSomething(); // other calculation to hide communication
PHIST_TASK_END(&iflag_)

PHIST_TASK_WAIT(Task1,&iflag_)

#pragma omp atomic capture
  order[3] = step++;

  std::vector<int> expectedOrder;
  expectedOrder.push_back(0);
  expectedOrder.push_back(1);
  expectedOrder.push_back(2);
  expectedOrder.push_back(3);
  ASSERT_EQ(expectedOrder, order);
}


// anonymous namespace
namespace
{
  const int maxLevel = 4;
  const int recursiveRespawns = 7;
  const int totalTasks = 847;
  bool globalErr = false;
  int totalTasksRun = 0;

  // helper class for stress test
  class SpawnTasks
  {
    private:
      PHIST_TASK_DECLARE(task_);
      int ierr_;
      int workId_, level_, nextLevel_;
      bool executed_;
      std::vector<SpawnTasks> childTasks_;

      // disallow copy constructor etc
#ifdef PHIST_HAVE_CXX11_MOVEDEFAULT
      SpawnTasks(const SpawnTasks&) = delete;
      SpawnTasks& operator=(const SpawnTasks&) = delete;
#endif

    public:

#ifdef PHIST_HAVE_CXX11_MOVEDEFAULT
      // move constructor
      SpawnTasks(SpawnTasks&&) = default;
      SpawnTasks& operator=(const SpawnTasks&) = default;
#else
      SpawnTasks(const SpawnTasks&) = default;
#endif
      SpawnTasks()
      {
        ierr_ = 0;
        executed_ = false;
        workId_ = -1;
        level_ = -1;

      }
      void runTask(int workId, int level)
      {
        workId_ = workId;
        level_ = level;
        nextLevel_ = level+1;
        if( level_ > maxLevel )
        {
          executed_ = true;
          return;
        }
        if( level_ == maxLevel && workId_ % 2 ) // TODO: for now only use blocking tasks at the highest level!
        {
          PHIST_TASK_BEGIN(task_);
          if( workId_ % 3 == 0)
          {
            childTasks_.resize(recursiveRespawns);
            for(int i = 0; i < recursiveRespawns; i++)
              childTasks_[i].runTask(i, nextLevel_);
          }
#pragma omp atomic update
          totalTasksRun++;
          if( workId_ % 4 == 1 && level_ % 4 == 2)
            usleep(1000);
          childTasks_.clear();
          executed_ = true;
          PHIST_TASK_END(&ierr_);
        }
        else
        {
          PHIST_TASK_BEGIN(task_);
          if( workId_ % 3 == 0)
          {
            childTasks_.resize(recursiveRespawns);
            for(int i = 0; i < recursiveRespawns; i++)
              childTasks_[i].runTask(i, nextLevel_);
          }
#pragma omp atomic update
          totalTasksRun++;
          if( workId_ % 4 == 0 && level_ % 2 == 0)
            usleep(1000);
          executed_ = true;

          // wait for child tasks
          childTasks_.clear();

          PHIST_TASK_END_NOWAIT(&ierr_);

          // wait for this task to run before leaving scope!
          PHIST_TASK_WAIT(task_,&ierr_);
        }
      }
      ~SpawnTasks()
      {
        if( ierr_ != 0 || !executed_ )
        {
#pragma omp atomic update
          globalErr |= true;
        }
      }
  };
}

#ifdef PHIST_HAVE_CXX11_MOVEDEFAULT
TEST_F(AsyncTaskTest, DISABLED_stressTest)
{
  // just start a crazy bunch of tasks and verify that all run
  globalErr = false;
  totalTasksRun = 0;
  PHIST_SOUT(PHIST_INFO,"Starting tasking stress test...\n");
  {
    std::vector<SpawnTasks> rootTasks;
    rootTasks.resize(recursiveRespawns);
    for(int i = 0; i < recursiveRespawns; i++)
      rootTasks[i].runTask(i, 0);
    rootTasks.clear();
  }
  PHIST_SOUT(PHIST_INFO,"Finished tasking stress test (%d tasks executed)\n", totalTasksRun);
  EXPECT_FALSE(globalErr);
  EXPECT_EQ(totalTasks,totalTasksRun);
}
#endif
