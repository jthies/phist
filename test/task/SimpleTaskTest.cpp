#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"

#include "phist_kernels.h"
#include "phist_macros.h"

#include <vector>

class SimpleTaskTest: public testing::Test
{
  public:

    virtual void SetUp()
    {
      iflag_ = 0;
    }

    virtual void TearDown()
    {
      ASSERT_EQ(0,iflag_);
    }

  protected:
    int iflag_;
};


// just executa some task and check wether they run
TEST_F(SimpleTaskTest, execute_some_tasks)
{
PHIST_TASK_DECLARE(Task1)
PHIST_TASK_DECLARE(Task2)
PHIST_TASK_DECLARE(Task3)

  std::vector<bool> taskRun(3, false);

PHIST_TASK_BEGIN(Task1)
  taskRun[0] = true;
PHIST_TASK_END(&iflag_)

  ASSERT_TRUE(taskRun[0]);
  ASSERT_FALSE(taskRun[1]);
  ASSERT_FALSE(taskRun[2]);

PHIST_TASK_BEGIN(Task2)
  taskRun[1] = true;
PHIST_TASK_END(&iflag_)

  ASSERT_TRUE(taskRun[0]);
  ASSERT_TRUE(taskRun[1]);
  ASSERT_FALSE(taskRun[2]);

PHIST_TASK_BEGIN(Task3)
  taskRun[2] = true;
PHIST_TASK_END(&iflag_)

  ASSERT_TRUE(taskRun[0]);
  ASSERT_TRUE(taskRun[1]);
  ASSERT_TRUE(taskRun[2]);
}


// check that nested blocks of tasks are executed in the expected order
TEST_F(SimpleTaskTest, nested_blocking_tasks)
{
PHIST_TASK_DECLARE(Task1);
PHIST_TASK_DECLARE(Task2);
PHIST_TASK_DECLARE(Task3);

  std::vector<int> order;

PHIST_TASK_BEGIN(Task1)
  order.push_back(1);
  PHIST_TASK_BEGIN(Task3)
  order.push_back(2);
  PHIST_TASK_END(&iflag_)
  order.push_back(3);
PHIST_TASK_END(&iflag_)

PHIST_TASK_BEGIN(Task2)
  order.push_back(4);
PHIST_TASK_END(&iflag_)

  std::vector<int> expectedOrder = {1,2,3,4};
  ASSERT_EQ(expectedOrder, order);
}


// just do some calculations in- and outside off individual tasks
TEST_F(SimpleTaskTest, execute_kernelfcn_in_task)
{
PHIST_TASK_DECLARE(Task1);
PHIST_TASK_DECLARE(Task2);

  comm_ptr_t comm;
  map_ptr_t map = NULL;
  Dmvec_ptr_t mvec = NULL;
  const int nv = 10;
  const int n = 100;
  int tflag = 0;
  const double sqrt_n = 10.;
  double nrm[10];

  phist_comm_create(&comm, &iflag_);
  ASSERT_EQ(0, iflag_);
  phist_map_create(&map,comm,n,&iflag_);
  ASSERT_EQ(0, iflag_);
  phist_Dmvec_create(&mvec,map,nv,&iflag_);

PHIST_TASK_BEGIN(Task1)
  phist_Dmvec_put_value(mvec,1.,&tflag);
PHIST_TASK_END(&iflag_)
  ASSERT_EQ(0, tflag);

  phist_Dmvec_norm2(mvec,nrm,&iflag_);
  ASSERT_EQ(0, iflag_);
  ASSERT_DOUBLE_EQ(sqrt_n, nrm[0]);

  phist_Dmvec_put_value(mvec,2.,&iflag_);
  ASSERT_EQ(0, iflag_);

PHIST_TASK_BEGIN(Task2)
  phist_Dmvec_norm2(mvec,nrm,&tflag);
PHIST_TASK_END(&iflag_)
  ASSERT_EQ(0, tflag);

  ASSERT_DOUBLE_EQ(2*sqrt_n, nrm[0]);
}


