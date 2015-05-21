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

  bool task1Run = false;
  bool task2Run = false;
  bool task3Run = false;

PHIST_TASK_BEGIN(Task1)
  task1Run = true;
PHIST_TASK_END(&iflag_)

  ASSERT_TRUE(task1Run);
  ASSERT_FALSE(task2Run);
  ASSERT_FALSE(task3Run);

PHIST_TASK_BEGIN(Task2)
  task2Run = true;
PHIST_TASK_END(&iflag_)

  ASSERT_TRUE(task1Run);
  ASSERT_TRUE(task2Run);
  ASSERT_FALSE(task3Run);

PHIST_TASK_BEGIN(Task3)
  task3Run = true;
PHIST_TASK_END(&iflag_)

  ASSERT_TRUE(task1Run);
  ASSERT_TRUE(task2Run);
  ASSERT_TRUE(task3Run);
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

  std::vector<int> expectedOrder;
  expectedOrder.push_back(1);
  expectedOrder.push_back(2);
  expectedOrder.push_back(3);
  expectedOrder.push_back(4);
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

  phist_Dmvec_delete(mvec,&iflag_);
  ASSERT_EQ(0, iflag_);
  phist_map_delete(map,&iflag_);
  ASSERT_EQ(0, iflag_);
  phist_comm_delete(comm,&iflag_);
  ASSERT_EQ(0, iflag_);
}


