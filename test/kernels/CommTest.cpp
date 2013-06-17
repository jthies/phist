
#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "phist_kernels.h"
#include "KernelTest.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

using namespace testing;

/*! Test fixure. */
class CommTest: public KernelTest {

public:

  /*! Set up routine.
   */
  virtual void SetUp()
    {
    KernelTest::SetUp();
#ifdef HAVE_MPI
    mpiComm_=MPI_COMM_WORLD;
    ierr_=MPI_Comm_rank(mpiComm_,&rank_);
    ierr_=MPI_Comm_size(mpiComm_,&size_);
#else
    mpiComm_=0;
    rank_=0;
    size_=1;
#endif
    }

  /*! Clean up.
   */
  virtual void TearDown() 
    {
    KernelTest::TearDown();
    }

#ifdef HAVE_MPI
  MPI_Comm mpiComm_;
#else
  int mpiComm_;
#endif
  int rank_, size_;

};

  /*! Test the comm_get_rank function. */
  TEST_F(CommTest, get_rank) {
        int rank;
        phist_comm_get_rank(comm_,&rank,&ierr_);
	ASSERT_EQ(0,ierr_);
	ASSERT_EQ(rank_, rank);
}

  /*! Test the comm_get_size function. */
  TEST_F(CommTest, get_size) {
        int size;
        phist_comm_get_size(comm_,&size,&ierr_);
	ASSERT_EQ(0,ierr_);
	ASSERT_EQ(size_, size);
	ASSERT_EQ(size_, 1); // parallel testing not supported, yet
}
