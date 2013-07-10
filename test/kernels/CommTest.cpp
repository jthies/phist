
#include "gtest/gtest.h"
//#include "gmock/gmock.h"

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
    std::cout << "in CommTest::SetUp"<<std::endl;
    KernelTest::SetUp();
    std::cout << "my rank: "<< mpi_rank_<<std::endl;
    }

  /*! Clean up.
   */
  virtual void TearDown() 
    {
    KernelTest::TearDown();
    }

};

  /*! Test the comm_get_rank function - is the comm in the kernel lib really MPI_COMM_WORLD?. */
  TEST_F(CommTest, get_rank) {
        int rank=42;
        phist_comm_get_rank(comm_,&rank,&ierr_);
	ASSERT_EQ(0,ierr_);
	ASSERT_EQ(mpi_rank_, rank);
}

  /*! Test the comm_get_size function. Is the comm in the kernel lib really MPI_COMM_WORLD? */
  TEST_F(CommTest, get_size) {
        int size=42;
        phist_comm_get_size(comm_,&size,&ierr_);
	ASSERT_EQ(0,ierr_);
	ASSERT_EQ(mpi_size_, size);
}
