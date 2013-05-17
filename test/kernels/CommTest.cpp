
#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "kernels/essex_kernels.h"
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
    ierr_=0;
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
    essex_comm_create(&essexComm_,&ierr_);
    }

  /*! Clean up.
   */
  virtual void TearDown() 
    {
    KernelTest::TearDown();
    }

  MPI_Comm mpiComm_;
  int rank_, size_;
  comm_ptr_t essexComm_;

};

  /*! Test the comm_get_rank function. */
  TEST_F(CommTest, get_rank) {
        int rank;
        essex_comm_get_rank(essexComm_,&rank,&ierr_);
	ASSERT_EQ(ierr_, 0);
	ASSERT_EQ(rank_, rank);
}

  /*! Test the comm_get_size function. */
  TEST_F(CommTest, get_size) {
        int size;
        essex_comm_get_size(essexComm_,&size,&ierr_);
	ASSERT_EQ(ierr_, 0);
	ASSERT_EQ(size_, size);
}
