#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"


#include "phist_kernels.h"
#include "KernelTest.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

using namespace testing;

/*! Test fixure. */
class XCommTest: public KernelTest {

public:

  /*! Set up routine.
   */
  virtual void SetUp()
    {
    KernelTest::SetUp();
    }

  /*! Clean up.
   */
  virtual void TearDown() 
    {
    KernelTest::TearDown();
    }

};

  /*! Test the comm_get_rank function - is the comm in the kernel lib really MPI_COMM_WORLD?. */
  TEST_F(XCommTest, get_rank) {
        int rank=42;
        phist_comm_get_rank(comm_,&rank,&iflag_);
	ASSERT_EQ(0,iflag_);
	ASSERT_EQ(mpi_rank_, rank);
}

  /*! Test the comm_get_size function. Is the comm in the kernel lib really MPI_COMM_WORLD? */
  TEST_F(XCommTest, get_size) {
        int size=42;
        phist_comm_get_size(comm_,&size,&iflag_);
	ASSERT_EQ(0,iflag_);
	ASSERT_EQ(mpi_size_, size);
}
