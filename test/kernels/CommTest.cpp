/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_config.h"

#include "phist_tools.h"
#include "kernels/phist_kernels.h"

#include "gtest/phist_gtest.h"

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

#ifdef PHIST_HAVE_MPI

/* obtain an MPI_Comm from a phist_comm */
TEST_F(XCommTest,get_mpi_comm)
{
  MPI_Comm mpi_comm;
  phist_comm_get_mpi_comm(comm_, &mpi_comm, &iflag_);
  ASSERT_EQ(0,iflag_);
  EXPECT_EQ(MPI_COMM_WORLD, mpi_comm);
}

  /*! Test the comm_get_rank function - is the comm in the kernel lib really MPI_COMM_WORLD?. */
  TEST_F(XCommTest, change_default_comm_to_SELF)
  {
    phist_comm_ptr comm=nullptr;
    phist_comm_create(&comm,&iflag_);
    ASSERT_EQ(0,iflag_);

    int rank, size;
    phist_comm_get_rank(comm,&rank,&iflag_);
    ASSERT_EQ(0,iflag_);
    phist_comm_get_size(comm,&size,&iflag_);
    ASSERT_EQ(0,iflag_);

    // the default comm should be MPI_COMM_WORLD
    ASSERT_EQ(mpi_rank_,rank);
    ASSERT_EQ(mpi_size_,size);

    phist_comm_delete(comm,&iflag_);

    // now change it to MPI_COMM_SELF
    phist_set_default_comm(MPI_COMM_SELF);

    comm=nullptr;
    phist_comm_create(&comm,&iflag_);
    ASSERT_EQ(0,iflag_);

    phist_comm_get_rank(comm,&rank,&iflag_);
    ASSERT_EQ(0,iflag_);
    phist_comm_get_size(comm,&size,&iflag_);
    ASSERT_EQ(0,iflag_);

    // this comm has size 1, so everyone is rank 0
    ASSERT_EQ(0,rank);
    ASSERT_EQ(1,size);

    phist_comm_delete(comm,&iflag_);    

    // reset the comm - otherwise subsequent tests fail...
    phist_set_default_comm(MPI_COMM_WORLD);
  }
#endif /* PHIST_HAVE_MPI */
