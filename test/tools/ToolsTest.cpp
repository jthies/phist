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
#include "gtest/phist_gtest.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

/*! Test fixure. Note the "X" prefix so we can filter on tests that are not
    specific to any data type.
 */
class XToolsTest: public ::testing::Test 
{

  public:

  int mpi_rank_, mpi_size_;

  /*! Set up routine.
   */
  virtual void SetUp()
  {
    mpi_rank_=0;
    mpi_size_=1;
    int iflag_;

#ifdef PHIST_HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank_);
    MPI_Comm_size(MPI_COMM_WORLD,&mpi_size_);
#endif
  }

  /*! Clean up.
   */
  virtual void TearDown() 
  {
  }
  

};

  /*! Test the comm_get_rank function - is the comm in the kernel lib really MPI_COMM_WORLD?. */
  TEST_F(XToolsTest, change_default_output_stream)
  {
    std::ostringstream oss, oss_expect;
    phist_set_CXX_output_stream(oss);
    PHIST_OUT(0,"Hello, World!");
    oss_expect << "PE"<<mpi_rank_<<": Hello, World!";
    ASSERT_EQ(oss_expect.str(),oss.str());
  }
