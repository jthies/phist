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

  int iflag_, mpi_rank_, mpi_size_;

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
    // it is crucial to reset the output stream because otherwise
    // the string gets deleted and we get a segfault in the PHIST_OUT.
    phist_set_CXX_output_stream(std::cout);
    ASSERT_EQ(oss_expect.str(),oss.str());
  }

#if defined(PHIST_KERNEL_LIB_GHOST)

#include "phist_ghost_internal.h"

/*! make sure that the PHIST_CHK_GERR macro propagates CUDA errors to all MPI processes if PHIST_GLOBALIZE_CUDA_ERRORS 
is defined.
 */
 TEST_F(XToolsTest, globalize_errors_with_CUDA)
 {
   int num_gpus = phist::ghost_internal::get_num_GPUs();
   iflag_= -mpi_rank_;
   int iflag_expected = iflag_;

#if defined(PHIST_GLOBALIZE_CUDA_ERRORS)
   if (num_gpus>0 && num_gpus<mpi_size_)
   {
     // iflag should be 'global MINned' if it was negative on some proc
     iflag_expected=-(mpi_size_-1);
   }
#endif
   phist::ghost_internal::globalize_cuda_errors(&iflag_);
   ASSERT_EQ(iflag_expected,iflag_);
 }

 TEST_F(XToolsTest, do_not_globalize_warnings_with_CUDA)
 {
   int num_gpus = phist::ghost_internal::get_num_GPUs();
   iflag_= mpi_rank_+1;
   int iflag_expected = iflag_;
   phist::ghost_internal::globalize_cuda_errors(&iflag_);
   ASSERT_EQ(iflag_expected,iflag_);
 }
#endif
