/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef PHIST_KERNEL_TEST_H
#define PHIST_KERNEL_TEST_H

#include "phist_config.h"

#include "phist_tools.h"
#include "phist_kernels.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#else
typedef int MPI_Comm;
#endif

#include "gtest/phist_gtest.h"


#ifdef PHIST_HAVE_GHOST
#include "ghost/util.h"
#include "ghost/taskq.h"
#endif
#ifdef PHIST_MVECS_ROW_MAJOR
#define VIDX(i,j,lda) ((j)+(i)*(lda))
#else
#define VIDX(i,j,lda) ((i)+(j)*(lda))
#endif
#ifdef PHIST_SDMATS_ROW_MAJOR
#define MIDX(i,j,lda) ((j)+(i)*(lda))
#else
#define MIDX(i,j,lda) ((i)+(j)*(lda))
#endif
/** Base for tests using kernel operations.
It calls the init and finalize routines of the kernel lib
and provides basic MPI support
 */
class KernelTest: public virtual testing::Test {

public:

 static phist_comm_ptr comm_;
 static bool haveS_, haveD_, haveC_, haveZ_;
 static MPI_Comm mpi_comm_;
 static unsigned int rseed_;//random number seed
 static int iflag_, mpi_rank_, mpi_size_;
 static bool isCuda_;
 KernelTest() : kernelTestSetupCounter_(0) {}

 //! these flags determine how to traverse arrays
 #ifdef PHIST_MVECS_ROW_MAJOR
 static const bool vflag_=true;
 #else
 static const bool vflag_=false;
 #endif
 #ifdef PHIST_SDMATS_ROW_MAJOR
 static const bool mflag_=true;
 #else
 static const bool mflag_=false;
 #endif
 
	/** static setup method for the complete test case
	 */
	static void SetUpTestCase() 
	{
    if( staticKernelTestSetupCounter_++ == 0 )
    {
#ifdef PHIST_HAVE_MPI
      mpi_comm_=MPI_COMM_WORLD;
      phist_comm_create(&comm_,&iflag_);
      ASSERT_EQ(0,iflag_);
      iflag_=MPI_Comm_rank(mpi_comm_,&mpi_rank_);
      ASSERT_EQ(0,iflag_);
      iflag_=MPI_Comm_size(mpi_comm_,&mpi_size_);
      ASSERT_EQ(0,iflag_);
#else
      mpi_comm_=-1;
      mpi_rank_=0;
      mpi_size_=1;
#endif

      haveS_=false;
      haveC_=false;
      haveD_=false;
      haveZ_=false;

#ifdef PHIST_HAVE_SP
      phist_Stype_avail(&iflag_); haveS_=(iflag_==0);
# ifdef PHIST_HAVE_CMPLX
      phist_Ctype_avail(&iflag_); haveC_=(iflag_==0);
# endif
#endif
      phist_Dtype_avail(&iflag_); haveD_=(iflag_==0);
#ifdef PHIST_HAVE_CMPLX
      phist_Ztype_avail(&iflag_); haveZ_=(iflag_==0);
#endif

      // initialize random number sequence in a reproducible way (yet
      // with a different seed on each MPI process)
      rseed_ = (unsigned int)(mpi_rank_*77+42);
      srand(rseed_);
      isCuda_=false;
#ifdef PHIST_KERNEL_LIB_GHOST
      ghost_type type;
      ghost_error gerr=ghost_type_get(&type);
      EXPECT_EQ(GHOST_SUCCESS,gerr);
      isCuda_ = (type==GHOST_TYPE_CUDA);
#endif
    }
	}

  /** dynamic setup method for every single test
   */
  virtual void SetUp()
  {
    kernelTestSetupCounter_++;
  }

  /** dynamic teardown method for every single test
   */
  virtual void TearDown()
  {
    if( --kernelTestSetupCounter_ == 0 )
    {
#ifdef PHIST_HAVE_GHOST
      ghost_taskq_waitall();
#endif  
    }
  }

  /** static teardown method for the complete test case
   */
  static void TearDownTestCase()
  {
  // should work if we count correctly
  if( --staticKernelTestSetupCounter_ == 0 )
    {
    phist_comm_delete(comm_,&iflag_);
    EXPECT_EQ(0,iflag_);
    comm_=NULL;
    }
  }
  
::testing::AssertionResult AssertNotNull(void* ptr)
  {
  if (ptr==NULL) return ::testing::AssertionFailure();
  return ::testing::AssertionSuccess();
  }

private:
  int kernelTestSetupCounter_;
  static int staticKernelTestSetupCounter_;
};

#endif
