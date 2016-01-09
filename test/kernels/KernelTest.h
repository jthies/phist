#ifndef PHIST_KERNEL_TEST_H
#define PHIST_KERNEL_TEST_H

#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#else
typedef int MPI_Comm;
#endif

#include "gtest/gtest.h"

#include "phist_typedefs.h"
#include "phist_kernels.h"


#ifdef PHIST_HAVE_GHOST
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
class KernelTest: public testing::Test {

public:

 comm_ptr_t comm_;
 bool haveS_, haveD_, haveC_, haveZ_;
 MPI_Comm mpi_comm_;
 unsigned int rseed_;//random number seed
 int iflag_, mpi_rank_, mpi_size_;
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
 
	/** Set up method.
	 */
	virtual void SetUp() 
	{
    if( kernelTestSetupCounter_++ == 0 )
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
#ifdef PHIST_HAVE_SP
      phist_Stype_avail(&iflag_); haveS_=(iflag_==0);
      phist_Ctype_avail(&iflag_); haveC_=(iflag_==0);
#else
      haveS_=false;
      haveC_=false;
#endif
      phist_Dtype_avail(&iflag_); haveD_=(iflag_==0);
      phist_Ztype_avail(&iflag_); haveZ_=(iflag_==0);

      // initialize random number sequence in a reproducible way (yet
      // with a different seed on each MPI process)
      rseed_ = (unsigned int)(mpi_rank_*77+42);
      srand(rseed_);
    }
	}

virtual void TearDown()
  {
  // should work if we count correctly
  if( --kernelTestSetupCounter_ == 0 )
    {
    phist_comm_delete(comm_,&iflag_);
    ASSERT_EQ(0,iflag_);
    comm_=NULL;
    }
#ifdef PHIST_HAVE_GHOST
    ghost_taskq_waitall();
#endif  
  }
  
::testing::AssertionResult AssertNotNull(void* ptr)
  {
  if (ptr==NULL) return ::testing::AssertionFailure();
  return ::testing::AssertionSuccess();
  }

private:
 int kernelTestSetupCounter_;
};

#endif
