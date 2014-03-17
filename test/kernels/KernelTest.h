#ifndef PHIST_KERNEL_TEST_H
#define PHIST_KERNEL_TEST_H

#include <iostream>
#include <iomanip>
#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#else
typedef int MPI_Comm;
#endif

#include "phist_typedefs.h"
#include "phist_kernels.h"

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
 int ierr_, mpi_rank_, mpi_size_;
 
 //! we store a pointer to the original stream buffer of std::cout,
 //! set it to NULL on rank!=0 at SetUp() and reset it at TearDown().
 std::streambuf *rdbuf_bak,*e_rdbuf_bak;
 
 std::ostream *cout;//! std::cout for everyone (std::cout is muted if rank!=0)
 std::ostream *cerr;//! std::cerr for everyone
 
	/** Set up method.
	 */
	virtual void SetUp() 
	{
#ifdef PHIST_HAVE_MPI	
	mpi_comm_=MPI_COMM_WORLD;
	phist_comm_create(&comm_,&ierr_);
	ASSERT_EQ(0,ierr_);
        ierr_=MPI_Comm_rank(mpi_comm_,&mpi_rank_);
	ASSERT_EQ(0,ierr_);
        ierr_=MPI_Comm_size(mpi_comm_,&mpi_size_);
	ASSERT_EQ(0,ierr_);
#else
        mpi_comm_=-1;
        mpi_rank_=0;
        mpi_size_=1;
#endif
#ifdef PHIST_HAVE_SP
	phist_Stype_avail(&ierr_); haveS_=(ierr_==0);
	phist_Ctype_avail(&ierr_); haveC_=(ierr_==0);
#else
        haveS_=false;
        haveC_=false;
#endif
	phist_Dtype_avail(&ierr_); haveD_=(ierr_==0);
	phist_Ztype_avail(&ierr_); haveZ_=(ierr_==0);
#if 0	
	rdbuf_bak = std::cout.rdbuf();
	e_rdbuf_bak = std::cerr.rdbuf();
	cout=new std::ostream(rdbuf_bak);
	cerr=new std::ostream(e_rdbuf_bak);
	if (mpi_rank_!=0) // this doesn't really seem to work as expected
	  {
	  std::cout.rdbuf(NULL);
	  std::cerr.rdbuf(NULL);
	  }
#else
        cout = &std::cout;
        cerr = &std::cerr;
#endif	
	// initialize random number sequence in a reproducible way (yet
	// with a different seed on each MPI process)
	rseed_ = (unsigned int)(mpi_rank_*77+42);
	srand(rseed_);
	}

virtual void TearDown()
  {
  if (false) // we do not delete the comm because it may be shared between
             // base classes of a derived class, so it is not clear when to
             // delete it without a smart pointer concept.
    {
    phist_comm_delete(comm_,&ierr_);
    ASSERT_EQ(0,ierr_);
    comm_=NULL;
    }
#if 0
	std::cout.rdbuf(rdbuf_bak);
	std::cerr.rdbuf(e_rdbuf_bak);
	delete cout;
	delete cerr;
#endif
  }
  
::testing::AssertionResult AssertNotNull(void* ptr)
  {
  if (ptr==NULL) return ::testing::AssertionFailure();
  return ::testing::AssertionSuccess();
  }

};

#endif
