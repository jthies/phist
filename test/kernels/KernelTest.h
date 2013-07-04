#ifndef PHIST_KERNEL_TEST_H
#define PHIST_KERNEL_TEST_H

#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#include "phist_typedefs.h"
#include "phist_kernels.h"

/** Base for calculation tests.
 * This class provides a basic test fixture for all calculation tests. It mainly provides a vector
 * with double precision values 1.0, 2.0 and 3.0 for usage as input to tests.
 * @note Derived classes need to call CalculationTest::SetUp() if they override it.
 */
class KernelTest: public testing::Test {
public:

 comm_ptr_t comm_;
 bool haveS_, haveD_, haveC_, haveZ_;
 int ierr_;
 
	/** Set up method.
	 * Fills internal data vector with values 1.0, 2.0 and 3.0.
	 */
	virtual void SetUp() 
	{
	phist_comm_create(&comm_,&ierr_);
	ASSERT_EQ(0,ierr_);
	phist_Stype_avail(&ierr_); haveS_=(ierr_==0);
	phist_Dtype_avail(&ierr_); haveD_=(ierr_==0);
	phist_Ctype_avail(&ierr_); haveC_=(ierr_==0);
	phist_Ztype_avail(&ierr_); haveZ_=(ierr_==0);
	}

::testing::AssertionResult AssertNotNull(void* ptr)
  {
  if (ptr==NULL) return ::testing::AssertionFailure();
  return ::testing::AssertionSuccess();
  }

};

/** This class is a mock object to test kernel implementation. */
class MockKernels {

public:
//    MOCK_CONST_METHOD0(evaluate, double());
};

#endif
