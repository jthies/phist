#ifndef ESSEX_KERNEL_TEST_H
#define ESSEX_KERNEL_TEST_H

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "kernels/essex_kernels.h"

/** Base for calculation tests.
 * This class provides a basic test fixture for all calculation tests. It mainly provides a vector
 * with double precision values 1.0, 2.0 and 3.0 for usage as input to tests.
 * @note Derived classes need to call CalculationTest::SetUp() if they override it.
 */
class KernelTest: public testing::Test {
public:
	/** Set up method.
	 * Fills internal data vector with values 1.0, 2.0 and 3.0.
	 */
	virtual void SetUp() {
	}

 int ierr_;
};


/** This class is a mock object to test kernel implementation. */
class MockKernels {

public:
//    MOCK_CONST_METHOD0(evaluate, double());
};

#endif /* CALCULATIONTEST_H_ */
