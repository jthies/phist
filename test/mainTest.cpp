/*
 * main_test.cpp
 *
 *  Created on: 09.02.2012
 *      Author: schlauch
 */

#include <iostream>

#include "gtest/gtest.h"
#include "kernels/phist_kernels.h"

GTEST_API_ int main(int argc, char **argv) {
    int ierr,test_result;
    phist_kernels_init(&argc,&argv,&ierr);
    //ASSERT_INT_EQ(ierr,0);
    std::cout << "Running main() from mainTest.cpp" << std::endl;

    testing::InitGoogleTest(&argc, argv);
    test_result=RUN_ALL_TESTS();
    phist_kernels_finalize(&ierr);
    //ASSERT_INT_EQ(ierr,0);
    return test_result;
    }
