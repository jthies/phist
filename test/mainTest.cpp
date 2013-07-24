/*
 * main_test.cpp
 *
 *  Created on: 09.02.2012
 *      Author: schlauch
 */

#include <iostream>

#include "gtest/gtest.h"
#include "kernels/phist_kernels.h"
#include "MpiRootOnlyPrinter.hpp"
#ifdef PHIST_HAVE_GHOST
#include "ghost.h"
#endif

GTEST_API_ int main(int argc, char **argv) {
    int ierr,test_result;
    testing::InitGoogleTest(&argc, argv);
    phist_kernels_init(&argc,&argv,&ierr);
#ifndef PHIST_KERNEL_LIB_GHOST
#ifdef PHIST_HAVE_GHOST
    ghost_init(argc,argv);
#endif
#endif
    // Gets hold of the event listener list.  
    ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();  
    // Adds a listener to the end.  Google Test takes the ownership.  
    listeners.Append(new MpiRootOnlyPrinter());

    test_result=RUN_ALL_TESTS();
    phist_kernels_finalize(&ierr);
#ifndef PHIST_KERNEL_LIB_GHOST
#ifdef PHIST_HAVE_GHOST
    ghost_finish();
#endif
#endif
    //ASSERT_INT_EQ(ierr,0);
    return test_result;
    }
