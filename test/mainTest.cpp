/*
 * main_test.cpp
 *
 *  Created on: 09.02.2012
 *      Author: schlauch
 */

#include "phist_config.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"
#include "kernels/phist_kernels.h"

GTEST_API_ int main(int argc, char **argv) {
    int ierr,test_result;
    testing::InitGoogleTest(&argc, argv);
    phist_kernels_init(&argc,&argv,&ierr);

    int rank = 0;
#ifdef PHIST_HAVE_MPI
    ierr = MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif

    // on MPI ranks != 0 remove the default output listeners if there are any
    if( rank != 0 )
    {
      ::testing::TestEventListener* defaultListener = ::testing::UnitTest::GetInstance()->listeners().default_result_printer();
      ::testing::UnitTest::GetInstance()->listeners().Release(defaultListener);
      delete defaultListener;

      ::testing::TestEventListener* defaultXMLListener = ::testing::UnitTest::GetInstance()->listeners().default_xml_generator();
      ::testing::UnitTest::GetInstance()->listeners().Release(defaultXMLListener);
      delete defaultXMLListener;
    }

    test_result=RUN_ALL_TESTS();

    phist_kernels_finalize(&ierr);
    //ASSERT_INT_EQ(ierr,0);
    return test_result;
}
