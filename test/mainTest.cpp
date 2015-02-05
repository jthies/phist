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
#ifdef PHIST_KERNEL_LIB_GHOST
#include "ghost/phist_ghost_macros.hpp"
#endif

GTEST_API_ int main(int argc, char **argv) {
    int iflag,test_result;
    phist_kernels_init(&argc,&argv,&iflag);
    PHIST_GHOST_TASK_BEGIN
    testing::InitGoogleTest(&argc, argv);

    int rank = 0;
#ifdef PHIST_HAVE_MPI
    iflag = MPI_Comm_rank(MPI_COMM_WORLD,&rank);
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

    PHIST_GHOST_TASK_END
    phist_kernels_finalize(&iflag);
    //ASSERT_INT_EQ(iflag,0);
    return test_result;
}
