/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
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

#include "gtest/phist_gtest.h"
#include "tools/phist_tools.h"
#include "kernels/phist_kernels.h"

GTEST_API_ int main(int argc, char **argv) {
    int iflag,test_result;
    iflag=PHIST_IFLAG_DEFAULT;
    phist_kernels_init(&argc,&argv,&iflag);
    PHIST_MAIN_TASK_BEGIN
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

    PHIST_MAIN_TASK_END
    // do not print timing information
    iflag=PHIST_KERNELS_QUIET;
    phist_kernels_finalize(&iflag);
    //ASSERT_INT_EQ(iflag,0);
    return test_result;
}
