/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_config.h"

#include "phist_tools.h"
#include "phist_kernels.h"

#include "gtest/phist_gtest.h"

#include "KernelTestWithMap.h"
#include "TestWithType.h"
#include "KernelTestWithVectors.h"
#include "../tools/MatrixIO.h"


using namespace ::testing;

#define _BASENAME_ MultipleMatrixTest
#define _USE_VIEWS_ 0
#define _N_ 125
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 125
#define _M_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 125
#define _M_ 7
#include "../phist_typed_test_gen.h"

