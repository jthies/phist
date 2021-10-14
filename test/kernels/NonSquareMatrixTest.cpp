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

#include "KernelTestWithSparseMat.h"
#include "KernelTestWithVectors.h"
#include "KernelTestWithSdMats.h"

#include "matfuncs/matfuncs.h"

using namespace ::testing;

// so far most kernel libs do not support non-square matrices
#ifndef PHIST_KERNEL_LIB_EPETRA
#define DISABLE_TESTCASE
#endif

#define _USE_VIEWS_ 0
#define CLASSFILE_DEF "NonSquareMatrixTest_def.hpp"


#define MATNAME MATNAME_IDFUNC
#define _BASENAME_ NonSquareMatrixTest_speye_aug

#define _N_ 30
#define _M_ 25
#define _K_ 1
#include "../phist_typed_test_gen.h"


#define _N_ 30
#define _M_ 25
#define _K_ 4
#include "../phist_typed_test_gen.h"


#define _N_ 25
#define _M_ 30
#define _K_ 1
#include "../phist_typed_test_gen.h"


#define _N_ 25
#define _M_ 30
#define _K_ 4
#include "../phist_typed_test_gen.h"

