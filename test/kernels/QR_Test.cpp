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

#include "KernelTestWithVectors.h"
#include "KernelTestWithSdMats.h"


using namespace testing;

// the _def file can also be used for testing chol_QR,
// but in the kernel tests we want to test the kernel 
// function mvec_QR (if implemented)
#define TEST_MVEC_QR

#define _BASENAME_ QR_Test

#define MVECS_VIEWED 0
#define SDMATS_VIEWED 0

#define _N_ 24
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 59
#define _M_ 5
#include "../phist_typed_test_gen.h"

// let's try something bigger... (slow in debug mode!)
#ifndef PHIST_TESTING
#define _N_ 9999
#define _M_ 65
#include "../phist_typed_test_gen.h"
#endif

/* small test but with viewed mvecs and sdMats */
#undef MVECS_VIEWED
#undef SDMATS_VIEWED
#define MVECS_VIEWED 1
#define SDMATS_VIEWED 1

#undef _BASENAME_
#define _BASENAME_ QR_TestWithUnalignedViews

#define _N_ 111
#define _M_ 8
#define CLASSFILE_DEF "QR_Test_def.hpp"
#include "../phist_typed_test_gen.h"

#undef MVECS_VIEWED
#undef SDMATS_VIEWED
#define MVECS_VIEWED 2
#define SDMATS_VIEWED 2

#undef _BASENAME_
#define _BASENAME_ QR_TestWithAlignedViews

#define _N_ 111
#define _M_ 8
#define CLASSFILE_DEF "QR_Test_def.hpp"
#include "../phist_typed_test_gen.h"
