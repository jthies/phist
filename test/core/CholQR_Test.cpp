/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_config.h"

#include "phist_tools.h"
#include "phist_kernels.h"

#include "gtest/phist_gtest.h"

#include "../kernels/KernelTestWithVectors.h"
#include "../kernels/KernelTestWithSdMats.h"

#include "phist_chol_QR.h"


using namespace testing;

// we use the same _def file as the QR_Test for the kernel
// function mvec_QR to test the core function chol_QR
#define TEST_CHOL_QR
#define CLASSFILE_DEF "../kernels/QR_Test_def.hpp"

#define _BASENAME_ CholQR_Test

#define MVECS_VIEWED 0
#define SDMATS_VIEWED 0

#define _N_ 24
#define _M_ 1
#define CLASSFILE_DEF "../kernels/QR_Test_def.hpp"
#include "../phist_typed_test_gen.h"

#define _N_ 59
#define _M_ 5
#define CLASSFILE_DEF "../kernels/QR_Test_def.hpp"
#include "../phist_typed_test_gen.h"

// let's try something bigger...
#ifndef PHIST_TESTING
#define _N_ 9999
#define _M_ 65
#define CLASSFILE_DEF "../kernels/QR_Test_def.hpp"
#include "../phist_typed_test_gen.h"
#endif

/* small test but with viewed mvecs and sdMats */
#undef MVECS_VIEWED
#undef SDMATS_VIEWED
#define MVECS_VIEWED 1
#define SDMATS_VIEWED 1

#undef _BASENAME_
#define _BASENAME_ CholQR_TestWithUnalignedViews

#define _N_ 111
#define _M_ 8
#define CLASSFILE_DEF "../kernels/QR_Test_def.hpp"
#include "../phist_typed_test_gen.h"

#undef MVECS_VIEWED
#undef SDMATS_VIEWED
#define MVECS_VIEWED 2
#define SDMATS_VIEWED 2

#undef _BASENAME_
#define _BASENAME_ CholQR_TestWithAlignedViews

#define _N_ 111
#define _M_ 8
#define CLASSFILE_DEF "../kernels/QR_Test_def.hpp"
#include "../phist_typed_test_gen.h"
