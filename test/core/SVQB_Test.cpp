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

#include "phist_svqb.h"

using namespace testing;

#define _BASENAME_ SVQB_Test

#define _N_ 24
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 59
#define _M_ 5
#include "../phist_typed_test_gen.h"

// let's try something bigger...
#ifndef PHIST_TESTING
#define _N_ 9999
#define _M_ 65
#include "../phist_typed_test_gen.h"
#endif
