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

#include "core/phist_sdFact.h"
#include "../kernels/KernelTestWithSdMats.h"

using namespace testing;

#define _BASENAME_ SdFactTest
#define CLASSFILE_DEF "SdFactTest_def.hpp"
#define _USE_VIEWS_ 0

#define _N_ 1
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 2
#define _M_ 2
#include "../phist_typed_test_gen.h"

#define _N_ 5
#define _M_ 5
#include "../phist_typed_test_gen.h"

#define _N_ 8
#define _M_ 8
#include "../phist_typed_test_gen.h"


#define _N_ 11
#define _M_ 16
#include "../phist_typed_test_gen.h"

#define _N_ 16
#define _M_ 11
#include "../phist_typed_test_gen.h"

