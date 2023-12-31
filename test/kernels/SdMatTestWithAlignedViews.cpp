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
#include "KernelTestWithSdMats.h"

using namespace testing;

#define _BASENAME_ SdMatTestWithAlignedViews
#define CLASSFILE_DEF "SdMatTest_def.hpp"
#define _USE_VIEWS_ 2

#define _N_ 1
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 2
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 1
#define _M_ 2
#include "../phist_typed_test_gen.h"

#define _N_ 3
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 1
#define _M_ 3
#include "../phist_typed_test_gen.h"

#define _N_ 2
#define _M_ 2
#include "../phist_typed_test_gen.h"

#define _N_ 3
#define _M_ 2
#include "../phist_typed_test_gen.h"

#define _N_ 2
#define _M_ 3
#include "../phist_typed_test_gen.h"

#define _N_ 5
#define _M_ 5
#include "../phist_typed_test_gen.h"

#define _N_ 6
#define _M_ 8
#include "../phist_typed_test_gen.h"

#define _N_ 8
#define _M_ 8
#include "../phist_typed_test_gen.h"

#define _N_ 25
#define _M_ 12
#include "../phist_typed_test_gen.h"

#define _N_ 12
#define _M_ 25
#include "../phist_typed_test_gen.h"

#define _N_ 20
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 1
#define _M_ 20
#include "../phist_typed_test_gen.h"
