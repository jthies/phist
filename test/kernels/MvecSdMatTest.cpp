/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_config.h"

#if !defined(PHIST_HIGH_PRECISION_KERNELS) && defined(PHIST_HIGH_PRECISION_KERNELS_FORCE)
#define PHIST_HIGH_PRECISION_KERNELS
#endif
#include "phist_tools.h"
#include "phist_kernels.h"

#include "gtest/phist_gtest.h"

#include "../kernels/KernelTestWithVectors.h"
#include "../kernels/KernelTestWithSdMats.h"


using namespace testing;

#define _BASENAME_ MvecSdMatTest
#define CLASSFILE_DEF "MvecSdMatTest_def.hpp"
#define _USE_VIEWS_V1_ 0
#define _USE_VIEWS_V2_ 0
#define _USE_VIEWS_M_ 0


#define _N_ 10
#define _M_ 1
#define _K_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 10
#define _M_ 1
#define _K_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 10
#define _M_ 4
#define _K_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 10
#define _M_ 1
#define _K_ 5
#include "../phist_typed_test_gen.h"

#define _N_ 64
#define _M_ 4
#define _K_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 25
#define _M_ 12
#define _K_ 12
#include "../phist_typed_test_gen.h"

#define _N_ 25
#define _M_ 1
#define _K_ 12
#include "../phist_typed_test_gen.h"

#define _N_ 25
#define _M_ 12
#define _K_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 64
#define _M_ 2
#define _K_ 2
#include "../phist_typed_test_gen.h"

#define _N_ 64
#define _M_ 1
#define _K_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 64
#define _M_ 4
#define _K_ 2
#include "../phist_typed_test_gen.h"

#define _N_ 64
#define _M_ 2
#define _K_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 64
#define _M_ 2
#define _K_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 64
#define _M_ 12
#define _K_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 64
#define _M_ 4
#define _K_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 64
#define _M_ 12
#define _K_ 5
#include "../phist_typed_test_gen.h"

// case with square sdMats (m=k=5)
#define _N_ 64
#define _M_ 5
#define _K_ 5
#include "../phist_typed_test_gen.h"

// this size caused problems in another test class (SymmMatTest) with CUDA and -np 12
#define _N_ 163
#define _M_ 30
#define _K_ 30
#include "../phist_typed_test_gen.h"

// case with 7, 11
#define _N_ 512
#define _M_ 7
#define _K_ 11
#include "../phist_typed_test_gen.h"


// large high precision tests (too slow in debug mode)
#if defined(PHIST_HIGH_PRECISION_KERNELS) && !defined(PHIST_TESTING)

#define _N_ 4000000
#define _M_ 1
#define _K_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 4000000
#define _M_ 1
#define _K_ 2
#include "../phist_typed_test_gen.h"

#define _N_ 4000000
#define _M_ 1
#define _K_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 4000000
#define _M_ 2
#define _K_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 4000000
#define _M_ 2
#define _K_ 2
#include "../phist_typed_test_gen.h"

#define _N_ 4000000
#define _M_ 2
#define _K_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 4000000
#define _M_ 4
#define _K_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 4000000
#define _M_ 4
#define _K_ 2
#include "../phist_typed_test_gen.h"

#define _N_ 4000000
#define _M_ 4
#define _K_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 4000000
#define _M_ 10
#define _K_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 4000000
#define _M_ 10
#define _K_ 2
#include "../phist_typed_test_gen.h"

#define _N_ 4000000
#define _M_ 10
#define _K_ 4
#include "../phist_typed_test_gen.h"

#endif /* PHIST_HIGH_PRECISION_KERNELS */
