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

#include "KernelTestWithSparseMat.h"
#include "KernelTestWithVectors.h"
#include "KernelTestWithSdMats.h"

using namespace ::testing;

//////////////////////////////////////////
// symmetric matrix tests               //
//////////////////////////////////////////
#define _BASENAME_ SymmMatTestWithAlignedViews
#define CLASSFILE_DEF "SymmMatTest_def.hpp"
#define _USE_VIEWS_ 2

#define MATNAME MATNAME_symmMat

#define _N_ 20
#define _M_ 8
#include "../phist_typed_test_gen.h"

#undef MATNAME
#define MATNAME MATNAME_sprandsym

#define _N_ 163
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 163
#define _M_ 2
#include "../phist_typed_test_gen.h"

#define _N_ 163
#define _M_ 3
#include "../phist_typed_test_gen.h"

#define _N_ 163
#define _M_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 163
#define _M_ 7
#include "../phist_typed_test_gen.h"

#define _N_ 163
#define _M_ 30
#include "../phist_typed_test_gen.h"
