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

#include "phist_tools.h"
#include "phist_kernels.h"
#include "phist_core.h"

#include "phist_jadaOp.h"

#include "../kernels/KernelTestWithSparseMat.h"
#include "../kernels/KernelTestWithMassMat.h"
#include "../kernels/KernelTestWithVectors.h"
#include "../kernels/KernelTestWithSdMats.h"

using namespace testing;

#define CLASSFILE_DEF "JadaOpTest_def.hpp"

#define _BASENAME_ JadaOpTest_speye
#define MATNAME MATNAME_speye

#define _N_ 25
#define _M_ 4
#define _K_ 10
#include "../phist_typed_test_gen.h"

#undef _BASENAME_
#undef MATNAME

#define _BASENAME_ JadaOpTest_sprandn
#define MATNAME MATNAME_sprandn

#define _N_ 25
#define _M_ 4
#define _K_ 10
#include "../phist_typed_test_gen.h"
