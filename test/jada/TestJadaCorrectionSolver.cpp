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

#include "phist_jadaOp.h"
#include "phist_jadaCorrectionSolver.h"
#include "phist_orthog.h"
#include "../kernels/KernelTestWithSparseMat.h"
#include "../kernels/KernelTestWithSdMats.h"
#include "../kernels/KernelTestWithVectors.h"
#include "JadaTestWithOpts.hpp"

using namespace testing;

#define CLASSFILE_DEF "TestJadaCorrectionSolver_def.hpp"

#define MATNAME MATNAME_sprandn
#define _BASENAME_ TestJadaCorrectionSolver_sprandn
#define _MAXBAS_ 5

#define _N_ 25
#define _M_ 10
#define _K_ 2
#include "../phist_typed_test_gen.h"

#undef _MAXBAS_
#undef MATNAME
#undef _BASENAME_

#define MATNAME MATNAME_jadaTestMat
#define _BASENAME_ TestJadaCorrectionSolver_jadaTestMat
#define _MAXBAS_ 10

#define _N_ 512
#define _M_ 11
#define _K_ 7
#include "../phist_typed_test_gen.h"

