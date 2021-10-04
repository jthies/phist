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

#include "phist_subspacejada.h"
#include "phist_orthog.h"
#include "../kernels/KernelTestWithSparseMat.h"
#include "../kernels/KernelTestWithMassMat.h"
#include "../kernels/KernelTestWithSdMats.h"
#include "../kernels/KernelTestWithVectors.h"
#include "JadaTestWithOpts.hpp"

using namespace testing;

#define CLASSFILE_DEF "TestSubspaceJaDa_def.hpp"

#define MATNAME MATNAME_jadaTestMat
#define INNER_SOLVTYPE phist_GMRES
#define _BASENAME_ TestSubspaceJaDa_jadaTestMat_GMRES

// A and B are N x N matrices,
// we use a block size of K and
// aim to find numEigs eigenpairs, with _M_=numEigs+_K_-1
#define _N_ 512
#define _K_ 2
#define _M_ 7
#include "../phist_typed_test_gen.h"

