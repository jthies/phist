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

#include "phist_precon.h"
#include "phist_jadaOp.h"
#include "phist_jadaCorrectionSolver.h"
#include "phist_orthog.h"
#include "../kernels/KernelTestWithSparseMat.h"
#include "../kernels/KernelTestWithSdMats.h"
#include "../kernels/KernelTestWithVectors.h"
#include "../jada/JadaTestWithOpts.hpp"

using namespace testing;

#define CLASSFILE_DEF "TestPreconOperators_def.hpp"

#define MATNAME MATNAME_hpd_tridiag
#define PRECNAME MATNAME_hpd_tridiag_ainv
#define _BASENAME_ TestPreconOperators_hpd_tridiag
#define _MAXBAS_ 20

// matrix size
#define _N_ 64
// number of iteration vectors _NV_
#define _M_ 1
// number of projection vectors in q
#define _NVP_ 1
#include "../phist_typed_test_gen.h"

#undef MATNAME
#undef PRECNAME
#undef _BASENAME_

#define _N_ 64
#define _M_ 1
#define _NVP_ 1
#define MATNAME MATNAME_nhpd_tridiag
#define PRECNAME MATNAME_nhpd_tridiag_ainv
#define _BASENAME_ TestPreconOperators_nhpd_tridiag
#include "../phist_typed_test_gen.h"
