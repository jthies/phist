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

#include "matfuncs/matfuncs.h"

using namespace ::testing;

#define _USE_VIEWS_ 0
#define CLASSFILE_DEF "SparseMatTest_def.hpp"


#define MATNAME MATNAME_speye
#define _BASENAME_ SparseMatTest_speye
// included source code will define row functions,
// after we had all data types (S,D,C,Z) we undef it
#define FIRST_TIME

#define _N_ 25
#define _M_ 1
#include "../phist_typed_test_gen.h"

#undef FIRST_TIME

#define _N_ 25
#define _M_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 25
#define _M_ 7
#include "../phist_typed_test_gen.h"

#define _N_ 25
#define _M_ 12
#include "../phist_typed_test_gen.h"

// larger tests that also make sense on hybrid CPU/GPU nodes

#define _N_ 388
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 388
#define _M_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 388
#define _M_ 7
#include "../phist_typed_test_gen.h"

#define _N_ 388
#define _M_ 12
#include "../phist_typed_test_gen.h"

#undef MATNAME
#undef _BASENAME_


#define MATNAME MATNAME_spzero
// the case A=0 does not work
// for ghost because the binCRS format
// doesn't handle it correctly, it seems
#ifdef PHIST_KERNEL_LIB_GHOST
#define DISABLE_TESTCASE
#endif
#define _BASENAME_ SparseMatTest_spzero

#define _N_ 25
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 25
#define _M_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 25
#define _M_ 7
#include "../phist_typed_test_gen.h"

#define _N_ 25
#define _M_ 12
#include "../phist_typed_test_gen.h"

// larger tests that also make sense on hybrid CPU/GPU nodes

#define _N_ 388
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 388
#define _M_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 388
#define _M_ 7
#include "../phist_typed_test_gen.h"

#define _N_ 388
#define _M_ 12
#include "../phist_typed_test_gen.h"

#ifdef DISABLE_TESTCASE
#undef DISABLE_TESTCASE
#endif
#undef MATNAME
#undef _BASENAME_


#define MATNAME MATNAME_sprandn
#define _BASENAME_ SparseMatTest_sprandn

#define _N_ 25
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 25
#define _M_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 25
#define _M_ 7
#include "../phist_typed_test_gen.h"

#define _N_ 25
#define _M_ 12
#include "../phist_typed_test_gen.h"

// larger tests that also make sense on hybrid CPU/GPU nodes

#define _N_ 388
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 388
#define _M_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 388
#define _M_ 7
#include "../phist_typed_test_gen.h"

#define _N_ 388
#define _M_ 12
#include "../phist_typed_test_gen.h"

#undef MATNAME
#undef _BASENAME_


#define MATNAME MATNAME_sprandn_nodiag
#ifdef PHIST_KERNEL_LIB_GHOST
#define DISABLE_TESTCASE
#endif
#define _BASENAME_ SparseMatTest_sprandn_nodiag

#define _N_ 25
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 25
#define _M_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 25
#define _M_ 7
#include "../phist_typed_test_gen.h"

#define _N_ 25
#define _M_ 12
#include "../phist_typed_test_gen.h"

// larger tests that also make sense on hybrid CPU/GPU nodes

#define _N_ 388
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 388
#define _M_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 388
#define _M_ 7
#include "../phist_typed_test_gen.h"

#define _N_ 388
#define _M_ 12
#include "../phist_typed_test_gen.h"

#ifdef DISABLE_TESTCASE
#undef DISABLE_TESTCASE
#endif
#undef MATNAME
#undef _BASENAME_


#define MATNAME MATNAME_spshift
#ifdef PHIST_KERNEL_LIB_GHOST
#define DISABLE_TESTCASE
#endif
#define _BASENAME_ SparseMatTest_spshift

#define _N_ 25
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 25
#define _M_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 25
#define _M_ 7
#include "../phist_typed_test_gen.h"

#define _N_ 25
#define _M_ 12
#include "../phist_typed_test_gen.h"

// larger tests that also make sense on hybrid CPU/GPU nodes

#define _N_ 388
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 388
#define _M_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 388
#define _M_ 7
#include "../phist_typed_test_gen.h"

#define _N_ 388
#define _M_ 12
#include "../phist_typed_test_gen.h"

#ifdef DISABLE_TESTCASE
#undef DISABLE_TESTCASE
#endif

#undef MATNAME
#undef _BASENAME_

#define MATNAME MATNAME_BENCH3D_8_A1
#define _BASENAME_ SparseMatTest_BENCH3D_8_A1

#define _N_ 512
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 512
#define _M_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 512
#define _M_ 7
#include "../phist_typed_test_gen.h"

#define _N_ 512
#define _M_ 12
#include "../phist_typed_test_gen.h"
