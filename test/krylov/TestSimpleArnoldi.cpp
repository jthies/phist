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

#include "phist_kernel_flags.h"
#include "phist_simple_arnoldi.h"

#include "../kernels/KernelTestWithSparseMat.h"
#include "../kernels/KernelTestWithMassMat.h"
#include "../kernels/KernelTestWithSdMats.h"
#include "../kernels/KernelTestWithVectors.h"

using namespace testing;

#define CLASSFILE_DEF "TestSimpleArnoldi_def.hpp"

// ==================== spzero matrix ====================
#ifdef PHIST_KERNEL_LIB_GHOST
#define DISABLE_TESTCASE
#endif

#define MATNAME MATNAME_spzero
#define BLOCK_SIZE 0
#define _BASENAME_ TestSimpleArnoldi_spzero

#define _N_ 25
#define _M_ 12
#include "../phist_typed_test_gen.h"


#undef _BASENAME_
#undef BLOCK_SIZE
#define BLOCK_SIZE 1
#define _BASENAME_ TestSimpleBlock1Arnoldi_spzero

#define _N_ 25
#define _M_ 12
#include "../phist_typed_test_gen.h"


#undef _BASENAME_
#undef BLOCK_SIZE
#define BLOCK_SIZE 2
#ifdef PHIST_KERNEL_LIB_GHOST
#define DISABLE_TESTCASE
#endif
#define _BASENAME_ TestSimpleBlock2Arnoldi_spzero

#define _N_ 25
#define _M_ 12
#include "../phist_typed_test_gen.h"


#undef _BASENAME_
#undef BLOCK_SIZE
#define BLOCK_SIZE 3
#define _BASENAME_ TestSimpleBlock3Arnoldi_spzero

#define _N_ 25
#define _M_ 12
#include "../phist_typed_test_gen.h"


#undef _BASENAME_
#undef BLOCK_SIZE
#define BLOCK_SIZE 4
#define _BASENAME_ TestSimpleBlock4Arnoldi_spzero

#define _N_ 25
#define _M_ 12
#include "../phist_typed_test_gen.h"

#undef _BASENAME_
#undef BLOCK_SIZE
#undef MATNAME
#ifdef DISABLE_TESTCASE
#undef DISABLE_TESTCASE
#endif


// ==================== speye matrix ====================
#define MATNAME MATNAME_speye
#define BLOCK_SIZE 0
#define _BASENAME_ TestSimpleArnoldi_speye

#define _N_ 25
#define _M_ 12
#include "../phist_typed_test_gen.h"


#undef _BASENAME_
#undef BLOCK_SIZE
#define BLOCK_SIZE 1
#define _BASENAME_ TestSimpleBlock1Arnoldi_speye

#define _N_ 25
#define _M_ 12
#include "../phist_typed_test_gen.h"


#undef _BASENAME_
#undef BLOCK_SIZE
#define BLOCK_SIZE 2
#ifdef PHIST_KERNEL_LIB_GHOST
#define DISABLE_TESTCASE
#endif
#define _BASENAME_ TestSimpleBlock2Arnoldi_speye

#define _N_ 25
#define _M_ 12
#include "../phist_typed_test_gen.h"


#undef _BASENAME_
#undef BLOCK_SIZE
#define BLOCK_SIZE 3
#define _BASENAME_ TestSimpleBlock3Arnoldi_speye

#define _N_ 25
#define _M_ 12
#include "../phist_typed_test_gen.h"


#undef _BASENAME_
#undef BLOCK_SIZE
#define BLOCK_SIZE 4
#define _BASENAME_ TestSimpleBlock4Arnoldi_speye

#define _N_ 25
#define _M_ 12
#include "../phist_typed_test_gen.h"

#undef _BASENAME_
#undef BLOCK_SIZE
#undef MATNAME
#ifdef DISABLE_TESTCASE
#undef DISABLE_TESTCASE
#endif

// ==================== sprandn matrix ====================
#define MATNAME MATNAME_sprandn
#define BLOCK_SIZE 0
#define _BASENAME_ TestSimpleArnoldi_sprandn

#define _N_ 25
#define _M_ 12
#include "../phist_typed_test_gen.h"


#undef _BASENAME_
#undef BLOCK_SIZE
#define BLOCK_SIZE 1
#define _BASENAME_ TestSimpleBlock1Arnoldi_sprandn

#define _N_ 25
#define _M_ 12
#include "../phist_typed_test_gen.h"


#undef _BASENAME_
#undef BLOCK_SIZE
#define BLOCK_SIZE 2
#ifdef PHIST_KERNEL_LIB_GHOST
#define DISABLE_TESTCASE
#endif
#define _BASENAME_ TestSimpleBlock2Arnoldi_sprandn

#define _N_ 25
#define _M_ 12
#include "../phist_typed_test_gen.h"


#undef _BASENAME_
#undef BLOCK_SIZE
#define BLOCK_SIZE 3
#define _BASENAME_ TestSimpleBlock3Arnoldi_sprandn

#define _N_ 25
#define _M_ 12
#include "../phist_typed_test_gen.h"


#undef _BASENAME_
#undef BLOCK_SIZE
#define BLOCK_SIZE 4
#define _BASENAME_ TestSimpleBlock4Arnoldi_sprandn

#define _N_ 25
#define _M_ 12
#include "../phist_typed_test_gen.h"

#undef _BASENAME_
#undef BLOCK_SIZE
#undef MATNAME
#ifdef DISABLE_TESTCASE
#undef DISABLE_TESTCASE
#endif

// ==================== sprandn_nodiag matrix ====================
#ifdef PHIST_KERNEL_LIB_GHOST
#define DISABLE_TESTCASE
#endif

#define MATNAME MATNAME_sprandn_nodiag
#define BLOCK_SIZE 0
#define _BASENAME_ TestSimpleArnoldi_sprandn_nodiag

#define _N_ 25
#define _M_ 12
#include "../phist_typed_test_gen.h"


#undef _BASENAME_
#undef BLOCK_SIZE
#define BLOCK_SIZE 1
#define _BASENAME_ TestSimpleBlock1Arnoldi_sprandn_nodiag

#define _N_ 25
#define _M_ 12
#include "../phist_typed_test_gen.h"


#undef _BASENAME_
#undef BLOCK_SIZE
#define BLOCK_SIZE 2
#ifdef PHIST_KERNEL_LIB_GHOST
#define DISABLE_TESTCASE
#endif
#define _BASENAME_ TestSimpleBlock2Arnoldi_sprandn_nodiag

#define _N_ 25
#define _M_ 12
#include "../phist_typed_test_gen.h"


#undef _BASENAME_
#undef BLOCK_SIZE
#define BLOCK_SIZE 3
#define _BASENAME_ TestSimpleBlock3Arnoldi_sprandn_nodiag

#define _N_ 25
#define _M_ 12
#include "../phist_typed_test_gen.h"


#undef _BASENAME_
#undef BLOCK_SIZE
#define BLOCK_SIZE 4
#define _BASENAME_ TestSimpleBlock4Arnoldi_sprandn_nodiag

#define _N_ 25
#define _M_ 12
#include "../phist_typed_test_gen.h"


#undef _BASENAME_
#undef BLOCK_SIZE
#undef MATNAME
#ifdef DISABLE_TESTCASE
#undef DISABLE_TESTCASE
#endif

