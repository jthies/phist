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

#include "phist_operator.h"
#include "phist_jadaOp.h"
#include "phist_orthog.h"

#include "../kernels/KernelTestWithSparseMat.h"
#include "../kernels/KernelTestWithVectors.h"
#include "../kernels/KernelTestWithSdMats.h"

#ifdef PHIST_HAVE_BELOS
# ifdef PHIST_KERNEL_LIB_TPETRA__disabled
/* note: compiling these tests with Tpetra+CUDA requires the nvcc_wrapper script,
    but we try to hide all Tpetra/Trilinos headers from the user to avoid this.
    So we disable them completely for now.
*/
# include "phist_tpetra_typedefs.hpp"
# include "BelosTpetraAdapter.hpp"
# elif defined(PHIST_KERNEL_LIB_EPETRA)
# include "Epetra_MultiVector.h"
# include "BelosEpetraAdapter.hpp"
# elif defined(PHIST_KERNEL_LIB_GHOST)
# include "Belos_GhostAdapter.hpp"
# endif
#endif

#ifdef PHIST_HAVE_BELOS
#include "phist_rcp_helpers.hpp"
#include "phist_BelosOperatorTraits.hpp"
#include "BelosMVOPTester.hpp"
#endif

using namespace testing;
using namespace phist::testing;

#define CLASSFILE_DEF "OpTest_def.hpp"

#define _BASENAME_ OpTest_sprandn
#define MATNAME MATNAME_sprandn

#define _N_ 25
#define _M_ 8
#include "../phist_typed_test_gen.h"

#undef _BASENAME_
#undef MATNAME


#define _BASENAME_ OpTest_sprandn_nodiag
#ifdef PHIST_KERNEL_LIB_GHOST
#define DISABLE_TESTCASE
#endif
#define MATNAME MATNAME_sprandn_nodiag

#define _N_ 25
#define _M_ 8
#include "../phist_typed_test_gen.h"

#ifdef DISABLE_TESTCASE
#undef DISABLE_TESTCASE
#endif
#undef _BASENAME_
#undef MATNAME
