#include "phist_config.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#include "phist_kernels.h"
#include "KernelTestWithMap.h"
#include "KernelTestWithVectors.h"

using namespace testing;

#define _BASENAME_ MvecTest
#define CLASSFILE_DEF "MvecTest_def.hpp"

#define _USE_VIEWS_ 0

// define MvecInitializer functions
#define FIRST_INSTANCE
#define _N_ 9
#define _M_ 1
#include "../phist_typed_test_gen.h"
#undef FIRST_INSTANCE

#define _N_ 16
#define _M_ 9
#include "../phist_typed_test_gen.h"

#define _N_ 25
#define _M_ 12
#include "../phist_typed_test_gen.h"

#define _N_ 100
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 100
#define _M_ 2
#include "../phist_typed_test_gen.h"

#define _N_ 100
#define _M_ 3
#include "../phist_typed_test_gen.h"

#define _N_ 100
#define _M_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 100
#define _M_ 5
#include "../phist_typed_test_gen.h"

#define _N_ 100
#define _M_ 6
#include "../phist_typed_test_gen.h"

#define _N_ 100
#define _M_ 7
#include "../phist_typed_test_gen.h"

#define _N_ 100
#define _M_ 8
#include "../phist_typed_test_gen.h"

#define _N_ 100
#define _M_ 9
#include "../phist_typed_test_gen.h"

#define _N_ 100
#define _M_ 10
#include "../phist_typed_test_gen.h"

#define _N_ 100
#define _M_ 17
#include "../phist_typed_test_gen.h"



// large high precision tests (too slow in debug mode)
#if defined(PHIST_HIGH_PRECISION_KERNELS) && !defined(TESTING)

#define _N_ 4000000
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 4000000
#define _M_ 2
#include "../phist_typed_test_gen.h"

#define _N_ 4000000
#define _M_ 4
#include "../phist_typed_test_gen.h"

#endif /* PHIST_HIGH_PRECISION_KERNELS */


// the Belos Tester selects the number of vectors itself,
// it is enough to test it for one vector length, I think
#ifdef PHIST_KERNEL_LIB_GHOST
#ifdef PHIST_HAVE_BELOS
// ghost views of mvecs are not general enough to pass these
// tests for row-major ordering
#ifndef PHIST_MVECS_ROW_MAJOR
#define DO_BELOS_TESTS
#include "phist_GhostMV.hpp"
#include "phist_rcp_helpers.hpp"
#include "Belos_GhostAdapter.hpp"
#include "BelosMVOPTester.hpp"
#include "BelosOutputManager.hpp"
#endif
#endif
#endif


#define _N_ 237
#define _M_ 4
#include "../phist_typed_test_gen.h"

// some tests on views of mvecs

#undef _BASENAME_
#define _BASENAME_ MvecTestWithUnalignedViews
#undef _USE_VIEWS_
#define _USE_VIEWS_ 1

#define _N_ 237
#define _M_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 237
#define _M_ 9
#include "../phist_typed_test_gen.h"

#undef _BASENAME_
#define _BASENAME_ MvecTestWithAlignedViews

#undef _USE_VIEWS_
#define _USE_VIEWS_ 2

#define _N_ 237
#define _M_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 237
#define _M_ 9
#include "../phist_typed_test_gen.h"

#ifdef DO_BELOS_TESTS
#undef DO_BELOS_TESTS
#endif

