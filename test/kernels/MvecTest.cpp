#include "phist_config.h"

#if !defined(PHIST_HIGH_PRECISION_KERNELS) && defined(PHIST_HIGH_PRECISION_KERNELS_FORCE)
#define PHIST_HIGH_PRECISION_KERNELS
#endif
#include "phist_tools.h"
#include "phist_kernels.h"

#include "gtest/phist_gtest.h"

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
#if defined(PHIST_HIGH_PRECISION_KERNELS) && !defined(PHIST_TESTING)

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
# ifdef PHIST_HAVE_BELOS
// ghost views of mvecs are not general enough to pass these
// tests for row-major ordering, however, in our own adapted
// version of the MVOPTester, we use only ascending col in- 
// dices, which should work with scattered views in GHOST.
//#  ifndef PHIST_MVECS_ROW_MAJOR
# ifndef GHOST_HAVE_CUDA
// these tests are a bit too hard for GHOST/CUDA up to now,
// even with our adjustments. We disable them for now because
// scattered views are not important for PHIST
#  define DO_BELOS_TESTS
# endif
#  include "phist_GhostMV.hpp"
#  include "phist_rcp_helpers.hpp"
#  include "Belos_GhostAdapter.hpp"
#  include "./BelosMVOPTester.hpp"
#  include "BelosOutputManager.hpp"
//      #  endif
# endif
#endif


#define _N_ 237
#define _M_ 4
#include "../phist_typed_test_gen.h"
