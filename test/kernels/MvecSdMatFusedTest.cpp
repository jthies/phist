#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#if !defined(PHIST_HIGH_PRECISION_KERNELS) && defined(PHIST_HIGH_PRECISION_KERNELS_FORCE)
#define PHIST_HIGH_PRECISION_KERNELS
#endif

#include "gtest/gtest.h"


#include "phist_kernels.h"
#include "../kernels/KernelTestWithVectors.h"
#include "../kernels/KernelTestWithSdMats.h"


using namespace testing;

#define _BASENAME_ MvecSdMatFusedTest
#define CLASSFILE_DEF "MvecSdMatFusedTest_def.hpp"
#define _USE_VIEWS_V_ 0
#define _USE_VIEWS_W_ 0
#define _USE_VIEWS_M_ 0
#define _USE_VIEWS_N_ 0


#define _N_ 10
#define _M_ 1
#define _K_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 64
#define _M_ 12
#define _K_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 10
#define _M_ 1
#define _K_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 64
#define _M_ 12
#define _K_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 64
#define _M_ 4
#define _K_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 64
#define _M_ 1
#define _K_ 2
#include "../phist_typed_test_gen.h"

#define _N_ 10
#define _M_ 1
#define _K_ 5
#include "../phist_typed_test_gen.h"

#define _N_ 64
#define _M_ 12
#define _K_ 5
#include "../phist_typed_test_gen.h"

#define _N_ 64
#define _M_ 5
#define _K_ 5
#include "../phist_typed_test_gen.h"

#define _N_ 512
#define _M_ 7
#define _K_ 11
#include "../phist_typed_test_gen.h"

// large high precision tests (too slow in debug mode)
#if defined(PHIST_HIGH_PRECISION_KERNELS) && !defined(TESTING)
#define _N_ 4000000
#define _M_ 1
#define _K_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 4000000
#define _M_ 1
#define _K_ 2
#include "../phist_typed_test_gen.h"

#define _N_ 4000000
#define _M_ 1
#define _K_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 4000000
#define _M_ 2
#define _K_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 4000000
#define _M_ 2
#define _K_ 2
#include "../phist_typed_test_gen.h"

#define _N_ 4000000
#define _M_ 2
#define _K_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 4000000
#define _M_ 4
#define _K_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 4000000
#define _M_ 4
#define _K_ 2
#include "../phist_typed_test_gen.h"

#define _N_ 4000000
#define _M_ 4
#define _K_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 4000000
#define _M_ 10
#define _K_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 4000000
#define _M_ 10
#define _K_ 2
#include "../phist_typed_test_gen.h"

#define _N_ 4000000
#define _M_ 10
#define _K_ 4
#include "../phist_typed_test_gen.h"
#endif /* PHIST_HIGH_PRECISION_KERNELS */

