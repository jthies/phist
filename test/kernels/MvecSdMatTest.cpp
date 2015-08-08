#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#include "phist_kernels.h"
#include "../kernels/KernelTest.h"
#include "../kernels/KernelTestWithMap.h"
#include "../kernels/KernelTestWithType.h"
#include "../kernels/KernelTestWithVectors.h"
#include "../kernels/KernelTestWithSdMats.h"


using namespace testing;

#define _BASESUITE_ kernels
#define _BASENAME_ MvecSdMatTest


#define _N_ 10
#define _M_ 1
#define _K_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 10
#define _M_ 1
#define _K_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 10
#define _M_ 4
#define _K_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 10
#define _M_ 1
#define _K_ 5
#include "../phist_typed_test_gen.h"

#define _N_ 64
#define _M_ 4
#define _K_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 64
#define _M_ 2
#define _K_ 2
#include "../phist_typed_test_gen.h"

#define _N_ 64
#define _M_ 1
#define _K_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 64
#define _M_ 4
#define _K_ 2
#include "../phist_typed_test_gen.h"

#define _N_ 64
#define _M_ 2
#define _K_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 64
#define _M_ 2
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
#define _M_ 12
#define _K_ 5
#include "../phist_typed_test_gen.h"

// case with square sdMats (m=k=5)
#define _N_ 64
#define _M_ 5
#define _K_ 5
#include "../phist_typed_test_gen.h"

// case with 7, 11
#define _N_ 512
#define _M_ 7
#define _K_ 11
#include "../phist_typed_test_gen.h"

// large tests for verifying high precision calculations
#ifdef PHIST_HIGH_PRECISION_KERNELS

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

