#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/phist_gtest.h"
//#include "gmock/gmock.h"

#include "phist_kernels.h"
#include "../kernels/KernelTest.h"
#include "../kernels/KernelTestWithMap.h"
#include "../kernels/KernelTestWithType.h"
#include "../kernels/KernelTestWithVectors.h"
#include "../kernels/KernelTestWithSdMats.h"


using namespace testing;

#define _BASESUITE_ kernels
#define CLASSFILE_DEF "MvecSdMatTest_def.hpp"
#define _BASENAME_ MvecSdMatTestWithAlignedViews
#define _USE_VIEWS_V1_ 2
#define _USE_VIEWS_V2_ 2
#define _USE_VIEWS_M_ 2


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

#define _N_ 25
#define _M_ 12
#define _K_ 12
#include "../phist_typed_test_gen.h"

#define _N_ 25
#define _M_ 1
#define _K_ 12
#include "../phist_typed_test_gen.h"

#define _N_ 25
#define _M_ 12
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