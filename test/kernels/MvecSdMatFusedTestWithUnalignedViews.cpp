#include "phist_config.h"

#include "phist_tools.h"
#include "phist_kernels.h"

#include "gtest/phist_gtest.h"

#include "../kernels/KernelTestWithVectors.h"
#include "../kernels/KernelTestWithSdMats.h"


using namespace testing;

#define _BASENAME_ MvecSdMatFusedTestWithUnalignedViews
#define CLASSFILE_DEF "MvecSdMatFusedTest_def.hpp"
#define _USE_VIEWS_V_ 1
#define _USE_VIEWS_W_ 1
#define _USE_VIEWS_M_ 1
#define _USE_VIEWS_N_ 1


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

