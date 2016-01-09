#include "phist_config.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"


#include "phist_kernels.h"
#include "KernelTestWithMap.h"
#include "KernelTestWithVectors.h"

using namespace testing;

#define _BASENAME_ MvecTestWithAlignedViews
#define CLASSFILE_DEF "MvecTest_def.hpp"

#define _USE_VIEWS_ 2

// define MvecInitializer functions

#define _N_ 9
#define _M_ 1
#include "../phist_typed_test_gen.h"


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
