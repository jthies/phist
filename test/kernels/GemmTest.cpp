#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#include "phist_kernels.h"
#include "KernelTestWithVectors.h"
#include "KernelTestWithSdMats.h"

using namespace testing;

#define _BASENAME_ GemmTest

#define _N_ 1
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 2
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 1
#define _M_ 2
#include "../phist_typed_test_gen.h"

#define _N_ 2
#define _M_ 2
#include "../phist_typed_test_gen.h"

#define _N_ 3
#define _M_ 2
#include "../phist_typed_test_gen.h"

#define _N_ 2
#define _M_ 3
#include "../phist_typed_test_gen.h"

#define _N_ 4
#define _M_ 3
#include "../phist_typed_test_gen.h"

#define _N_ 4
#define _M_ 2
#include "../phist_typed_test_gen.h"

#define _N_ 4
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 1
#define _M_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 2
#define _M_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 5
#define _M_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 4
#define _M_ 5
#include "../phist_typed_test_gen.h"

#define _N_ 12
#define _M_ 9
#include "../phist_typed_test_gen.h"

#define _N_ 13
#define _M_ 18
#include "../phist_typed_test_gen.h"

#define _N_ 1
#define _M_ 19
#include "../phist_typed_test_gen.h"

#define _N_ 4
#define _M_ 15
#include "../phist_typed_test_gen.h"

