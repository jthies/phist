#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#include "kernels/phist_kernels.h"
#include "KernelTestWithSdMats.h"

using namespace testing;

#define _BASENAME_ SdMatTest

#define _N_ 5
#define _M_ 5
#include "../phist_typed_test_gen.h"

#define _N_ 6
#define _M_ 8
#include "../phist_typed_test_gen.h"

#define _N_ 8
#define _M_ 8
#include "../phist_typed_test_gen.h"

#define _N_ 25
#define _M_ 12
#include "../phist_typed_test_gen.h"

#define _N_ 12
#define _M_ 25
#include "../phist_typed_test_gen.h"
