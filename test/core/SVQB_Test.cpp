#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif


#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#include "../kernels/KernelTestWithVectors.h"
#include "../kernels/KernelTestWithSdMats.h"

#include "phist_kernels.h"
#include "phist_svqb.h"

using namespace testing;

#define _BASENAME_ SVQB_Test

#define _N_ 24
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 59
#define _M_ 5
#include "../phist_typed_test_gen.h"

// let's try something bigger...
#define _N_ 9999
#define _M_ 65
#include "../phist_typed_test_gen.h"
