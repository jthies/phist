#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/phist_gtest.h"


#include "kernels/phist_kernels.h"
#include "core/phist_sdFact.h"
#include "../kernels/KernelTestWithSdMats.h"

using namespace testing;

#define _BASENAME_ SdFactTestWithAlignedViews
#define CLASSFILE_DEF "SdFactTest_def.hpp"
#define _USE_VIEWS_ 2

#define _N_ 1
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 2
#define _M_ 2
#include "../phist_typed_test_gen.h"

#define _N_ 5
#define _M_ 5
#include "../phist_typed_test_gen.h"

#define _N_ 8
#define _M_ 8
#include "../phist_typed_test_gen.h"

