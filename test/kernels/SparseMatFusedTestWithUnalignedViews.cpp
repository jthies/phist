#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/phist_gtest.h"


#include "kernels/phist_kernels.h"
#include "KernelTestWithSparseMat.h"
#include "KernelTestWithVectors.h"
#include "KernelTestWithSdMats.h"

using namespace ::testing;

#define _USE_VIEWS_ 1
#define CLASSFILE_DEF "SparseMatFusedTest_def.hpp"


#define MATNAME MATNAME_speye
#define _BASENAME_ SparseMatFusedTestWithUnalignedViews_speye
// included source code will define row functions,
// after we had all data types (S,D,C,Z) we undef it
#define FIRST_TIME

#define _N_ 25
#define _M_ 1
#include "../phist_typed_test_gen.h"

#undef FIRST_TIME

#define _N_ 25
#define _M_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 25
#define _M_ 7
#include "../phist_typed_test_gen.h"

#define _N_ 25
#define _M_ 12
#include "../phist_typed_test_gen.h"

// larger tests that also make sense on hybrid CPU/GPU nodes

#define _N_ 388
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 388
#define _M_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 388
#define _M_ 7
#include "../phist_typed_test_gen.h"

#define _N_ 388
#define _M_ 12
#include "../phist_typed_test_gen.h"

#undef MATNAME
#undef _BASENAME_


#define MATNAME MATNAME_sprandn
#define _BASENAME_ SparseMatFusedTestWithUnalignedViews_sprandn

#define _N_ 25
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 25
#define _M_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 25
#define _M_ 7
#include "../phist_typed_test_gen.h"

#define _N_ 25
#define _M_ 12
#include "../phist_typed_test_gen.h"

// larger tests that also make sense on hybrid CPU/GPU nodes

#define _N_ 388
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 388
#define _M_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 388
#define _M_ 7
#include "../phist_typed_test_gen.h"

#define _N_ 388
#define _M_ 12
#include "../phist_typed_test_gen.h"

#undef MATNAME
#undef _BASENAME_

