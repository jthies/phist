#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#include "kernels/phist_kernels.h"
#include "KernelTestWithVectors.h"
#include "KernelTestWithSdMats.h"
#include "../tools/MatrixIO.h"

using namespace ::testing;

//////////////////////////////////////////
// symmetric matrix tests               //
//////////////////////////////////////////
#define _BASENAME_ SymmMatTestWithAlignedViews
#define CLASSFILE_DEF "SymmMatTest_def.hpp"
#define _USE_VIEWS_ 2

#define _N_ 20
#define _M_ 8
#include "../phist_typed_test_gen.h"

#define _N_ 163
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 163
#define _M_ 2
#include "../phist_typed_test_gen.h"

#define _N_ 163
#define _M_ 3
#include "../phist_typed_test_gen.h"

#define _N_ 163
#define _M_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 163
#define _M_ 7
#include "../phist_typed_test_gen.h"

#define _N_ 163
#define _M_ 30
#include "../phist_typed_test_gen.h"
