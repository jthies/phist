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

#define _BASENAME_ SparseMatTest
#define _USE_VIEWS_ 0
#define CLASSFILE_DEF "SparseMatTest_def.hpp"

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

#undef _USE_VIEWS_
#undef _BASENAME_


#define _BASENAME_ SparseMatTestWithUnalignedViews
#define _USE_VIEWS_ 1

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

#undef _USE_VIEWS_
#undef _BASENAME_


#define _BASENAME_ SparseMatTestWithAlignedViews
#define _USE_VIEWS_ 2

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


#undef _BASENAME_
#undef CLASSFILE_DEF
#undef _USE_VIEWS_

//////////////////////////////////////////
// symmetric matrix tests               //
//////////////////////////////////////////
#define _BASENAME_ SymmMatTest
#define CLASSFILE_DEF "SymmMatTest_def.hpp"
#define _USE_VIEWS_ 0

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

#undef _USE_VIEWS_
#undef _BASENAME_


#define _BASENAME_ SymmMatTestWithUnalignedViews
#define _USE_VIEWS_ 1

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

#undef _USE_VIEWS_
#undef _BASENAME_


#define _BASENAME_ SymmMatTestWithAlignedViews
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

