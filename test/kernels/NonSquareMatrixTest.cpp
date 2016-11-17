#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/phist_gtest.h"

#include "kernels/phist_kernels.h"
#include "KernelTestWithSparseMat.h"
#include "KernelTestWithVectors.h"
#include "KernelTestWithSdMats.h"

#include "matfuncs/matfuncs.h"

using namespace ::testing;

// so far most kernel libs do not support non-square matrices
#ifndef PHIST_KERNEL_LIB_EPETRA
#define DISABLE_TESTCASE
#endif

#define _USE_VIEWS_ 0
#define CLASSFILE_DEF "NonSquareMatrixTest_def.hpp"


#define MATNAME MATNAME_IDFUNC
#define _BASENAME_ NonSquareMatrixTest_speye_aug

#define _N_ 30
#define _M_ 25
#define _K_ 1
#include "../phist_typed_test_gen.h"


#define _N_ 30
#define _M_ 25
#define _K_ 4
#include "../phist_typed_test_gen.h"


#define _N_ 25
#define _M_ 30
#define _K_ 1
#include "../phist_typed_test_gen.h"


#define _N_ 25
#define _M_ 30
#define _K_ 4
#include "../phist_typed_test_gen.h"

