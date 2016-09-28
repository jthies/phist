#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/phist_gtest.h"

#include "../tools/MatrixIO.h"

#include "phist_orthog.h"
#include "phist_kernels.h"

#include "../kernels/KernelTestWithMassMat.h"
#include "../kernels/KernelTestWithVectors.h"
#include "../kernels/KernelTestWithSdMats.h"


using namespace testing;

#define CLASSFILE_DEF "OrthogTest_def.hpp"
#define _BASENAME_ OrthogTestWithB

#define ORTHOG_WITH_HPD_B

#define _N_ 17
#define _M_ 1
#define _K_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 17
#define _M_ 6
#define _K_ 1
#include "../phist_typed_test_gen.h"

// larger block size for W
#define _N_ 64
#define _M_ 5
#define _K_ 3
#include "../phist_typed_test_gen.h"

