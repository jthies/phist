#include "phist_config.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"


#include "phist_macros.h"
#include "phist_kernels.h"
#include "phist_enums.h"
#include "phist_jadaOp.hpp"
#include "phist_jadaCorrectionSolver.h"
#include "phist_ScalarTraits.hpp"
#include "phist_orthog.h"
#include "../kernels/KernelTestWithSparseMat.h"
#include "../kernels/KernelTestWithSdMats.h"
#include "../kernels/KernelTestWithVectors.h"
#include "JadaTestWithOpts.hpp"

using namespace testing;

#define CLASSFILE_DEF "TestJadaCorrectionSolver_def.hpp"

#define MATNAME MATNAME_sprandn
#define _BASENAME_ TestJadaCorrectionSolver_sprandn
#define _MAXBAS_ 5

#define _N_ 25
#define _M_ 10
#define _K_ 2
#include "../phist_typed_test_gen.h"

#undef _MAXBAS_
#undef MATNAME
#undef _BASENAME_

#define MATNAME MATNAME_jadaTestMat
#define _BASENAME_ TestJadaCorrectionSolver_jadaTestMat
#define _MAXBAS_ 10

#define _N_ 512
#define _M_ 11
#define _K_ 7
#include "../phist_typed_test_gen.h"

