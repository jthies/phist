#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"


#include "phist_macros.h"
#include "phist_kernels.h"
#include "phist_enums.h"
#include "phist_blockedgmres.h"
#include "phist_blockedminres.h"
#include "phist_ScalarTraits.hpp"
#include "../kernels/KernelTestWithSparseMat.h"
#include "../kernels/KernelTestWithVectors.h"

using namespace testing;

#define CLASSFILE_DEF "TestBlockedGMRES_def.hpp"

#define _BASENAME_ TestBlockedGMRES
#define MATNAME MATNAME_jadaTestMat
#define MAXBAS 25
#define TOLA 2.5e-2
#define TOLB 1.0e-4

#define _N_ 512
#define _M_ 6
#include "../phist_typed_test_gen.h"

#undef MATNAME
#undef MAXBAS
#undef TOLA
#undef TOLB
#undef _BASENAME_



#define _BASENAME_ TestBlockedGMRES_symmMat
#define MATNAME MATNAME_symmMat
#define MAXBAS 12
#define TOLA 2.5e-2
#define TOLB 1.0e-3
#define MATSYMMETRIC

#define _N_ 20
#define _M_ 5
#include "../phist_typed_test_gen.h"
