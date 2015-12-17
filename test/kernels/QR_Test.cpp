#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"

#ifdef PHIST_KERNEL_LIB_GHOST
#include "ghost/config.h"
#if defined(PHIST_HAVE_TEUCHOS)&&defined(PHIST_HAVE_KOKKOS)&&(!defined(GHOST_HAVE_CUDA))
#define HAVE_MVEC_QR
#endif
#else
#define HAVE_MVEC_QR
#endif

#include "KernelTestWithVectors.h"
#include "KernelTestWithSdMats.h"

#include "phist_kernels.h"


using namespace testing;

#define _BASENAME_ QR_Test

#define MVECS_VIEWED 0
#define SDMATS_VIEWED 0

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

/* small test but with viewed mvecs and sdMats */
#undef MVECS_VIEWED
#undef SDMATS_VIEWED
#define MVECS_VIEWED 1
#define SDMATS_VIEWED 1

#undef _BASENAME_
#define _BASENAME_ QR_TestWithUnalignedViews

#define _N_ 111
#define _M_ 8
#define CLASSFILE_DEF "QR_Test_def.hpp"
#include "../phist_typed_test_gen.h"

#undef MVECS_VIEWED
#undef SDMATS_VIEWED
#define MVECS_VIEWED 2
#define SDMATS_VIEWED 2

#undef _BASENAME_
#define _BASENAME_ QR_TestWithAlignedViews

#define _N_ 111
#define _M_ 8
#define CLASSFILE_DEF "QR_Test_def.hpp"
#include "../phist_typed_test_gen.h"
