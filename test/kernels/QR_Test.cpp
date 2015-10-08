#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#ifdef PHIST_KERNEL_LIB_GHOST
#if defined(PHIST_HAVE_TEUCHOS)&&defined(PHIST_HAVE_KOKKOS)
#define HAVE_MVEC_QR
#endif
#else
#define HAVE_MVEC_QR
#endif

#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#include "KernelTestWithVectors.h"
#include "KernelTestWithSdMats.h"

#include "phist_kernels.h"


using namespace testing;

#define _BASENAME_ QR_Test

#define MVECS_VIEWED false
#define SDMATS_VIEWED false

#define _N_ 24
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 59
#define _M_ 5
#include "../phist_typed_test_gen.h"


/* same test but with viewed mvecs and sdMats */
#undef MVECS_VIEWED
#undef SDMATS_VIEWED
#define MVECS_VIEWED true
#define SDMATS_VIEWED true

#define _N_ 59
#define _M_ 7
#include "../phist_typed_test_gen.h"

#define _N_ 25
#define _M_ 12
#include "../phist_typed_test_gen.h"

#undef MVECS_VIEWED
#undef SDMATS_VIEWED
#define MVECS_VIEWED false
#define SDMATS_VIEWED false

// let's try something bigger...
#define _N_ 9999
#define _M_ 65
#include "../phist_typed_test_gen.h"

