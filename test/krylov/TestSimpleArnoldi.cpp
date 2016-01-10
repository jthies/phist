#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"


#include "phist_macros.h"
#include "phist_kernels.h"
#include "phist_enums.h"
#include "phist_simple_arnoldi.h"
#include "../kernels/KernelTest.h"
#include "../kernels/KernelTestWithMap.h"
#include "../kernels/TestWithType.h"
#include "../kernels/KernelTestWithSdMats.h"
#include "../kernels/KernelTestWithVectors.h"

using namespace testing;

#ifdef PHIST_KERNEL_LIB_GHOST
// we can't easily store an all-zero matrix in ghost binCRS, so skip these tests
#define SKIP_ZERO_MAT
/* mvec_QR is not (yet) implemented in GHOST unless TSQR is available */
#if !defined(PHIST_HAVE_KOKKOS)||!defined(PHIST_HAVE_BELOS)
# define SKIP_BLOCK_ARNOLDI
#endif
#endif

#ifdef CLASSNAME
#undef CLASSNAME
#endif

#define _N_ 25
#define _M_ 12
#define BLOCK_SIZE 0

#ifdef PHIST_HAVE_SP

#define CLASSNAME STestSimpleArnoldi_25_12
#include "phist_gen_s.h"
#include "../tools/MatrixIO_def.hpp"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME

#define CLASSNAME CTestSimpleArnoldi_25_12
#include "phist_gen_c.h"
#include "../tools/MatrixIO_def.hpp"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME

#endif

#define CLASSNAME DTestSimpleArnoldi_25_12
#include "phist_gen_d.h"
#include "../tools/MatrixIO_def.hpp"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME

#define CLASSNAME ZTestSimpleArnoldi_25_12
#include "phist_gen_z.h"
#include "../tools/MatrixIO_def.hpp"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME


#undef BLOCK_SIZE
#define BLOCK_SIZE 1

#ifdef PHIST_HAVE_SP

#define CLASSNAME STestSimpleBlock1Arnoldi_25_12
#include "phist_gen_s.h"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME

#define CLASSNAME CTestSimpleBlock1Arnoldi_25_12
#include "phist_gen_c.h"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME

#endif

#define CLASSNAME DTestSimpleBlock1Arnoldi_25_12
#include "phist_gen_d.h"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME

#define CLASSNAME ZTestSimpleBlock1Arnoldi_25_12
#include "phist_gen_z.h"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME


#undef BLOCK_SIZE
#define BLOCK_SIZE 2

/* the tests below do not work unless mvec_QR is implemented! */
#ifdef SKIP_BLOCK_ARNOLDI
#define BLOCK_PREFIX(_T,_NAME,_BS) DISABLED_ ## _T ## Block ## _BS ## _NAME
#else
#define BLOCK_PREFIX(_T,_NAME,_BS) _T ## Block ## _BS ## _NAME
#endif

#ifdef PHIST_HAVE_SP

#define CLASSNAME BLOCK_PREFIX(S,Arnoldi_25_12,2)
#include "phist_gen_s.h"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME

#define CLASSNAME BLOCK_PREFIX(C,Arnoldi_25_12,2)
#include "phist_gen_c.h"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME

#endif

#define CLASSNAME BLOCK_PREFIX(D,Arnoldi_25_12,2)
#include "phist_gen_d.h"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME

#define CLASSNAME BLOCK_PREFIX(Z,Arnoldi_25_12,2)
#include "phist_gen_z.h"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME


#undef BLOCK_SIZE
#define BLOCK_SIZE 3

#ifdef PHIST_HAVE_SP

#define CLASSNAME BLOCK_PREFIX(S,Arnoldi_25_12,3)
#include "phist_gen_s.h"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME

#define CLASSNAME BLOCK_PREFIX(C,Arnoldi_25_12,3)
#include "phist_gen_c.h"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME

#endif

#define CLASSNAME BLOCK_PREFIX(D,Arnoldi_25_12,3)
#include "phist_gen_d.h"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME

#define CLASSNAME BLOCK_PREFIX(Z,Arnoldi_25_12,3)
#include "phist_gen_z.h"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME


#undef BLOCK_SIZE
#define BLOCK_SIZE 4

#ifdef PHIST_HAVE_SP

#define CLASSNAME BLOCK_PREFIX(S,Arnoldi_25_12,4)
#include "phist_gen_s.h"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME

#define CLASSNAME BLOCK_PREFIX(C,Arnoldi_25_12,4)
#include "phist_gen_c.h"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME

#endif

#define CLASSNAME BLOCK_PREFIX(D,Arnoldi_25_12,4)
#include "phist_gen_d.h"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME

#define CLASSNAME BLOCK_PREFIX(Z,Arnoldi_25_12,4)
#include "phist_gen_z.h"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME


