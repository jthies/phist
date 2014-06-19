#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#include "phist_macros.h"
#include "phist_kernels.h"
#include "phist_enums.h"
#include "phist_simple_arnoldi.h"
#include "../kernels/KernelTest.h"
#include "../kernels/KernelTestWithMap.h"
#include "../kernels/KernelTestWithType.h"
#include "../kernels/KernelTestWithSdMats.h"
#include "../kernels/KernelTestWithVectors.h"

using namespace testing;

// we can't easily store an all-zero matrix in ghost binCRS, so skip these tests
#ifdef PHIST_KERNEL_LIB_GHOST
#define SKIP_ZERO_MAT
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

#ifdef PHIST_HAVE_SP

#define CLASSNAME STestSimpleBlock2Arnoldi_25_12
#include "phist_gen_s.h"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME

#define CLASSNAME CTestSimpleBlock2Arnoldi_25_12
#include "phist_gen_c.h"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME

#endif

#define CLASSNAME DTestSimpleBlock2Arnoldi_25_12
#include "phist_gen_d.h"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME

#define CLASSNAME ZTestSimpleBlock2Arnoldi_25_12
#include "phist_gen_z.h"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME


#undef BLOCK_SIZE
#define BLOCK_SIZE 3

#ifdef PHIST_HAVE_SP

#define CLASSNAME STestSimpleBlock3Arnoldi_25_12
#include "phist_gen_s.h"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME

#define CLASSNAME CTestSimpleBlock3Arnoldi_25_12
#include "phist_gen_c.h"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME

#endif

#define CLASSNAME DTestSimpleBlock3Arnoldi_25_12
#include "phist_gen_d.h"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME

#define CLASSNAME ZTestSimpleBlock3Arnoldi_25_12
#include "phist_gen_z.h"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME


#undef BLOCK_SIZE
#define BLOCK_SIZE 4

#ifdef PHIST_HAVE_SP

#define CLASSNAME STestSimpleBlock4Arnoldi_25_12
#include "phist_gen_s.h"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME

#define CLASSNAME CTestSimpleBlock4Arnoldi_25_12
#include "phist_gen_c.h"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME

#endif

#define CLASSNAME DTestSimpleBlock4Arnoldi_25_12
#include "phist_gen_d.h"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME

#define CLASSNAME ZTestSimpleBlock4Arnoldi_25_12
#include "phist_gen_z.h"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME


