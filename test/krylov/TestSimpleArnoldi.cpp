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
#define _M_ 10

#ifdef PHIST_HAVE_SP

#define CLASSNAME STestSimpleArnoldi_25_10
#include "phist_gen_s.h"
#include "../tools/MatrixIO_def.hpp"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME

#define CLASSNAME CTestSimpleArnoldi_25_10
#include "phist_gen_c.h"
#include "../tools/MatrixIO_def.hpp"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME

#endif

#define CLASSNAME DTestSimpleArnoldi_25_10
#include "phist_gen_d.h"
#include "../tools/MatrixIO_def.hpp"
#include "TestSimpleArnoldi_def.hpp"
#undef CLASSNAME

#define CLASSNAME ZTestSimpleArnoldi_25_10
#include "phist_gen_z.h"
#include "../tools/MatrixIO_def.hpp"
#include "TestSimpleArnoldi_def.hpp"

