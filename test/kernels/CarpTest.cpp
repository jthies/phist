#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/phist_gtest.h"


#include "kernels/phist_kernels.h"
#include "kernels/phist_kernel_flags.h"
#include "KernelTestWithSparseMat.h"
#include "KernelTestWithVectors.h"

#include <complex>

using namespace ::testing;

#define CLASSFILE_DEF "CarpTest_def.hpp"

#define MATNAME MATNAME_BENCH3D_8_A1
#define _BASENAME_ CarpTest_A
/* the tests are only executed if the carp kernel is implemented,
if the setup routine returns -99, further tests will not be run
*/

#define _N_ 512
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 512
#define _M_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 512
#define _M_ 7
#include "../phist_typed_test_gen.h"

#undef MATNAME
#define MATNAME MATNAME_IDFUNC
#undef _BASENAME_
#define _BASENAME_ CarpTest_I

#define _N_ 99
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 99
#define _M_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 99
#define _M_ 7
#include "../phist_typed_test_gen.h"
