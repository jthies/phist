#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/phist_gtest.h"

#include "kernels/phist_kernels.h"
#include "KernelTestWithMap.h"
#include "KernelTestWithVectors.h"
#include "../tools/MatrixIO.h"


using namespace ::testing;

#define _BASENAME_ MultipleMatrixTest
#define _USE_VIEWS_ 0
#define _N_ 125
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 125
#define _M_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 125
#define _M_ 7
#include "../phist_typed_test_gen.h"

