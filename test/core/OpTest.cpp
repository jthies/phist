#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"

#include "phist_operator.h"
#include "phist_jadaOp.hpp"
#include "phist_orthog.h"
#include "../tools/MatrixIO.h"

#include "../kernels/KernelTest.h"
#include "../kernels/KernelTestWithMap.h"
#include "../kernels/KernelTestWithType.h"
#include "../kernels/KernelTestWithVectors.h"

#ifdef PHIST_KERNEL_LIB_TPETRA
#include "phist_tpetra_typedefs.hpp"
#include "BelosTpetraAdapter.hpp"
#elif defined(PHIST_KERNEL_LIB_EPETRA)
#include "Epetra_MultiVector.h"
#include "BelosEpetraAdapter.hpp"
#elif defined(PHIST_KERNEL_LIB_GHOST)
#ifdef PHIST_HAVE_BELOS
#include "Belos_GhostAdapter.hpp"
#endif
#endif

#ifdef PHIST_HAVE_BELOS
#include "phist_rcp_helpers.hpp"
#include "phist_BelosOperatorTraits.hpp"
#include "BelosMVOPTester.hpp"
#endif

using namespace testing;

#define _BASENAME_ OpTest

#define _N_ 25
#define _M_ 8
#include "../phist_typed_test_gen.h"

