
#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#include "phist_operator.h"
#include "phist_jadaOp.h"

#include "../kernels/KernelTest.h"
#include "../kernels/KernelTestWithMap.h"
#include "../kernels/KernelTestWithType.h"
#include "../kernels/KernelTestWithVectors.h"

#ifdef PHIST_KERNEL_LIB_TPETRA
#include "phist_tpetra_typedefs.hpp"
#include "BelosTpetraAdapter.hpp"
#elif defined(PHIST_KERNEL_LIB_EPETRA)
#include "BelosEpetraAdapter.hpp"
#elif defined(PHIST_KERNEL_LIB_GHOST)
#include "Belos_GhostAdapter.hpp"
#endif
#include "phist_rcp_helpers.hpp"
#include "phist_BelosOperatorTraits.hpp"
#include "BelosMVOPTester.hpp"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

using namespace testing;

#ifdef CLASSNAME
#undef CLASSNAME
#endif

#define _N_ 25
#define _NV_ 8

#define CLASSNAME SOpTest_25_8
#include "phist_gen_s.h"
#include "OpTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME DOpTest_25_8

#include "phist_gen_d.h"
#include "OpTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME COpTest_25_8

#include "phist_gen_c.h"
#include "OpTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZOpTest_25_8

#include "phist_gen_d.h"
#include "OpTest_def.hpp"

