
#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#include "phist_kernels.h"
#include "KernelTestWithMap.h"
#include "KernelTestWithVectors.h"

#ifdef PHIST_KERNEL_LIB_GHOST
#define DO_BELOS_TESTS
#include "phist_GhostMV.hpp"
#include "phist_rcp_helpers.hpp"
#include "Belos_GhostAdapter.hpp"
#include "BelosMVOPTester.hpp"
#include "BelosOutputManager.hpp"
#endif

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

using namespace testing;

#define _N_ 9
#define _NV_ 1
#define DO_BELOS_TESTS
#ifdef CLASSNAME
#undef CLASSNAME
#endif
#define CLASSNAME SMvecTest_9_1
#include "phist_gen_s.h"
#include "MvecTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME DMvecTest_9_1

#include "phist_gen_d.h"
#include "MvecTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecTest_9_1

#include "phist_gen_c.h"
#include "MvecTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecTest_9_1

#include "phist_gen_z.h"
#include "MvecTest_def.hpp"

// the Belos Tester selects the number of vectors itself,
// it is enough to test it for one vector length, I think
#undef DO_BELOS_TESTS

#undef _N_
#define _N_ 6
#undef _NV_
#define _NV_ 3

#undef CLASSNAME
#define CLASSNAME SMvecTest_6_3

#include "phist_gen_s.h"
#include "MvecTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME DMvecTest_6_3

#include "phist_gen_d.h"
#include "MvecTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecTest_6_3

#include "phist_gen_c.h"
#include "MvecTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecTest_6_3

#include "phist_gen_z.h"
#include "MvecTest_def.hpp"
