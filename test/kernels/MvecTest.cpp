#include "phist_config.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#include "phist_kernels.h"
#include "KernelTestWithMap.h"
#include "KernelTestWithVectors.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

using namespace testing;

#define _N_ 9
#define _NV_ 1
#ifdef CLASSNAME
#undef CLASSNAME
#endif
#ifdef PHIST_HAVE_SP
#define CLASSNAME SMvecTest_9_1
#include "phist_gen_s.h"
#include "MvecTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecTest_9_1

#include "phist_gen_c.h"
#include "MvecTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecTest_9_1

#include "phist_gen_d.h"
#include "MvecTest_def.hpp"


#undef CLASSNAME
#define CLASSNAME ZMvecTest_9_1

#include "phist_gen_z.h"
#include "MvecTest_def.hpp"

#undef _N_
#define _N_ 16
#undef _NV_
#define _NV_ 9

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecTest_16_9

#include "phist_gen_s.h"
#include "MvecTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecTest_16_9

#include "phist_gen_c.h"
#include "MvecTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecTest_16_9

#include "phist_gen_d.h"
#include "MvecTest_def.hpp"


#undef CLASSNAME
#define CLASSNAME ZMvecTest_16_9

#include "phist_gen_z.h"
#include "MvecTest_def.hpp"


#undef _N_
/* 4000000 are required for my accuracy test ! */
#define _N_ 4000000
#undef _NV_
#define _NV_ 1

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecTest_4M_1

#include "phist_gen_s.h"
#include "MvecTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecTest_4M_1

#include "phist_gen_c.h"
#include "MvecTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecTest_4M_1

#include "phist_gen_d.h"
#include "MvecTest_def.hpp"


#undef CLASSNAME
#define CLASSNAME ZMvecTest_4M_1

#include "phist_gen_z.h"
#include "MvecTest_def.hpp"


#undef _N_
#define _N_ 4000000
#undef _NV_
#define _NV_ 2

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecTest_4M_2

#include "phist_gen_s.h"
#include "MvecTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecTest_4M_2

#include "phist_gen_c.h"
#include "MvecTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecTest_4M_2

#include "phist_gen_d.h"
#include "MvecTest_def.hpp"


#undef CLASSNAME
#define CLASSNAME ZMvecTest_4M_2

#include "phist_gen_z.h"
#include "MvecTest_def.hpp"


#undef _N_
#define _N_ 4000000
#undef _NV_
#define _NV_ 4

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecTest_4M_4

#include "phist_gen_s.h"
#include "MvecTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecTest_4M_4

#include "phist_gen_c.h"
#include "MvecTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecTest_4M_4

#include "phist_gen_d.h"
#include "MvecTest_def.hpp"


#undef CLASSNAME
#define CLASSNAME ZMvecTest_4M_4

#include "phist_gen_z.h"
#include "MvecTest_def.hpp"


// the Belos Tester selects the number of vectors itself,
// it is enough to test it for one vector length, I think
#ifdef PHIST_KERNEL_LIB_GHOST
#ifdef PHIST_HAVE_BELOS
// ghost views of mvecs are not general enough to pass these
// tests for row-major ordering
#ifndef PHIST_MVECS_ROW_MAJOR
#define DO_BELOS_TESTS
#include "phist_GhostMV.hpp"
#include "phist_rcp_helpers.hpp"
#include "Belos_GhostAdapter.hpp"
#include "BelosMVOPTester.hpp"
#include "BelosOutputManager.hpp"
#endif
#endif
#endif


#undef _N_
#define _N_ 237
#undef _NV_
#define _NV_ 4

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecTest_237_4

#include "phist_gen_s.h"
#include "MvecTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecTest_237_4

#include "phist_gen_c.h"
#include "MvecTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecTest_237_4

#include "phist_gen_d.h"
#include "MvecTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecTest_237_4

#include "phist_gen_z.h"
#include "MvecTest_def.hpp"

#ifdef DO_BELOS_TESTS
#undef DO_BELOS_TESTS
#endif
