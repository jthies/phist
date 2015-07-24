#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#include "phist_kernels.h"
#include "../kernels/KernelTest.h"
#include "../kernels/KernelTestWithMap.h"
#include "../kernels/KernelTestWithType.h"
#include "../kernels/KernelTestWithVectors.h"
#include "../kernels/KernelTestWithSdMats.h"


#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

using namespace testing;

// Tests with V NxM, W NxK, C MxK, D KxK

///////////////////////////////////////////////
// case k=1                                  //
///////////////////////////////////////////////
#ifdef CLASSNAME
#undef CLASSNAME
#endif

#define _N_ 10
#define _M_ 1
#define _K_ 1

#ifdef PHIST_HAVE_SP
#define CLASSNAME SMvecSdMatFusedTest_10_1_1
#include "phist_gen_s.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatFusedTest_10_1_1

#include "phist_gen_c.h"
#include "MvecSdMatFusedTest_def.hpp"
#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatFusedTest_10_1_1

#include "phist_gen_d.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatFusedTest_10_1_1

#include "phist_gen_z.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef _N_
#define _N_ 64
#undef _M_
#define _M_ 12

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatFusedTest_64_12_1

#include "phist_gen_s.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatFusedTest_64_12_1

#include "phist_gen_c.h"
#include "MvecSdMatFusedTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatFusedTest_64_12_1

#include "phist_gen_d.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatFusedTest_64_12_1

#include "phist_gen_z.h"
#include "MvecSdMatFusedTest_def.hpp"

///////////////////////////////////////////////
// case k=4                                  //
///////////////////////////////////////////////
#ifdef CLASSNAME
#undef CLASSNAME
#endif

#undef _N_
#undef _M_
#undef _K_

#define _N_ 10
#define _M_ 1
#define _K_ 4

#ifdef PHIST_HAVE_SP
#define CLASSNAME SMvecSdMatFusedTest_10_1_4
#include "phist_gen_s.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatFusedTest_10_1_4

#include "phist_gen_c.h"
#include "MvecSdMatFusedTest_def.hpp"
#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatFusedTest_10_1_4

#include "phist_gen_d.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatFusedTest_10_1_4

#include "phist_gen_z.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef _N_
#define _N_ 64
#undef _M_
#define _M_ 12

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatFusedTest_64_12_4

#include "phist_gen_s.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatFusedTest_64_12_4

#include "phist_gen_c.h"
#include "MvecSdMatFusedTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatFusedTest_64_12_4

#include "phist_gen_d.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatFusedTest_64_12_4

#include "phist_gen_z.h"
#include "MvecSdMatFusedTest_def.hpp"

// case with square sdMats (m==k)

#undef _M_
#define _M_ 4

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatFusedTest_64_4_4

#include "phist_gen_s.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatFusedTest_64_4_4

#include "phist_gen_c.h"
#include "MvecSdMatFusedTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatFusedTest_64_4_4

#include "phist_gen_d.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatFusedTest_64_4_4

#include "phist_gen_z.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef _M_
#define _M_ 1
#undef _K_
#define _K_ 2

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatFusedTest_64_1_2

#include "phist_gen_s.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatFusedTest_64_1_2

#include "phist_gen_c.h"
#include "MvecSdMatFusedTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatFusedTest_64_1_2

#include "phist_gen_d.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatFusedTest_64_1_2

#include "phist_gen_z.h"
#include "MvecSdMatFusedTest_def.hpp"


///////////////////////////////////////////////
// case k=5                                  //
///////////////////////////////////////////////
#ifdef CLASSNAME
#undef CLASSNAME
#endif

#undef _N_
#undef _M_
#undef _K_

#define _N_ 10
#define _M_ 1
#define _K_ 5

#ifdef PHIST_HAVE_SP
#define CLASSNAME SMvecSdMatFusedTest_10_1_5
#include "phist_gen_s.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatFusedTest_10_1_5

#include "phist_gen_c.h"
#include "MvecSdMatFusedTest_def.hpp"
#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatFusedTest_10_1_5

#include "phist_gen_d.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatFusedTest_10_1_5

#include "phist_gen_z.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef _N_
#define _N_ 64
#undef _M_
#define _M_ 12

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatFusedTest_64_12_5

#include "phist_gen_s.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatFusedTest_64_12_5

#include "phist_gen_c.h"
#include "MvecSdMatFusedTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatFusedTest_64_12_5

#include "phist_gen_d.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatFusedTest_64_12_5

#include "phist_gen_z.h"
#include "MvecSdMatFusedTest_def.hpp"

// case with square sdMats (m=k=5)

#undef _M_
#define _M_ 5

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatFusedTest_64_5_5

#include "phist_gen_s.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatFusedTest_64_5_5

#include "phist_gen_c.h"
#include "MvecSdMatFusedTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatFusedTest_64_5_5

#include "phist_gen_d.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatFusedTest_64_5_5

#include "phist_gen_z.h"
#include "MvecSdMatFusedTest_def.hpp"

// case with 7, 11

#undef _N_
#define _N_ 512
#undef _M_
#define _M_ 7
#undef _K_
#define _K_ 11

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatFusedTest_512_7_11

#include "phist_gen_s.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatFusedTest_512_7_11

#include "phist_gen_c.h"
#include "MvecSdMatFusedTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatFusedTest_512_7_11

#include "phist_gen_d.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatFusedTest_512_7_11

#include "phist_gen_z.h"
#include "MvecSdMatFusedTest_def.hpp"


// large tests for verifying high precision calculations
#undef _N_
#ifdef PHIST_HIGH_PRECISION_KERNELS
#define _N_ 4000000
#undef _M_
#define _M_ 1
#undef _K_
#define _K_ 1

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatFusedTest_4M_1_1

#include "phist_gen_s.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatFusedTest_4M_1_1

#include "phist_gen_c.h"
#include "MvecSdMatFusedTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatFusedTest_4M_1_1

#include "phist_gen_d.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatFusedTest_4M_1_1

#include "phist_gen_z.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef _K_
#define _K_ 2

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatFusedTest_4M_1_2

#include "phist_gen_s.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatFusedTest_4M_1_2

#include "phist_gen_c.h"
#include "MvecSdMatFusedTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatFusedTest_4M_1_2

#include "phist_gen_d.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatFusedTest_4M_1_2

#include "phist_gen_z.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef _K_
#define _K_ 4

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatFusedTest_4M_1_4

#include "phist_gen_s.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatFusedTest_4M_1_4

#include "phist_gen_c.h"
#include "MvecSdMatFusedTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatFusedTest_4M_1_4

#include "phist_gen_d.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatFusedTest_4M_1_4

#include "phist_gen_z.h"
#include "MvecSdMatFusedTest_def.hpp"


#undef _M_
#define _M_ 2
#undef _K_
#define _K_ 1

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatFusedTest_4M_2_1

#include "phist_gen_s.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatFusedTest_4M_2_1

#include "phist_gen_c.h"
#include "MvecSdMatFusedTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatFusedTest_4M_2_1

#include "phist_gen_d.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatFusedTest_4M_2_1

#include "phist_gen_z.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef _K_
#define _K_ 2

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatFusedTest_4M_2_2

#include "phist_gen_s.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatFusedTest_4M_2_2

#include "phist_gen_c.h"
#include "MvecSdMatFusedTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatFusedTest_4M_2_2

#include "phist_gen_d.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatFusedTest_4M_2_2

#include "phist_gen_z.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef _K_
#define _K_ 4

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatFusedTest_4M_2_4

#include "phist_gen_s.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatFusedTest_4M_2_4

#include "phist_gen_c.h"
#include "MvecSdMatFusedTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatFusedTest_4M_2_4

#include "phist_gen_d.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatFusedTest_4M_2_4

#include "phist_gen_z.h"
#include "MvecSdMatFusedTest_def.hpp"


#undef _M_
#define _M_ 4
#undef _K_
#define _K_ 1

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatFusedTest_4M_4_1

#include "phist_gen_s.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatFusedTest_4M_4_1

#include "phist_gen_c.h"
#include "MvecSdMatFusedTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatFusedTest_4M_4_1

#include "phist_gen_d.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatFusedTest_4M_4_1

#include "phist_gen_z.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef _K_
#define _K_ 2

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatFusedTest_4M_4_2

#include "phist_gen_s.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatFusedTest_4M_4_2

#include "phist_gen_c.h"
#include "MvecSdMatFusedTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatFusedTest_4M_4_2

#include "phist_gen_d.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatFusedTest_4M_4_2

#include "phist_gen_z.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef _K_
#define _K_ 4

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatFusedTest_4M_4_4

#include "phist_gen_s.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatFusedTest_4M_4_4

#include "phist_gen_c.h"
#include "MvecSdMatFusedTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatFusedTest_4M_4_4

#include "phist_gen_d.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatFusedTest_4M_4_4

#include "phist_gen_z.h"
#include "MvecSdMatFusedTest_def.hpp"



#undef _M_
#define _M_ 10
#undef _K_
#define _K_ 1

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatFusedTest_4M_10_1

#include "phist_gen_s.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatFusedTest_4M_10_1

#include "phist_gen_c.h"
#include "MvecSdMatFusedTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatFusedTest_4M_10_1

#include "phist_gen_d.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatFusedTest_4M_10_1

#include "phist_gen_z.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef _K_
#define _K_ 2

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatFusedTest_4M_10_2

#include "phist_gen_s.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatFusedTest_4M_10_2

#include "phist_gen_c.h"
#include "MvecSdMatFusedTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatFusedTest_4M_10_2

#include "phist_gen_d.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatFusedTest_4M_10_2

#include "phist_gen_z.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef _K_
#define _K_ 4

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatFusedTest_4M_10_4

#include "phist_gen_s.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatFusedTest_4M_10_4

#include "phist_gen_c.h"
#include "MvecSdMatFusedTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatFusedTest_4M_10_4

#include "phist_gen_d.h"
#include "MvecSdMatFusedTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatFusedTest_4M_10_4

#include "phist_gen_z.h"
#include "MvecSdMatFusedTest_def.hpp"

#endif /* PHIST_HIGH_PRECISION_KERNELS */

