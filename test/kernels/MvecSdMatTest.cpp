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

// Tests with V NxM, W NxK, C MxK

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
#define CLASSNAME SMvecSdMatTest_10_1_1
#include "phist_gen_s.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatTest_10_1_1

#include "phist_gen_c.h"
#include "MvecSdMatTest_def.hpp"
#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatTest_10_1_1

#include "phist_gen_d.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatTest_10_1_1

#include "phist_gen_z.h"
#include "MvecSdMatTest_def.hpp"

#undef _N_
#define _N_ 64
#undef _M_
#define _M_ 12

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatTest_64_12_1

#include "phist_gen_s.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatTest_64_12_1

#include "phist_gen_c.h"
#include "MvecSdMatTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatTest_64_12_1

#include "phist_gen_d.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatTest_64_12_1

#include "phist_gen_z.h"
#include "MvecSdMatTest_def.hpp"

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
#define CLASSNAME SMvecSdMatTest_10_1_4
#include "phist_gen_s.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatTest_10_1_4

#include "phist_gen_c.h"
#include "MvecSdMatTest_def.hpp"
#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatTest_10_1_4

#include "phist_gen_d.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatTest_10_1_4

#include "phist_gen_z.h"
#include "MvecSdMatTest_def.hpp"

#undef _N_
#define _N_ 64
#undef _M_
#define _M_ 12

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatTest_64_12_4

#include "phist_gen_s.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatTest_64_12_4

#include "phist_gen_c.h"
#include "MvecSdMatTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatTest_64_12_4

#include "phist_gen_d.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatTest_64_12_4

#include "phist_gen_z.h"
#include "MvecSdMatTest_def.hpp"

// case with square sdMats (m==k)

#undef _M_
#define _M_ 4

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatTest_64_4_4

#include "phist_gen_s.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatTest_64_4_4

#include "phist_gen_c.h"
#include "MvecSdMatTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatTest_64_4_4

#include "phist_gen_d.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatTest_64_4_4

#include "phist_gen_z.h"
#include "MvecSdMatTest_def.hpp"

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
#define CLASSNAME SMvecSdMatTest_10_1_5
#include "phist_gen_s.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatTest_10_1_5

#include "phist_gen_c.h"
#include "MvecSdMatTest_def.hpp"
#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatTest_10_1_5

#include "phist_gen_d.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatTest_10_1_5

#include "phist_gen_z.h"
#include "MvecSdMatTest_def.hpp"

#undef _N_
#define _N_ 64
#undef _M_
#define _M_ 12

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatTest_64_12_5

#include "phist_gen_s.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatTest_64_12_5

#include "phist_gen_c.h"
#include "MvecSdMatTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatTest_64_12_5

#include "phist_gen_d.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatTest_64_12_5

#include "phist_gen_z.h"
#include "MvecSdMatTest_def.hpp"

// case with square sdMats (m=k=5)

#undef _M_
#define _M_ 5

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatTest_64_5_5

#include "phist_gen_s.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatTest_64_5_5

#include "phist_gen_c.h"
#include "MvecSdMatTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatTest_64_5_5

#include "phist_gen_d.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatTest_64_5_5

#include "phist_gen_z.h"
#include "MvecSdMatTest_def.hpp"


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
#define CLASSNAME SMvecSdMatTest_4M_1_1

#include "phist_gen_s.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatTest_4M_1_1

#include "phist_gen_c.h"
#include "MvecSdMatTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatTest_4M_1_1

#include "phist_gen_d.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatTest_4M_1_1

#include "phist_gen_z.h"
#include "MvecSdMatTest_def.hpp"

#undef _K_
#define _K_ 2

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatTest_4M_1_2

#include "phist_gen_s.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatTest_4M_1_2

#include "phist_gen_c.h"
#include "MvecSdMatTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatTest_4M_1_2

#include "phist_gen_d.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatTest_4M_1_2

#include "phist_gen_z.h"
#include "MvecSdMatTest_def.hpp"

#undef _K_
#define _K_ 4

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatTest_4M_1_4

#include "phist_gen_s.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatTest_4M_1_4

#include "phist_gen_c.h"
#include "MvecSdMatTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatTest_4M_1_4

#include "phist_gen_d.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatTest_4M_1_4

#include "phist_gen_z.h"
#include "MvecSdMatTest_def.hpp"


#undef _M_
#define _M_ 2
#undef _K_
#define _K_ 1

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatTest_4M_2_1

#include "phist_gen_s.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatTest_4M_2_1

#include "phist_gen_c.h"
#include "MvecSdMatTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatTest_4M_2_1

#include "phist_gen_d.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatTest_4M_2_1

#include "phist_gen_z.h"
#include "MvecSdMatTest_def.hpp"

#undef _K_
#define _K_ 2

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatTest_4M_2_2

#include "phist_gen_s.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatTest_4M_2_2

#include "phist_gen_c.h"
#include "MvecSdMatTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatTest_4M_2_2

#include "phist_gen_d.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatTest_4M_2_2

#include "phist_gen_z.h"
#include "MvecSdMatTest_def.hpp"

#undef _K_
#define _K_ 4

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatTest_4M_2_4

#include "phist_gen_s.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatTest_4M_2_4

#include "phist_gen_c.h"
#include "MvecSdMatTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatTest_4M_2_4

#include "phist_gen_d.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatTest_4M_2_4

#include "phist_gen_z.h"
#include "MvecSdMatTest_def.hpp"


#undef _M_
#define _M_ 4
#undef _K_
#define _K_ 1

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatTest_4M_4_1

#include "phist_gen_s.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatTest_4M_4_1

#include "phist_gen_c.h"
#include "MvecSdMatTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatTest_4M_4_1

#include "phist_gen_d.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatTest_4M_4_1

#include "phist_gen_z.h"
#include "MvecSdMatTest_def.hpp"

#undef _K_
#define _K_ 2

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatTest_4M_4_2

#include "phist_gen_s.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatTest_4M_4_2

#include "phist_gen_c.h"
#include "MvecSdMatTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatTest_4M_4_2

#include "phist_gen_d.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatTest_4M_4_2

#include "phist_gen_z.h"
#include "MvecSdMatTest_def.hpp"

#undef _K_
#define _K_ 4

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatTest_4M_4_4

#include "phist_gen_s.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatTest_4M_4_4

#include "phist_gen_c.h"
#include "MvecSdMatTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatTest_4M_4_4

#include "phist_gen_d.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatTest_4M_4_4

#include "phist_gen_z.h"
#include "MvecSdMatTest_def.hpp"



#undef _M_
#define _M_ 10
#undef _K_
#define _K_ 1

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatTest_4M_10_1

#include "phist_gen_s.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatTest_4M_10_1

#include "phist_gen_c.h"
#include "MvecSdMatTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatTest_4M_10_1

#include "phist_gen_d.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatTest_4M_10_1

#include "phist_gen_z.h"
#include "MvecSdMatTest_def.hpp"

#undef _K_
#define _K_ 2

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatTest_4M_10_2

#include "phist_gen_s.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatTest_4M_10_2

#include "phist_gen_c.h"
#include "MvecSdMatTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatTest_4M_10_2

#include "phist_gen_d.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatTest_4M_10_2

#include "phist_gen_z.h"
#include "MvecSdMatTest_def.hpp"

#undef _K_
#define _K_ 4

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SMvecSdMatTest_4M_10_4

#include "phist_gen_s.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatTest_4M_10_4

#include "phist_gen_c.h"
#include "MvecSdMatTest_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DMvecSdMatTest_4M_10_4

#include "phist_gen_d.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatTest_4M_10_4

#include "phist_gen_z.h"
#include "MvecSdMatTest_def.hpp"

#endif /* PHIST_HIGH_PRECISION_KERNELS */

