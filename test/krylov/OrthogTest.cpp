
#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#include "phist_kernels.h"
#include "../kernels/KernelTest.h"
#include "../kernels/KernelTest.h"
#include "../kernels/KernelTestWithMap.h"
#include "../kernels/KernelTestWithType.h"
#include "../kernels/KernelTestWithVectors.h"
#include "../kernels/KernelTestWithSdMats.h"


#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

using namespace testing;

#ifdef CLASSNAME
#undef CLASSNAME
#endif

#define _N_ 8
#define _M_ 1
#define _K_ 1

#define CLASSNAME SOrthogTest_8_1_1
#include "phist_gen_s.h"
#include "OrthogTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME DOrthogTest_8_1_1

#include "phist_gen_d.h"
#include "OrthogTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME COrthogTest_8_1_1

#include "phist_gen_c.h"
#include "OrthogTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZOrthogTest_8_1_1

#include "phist_gen_z.h"
#include "OrthogTest_def.hpp"

#undef _N_
#define _N_ 15
#undef _M_
#define _M_ 5
#undef _K_
#define _K_ 3

#undef CLASSNAME
#define CLASSNAME SOrthogTest_15_5_3

#include "phist_gen_s.h"
#include "OrthogTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME DOrthogTest_15_5_3

#include "phist_gen_d.h"
#include "OrthogTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME COrthogTest_15_5_3

#include "phist_gen_c.h"
#include "OrthogTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZOrthogTest_15_5_3

#include "phist_gen_z.h"
#include "OrthogTest_def.hpp"

