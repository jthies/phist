
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

#ifdef CLASSNAME
#undef CLASSNAME
#endif

#define _N_ 10
#define _M_ 1
#define _K_ 1

#define CLASSNAME SMvecSdMatTest_10_1
#include "phist_gen_s.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME DMvecSdMatTest_10_1

#include "phist_gen_d.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatTest_10_1

#include "phist_gen_c.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatTest_10_1

#include "phist_gen_z.h"
#include "MvecSdMatTest_def.hpp"

#undef _N_
#define _N_ 64
#undef _M_
#define _M_ 12

#undef CLASSNAME
#define CLASSNAME SMvecSdMatTest_64_12

#include "phist_gen_s.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME DMvecSdMatTest_64_12

#include "phist_gen_d.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecSdMatTest_64_12

#include "phist_gen_c.h"
#include "MvecSdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecSdMatTest_64_12

#include "phist_gen_z.h"
#include "MvecSdMatTest_def.hpp"

