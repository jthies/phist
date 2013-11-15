
#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#include "phist_macros.h"
#include "phist_kernels.h"
#include "phist_enums.h"
#include "phist_simple_arnoldi.h"
#include "../kernels/KernelTest.h"
#include "../kernels/KernelTestWithMap.h"
#include "../kernels/KernelTestWithType.h"
#include "../kernels/KernelTestWithSdMats.h"
#include "../kernels/KernelTestWithVectors.h"

using namespace testing;

#ifdef CLASSNAME
#undef CLASSNAME
#endif

#define _N_ 25
#define _M_ 10

#define CLASSNAME STestSimpleArnoldi_25_10
#include "phist_gen_s.h"
#include "TestSimpleArnoldi_def.hpp"

#undef CLASSNAME
#define CLASSNAME DTestSimpleArnoldi_25_10

#include "phist_gen_d.h"
#include "TestSimpleArnoldi_def.hpp"

#undef CLASSNAME
#define CLASSNAME CTestSimpleArnoldi_25_10

#include "phist_gen_c.h"
#include "TestSimpleArnoldi_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZTestSimpleArnoldi_25_10

#include "phist_gen_z.h"
#include "TestSimpleArnoldi_def.hpp"

