
#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#include "phist_macros.h"
#include "phist_kernels.h"
#include "phist_enums.h"
#include "phist_schur_decomp.h"
#include "../kernels/KernelTest.h"
#include "../kernels/KernelTestWithType.h"
#include "../kernels/KernelTestWithSdMats.h"

using namespace testing;

#ifdef CLASSNAME
#undef CLASSNAME
#endif

#define _N_ 10

#define CLASSNAME STestSchurDecomp_10
#include "phist_gen_s.h"
#include "TestSchurDecomp_def.hpp"

#undef CLASSNAME
#define CLASSNAME DTestSchurDecomp_10

#include "phist_gen_d.h"
#include "TestSchurDecomp_def.hpp"

#undef CLASSNAME
#define CLASSNAME CTestSchurDecomp_10

#include "phist_gen_c.h"
#include "TestSchurDecomp_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZTestSchurDecomp_10

#include "phist_gen_z.h"
#include "TestSchurDecomp_def.hpp"
