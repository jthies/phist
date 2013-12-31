#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#include "phist_macros.h"
#include "phist_kernels.h"
#include "phist_enums.h"
#include "phist_jadaInnerGmres.h"
#include "phist_jadaOp.h"
#include "phist_ScalarTraits.hpp"
#include "../kernels/KernelTest.h"
#include "../kernels/KernelTestWithMap.h"
#include "../kernels/KernelTestWithType.h"
#include "../kernels/KernelTestWithSdMats.h"
#include "../kernels/KernelTestWithVectors.h"

#include "../tools/MatrixIO.h"

using namespace testing;

#ifdef CLASSNAME
#undef CLASSNAME
#endif

#ifdef MATNAME
#undef MATNAME
#endif

#define MATNAME "jadaTestMat"

#define _N_ 512
#define _NV_ 3
#define _NVP_ 7
#define _MAXBAS_ 40

#define CLASSNAME SJadaInnerGmresTest_512_3_7
#include "phist_gen_s.h"
#include "JadaInnerGmresTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME DJadaInnerGmresTest_512_3_7

#include "phist_gen_d.h"
#include "JadaInnerGmresTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CJadaInnerGmresTest_512_3_7

#include "phist_gen_c.h"
#include "JadaInnerGmresTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZJadaInnerGmresTest_512_3_7

#include "phist_gen_z.h"
#include "JadaInnerGmresTest_def.hpp"


