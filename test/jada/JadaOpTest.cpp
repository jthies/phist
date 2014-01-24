#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#include "phist_jadaOp.hpp"
#include "../tools/MatrixIO.h"

#include "../kernels/KernelTest.h"
#include "../kernels/KernelTestWithMap.h"
#include "../kernels/KernelTestWithType.h"
#include "../kernels/KernelTestWithVectors.h"
#include "../kernels/KernelTestWithSdMats.h"

using namespace testing;

#ifdef CLASSNAME
#undef CLASSNAME
#endif

#define _N_ 25
#define _NV_ 4
#define _NVP_ 10

#define CLASSNAME SJadaOpTest_25_4_10
#include "phist_gen_s.h"
#include "JadaOpTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME DJadaOpTest_25_4_10

#include "phist_gen_d.h"
#include "JadaOpTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CJadaOpTest_25_4_10

#include "phist_gen_c.h"
#include "JadaOpTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZJadaOpTest_25_4_10

#include "phist_gen_d.h"
#include "JadaOpTest_def.hpp"


