#include "phist_config.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/phist_gtest.h"
//#include "gmock/gmock.h"

#include "phist_macros.h"
#include "phist_kernels.h"
#include "phist_enums.h"
#include "phist_jadaOp.hpp"
#include "phist_jadaCorrectionSolver.h"
#include "phist_ScalarTraits.hpp"
#include "phist_orthog.h"
#include "../kernels/KernelTest.h"
#include "../kernels/KernelTestWithMap.h"
#include "../kernels/KernelTestWithType.h"
#include "../kernels/KernelTestWithSdMats.h"
#include "../kernels/KernelTestWithVectors.h"
#include "JadaTestWithOpts.hpp"

#include "../tools/MatrixIO.h"

using namespace testing;

#ifdef CLASSNAME
#undef CLASSNAME
#endif

#ifdef MATNAME
#undef MATNAME
#endif

#define MATNAME "sprandn"

#define _N_ 25
#define _NV_ 10
#define _NVP_ 2
#define _MAXBAS_ 5

#ifdef PHIST_HAVE_SP

#define CLASSNAME STestJadaCorrectionSolver_25_10_2
#include "phist_gen_s.h"
#include "TestJadaCorrectionSolver_def.hpp"

#undef CLASSNAME
#define CLASSNAME CTestJadaCorrectionSolver_25_10_2

#include "phist_gen_c.h"
#include "TestJadaCorrectionSolver_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DTestJadaCorrectionSolver_25_10_2

#include "phist_gen_d.h"
#include "TestJadaCorrectionSolver_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZTestJadaCorrectionSolver_25_10_2

#include "phist_gen_z.h"
#include "TestJadaCorrectionSolver_def.hpp"

#undef MATNAME
#define MATNAME "jadaTestMat"

#undef _N_
#undef _NV_
#undef _NVP_
#undef _MAXBAS_

#define _N_ 512
#define _NV_ 11
#define _NVP_ 7
#define _MAXBAS_ 10

#ifdef PHIST_HAVE_SP

#undef CLASSNAME
#define CLASSNAME STestJadaCorrectionSolver_512_11_7
#include "phist_gen_s.h"
#include "TestJadaCorrectionSolver_def.hpp"

#undef CLASSNAME
#define CLASSNAME CTestJadaCorrectionSolver_512_11_7

#include "phist_gen_c.h"
#include "TestJadaCorrectionSolver_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DTestJadaCorrectionSolver_512_11_7

#include "phist_gen_d.h"
#include "TestJadaCorrectionSolver_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZTestJadaCorrectionSolver_512_11_7

#include "phist_gen_z.h"
#include "TestJadaCorrectionSolver_def.hpp"


