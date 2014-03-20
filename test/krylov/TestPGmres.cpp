
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#include "phist_macros.h"
#include "phist_kernels.h"
#include "phist_enums.h"
#include "phist_pgmres.h"
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
#define _M_ 6

#define MAXBAS 25
//TODO - tune these settings
#define TOLA 2.5e-2
#define TOLB 1.0e-4

#ifdef PHIST_HAVE_SP

#define CLASSNAME STestPGmres
#include "phist_gen_s.h"
#include "TestPGmres_def.hpp"
#undef CLASSNAME

#endif

#define CLASSNAME DTestPGmres
#include "phist_gen_d.h"
#include "TestPGmres_def.hpp"
#undef CLASSNAME

#undef TOLA
#undef TOLB
//TODO - tune these settings
#define TOLA 2.5e-2
#define TOLB 1.0e-4

#ifdef PHIST_HAVE_SP

#define CLASSNAME CTestPGmres
#include "phist_gen_c.h"
#include "TestPGmres_def.hpp"
#undef CLASSNAME

#endif

#define CLASSNAME ZTestPGmres
#include "phist_gen_z.h"
#include "TestPGmres_def.hpp"

