
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#include "phist_macros.h"
#include "phist_kernels.h"
#include "phist_enums.h"
#include "phist_gmres.h"
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

#define MAXBAS 20
//TODO - tune these settings
#define TOLA 5.0e-2
#define TOLB 1.0e-4

#define CLASSNAME STestGmres
#include "phist_gen_s.h"
#include "TestGmres_def.hpp"

#undef CLASSNAME
#define CLASSNAME DTestGmres

#include "phist_gen_d.h"
#include "TestGmres_def.hpp"

#undef CLASSNAME
#undef TOLA
#undef TOLB
//TODO - tune these settings
#define TOLA 1.0e-3
#define TOLB 1.0e-5
#define CLASSNAME CTestGmres

#include "phist_gen_c.h"
#include "TestGmres_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZTestGmres

#include "phist_gen_z.h"
#include "TestGmres_def.hpp"

