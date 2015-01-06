#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#include "phist_macros.h"
#include "phist_kernels.h"
#include "phist_enums.h"
#include "phist_carp_cg.h"
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

//TODO - tune these settings
#define TOLA 1.0e-3
#define TOLB 1.0e-4

#ifdef PHIST_HAVE_SP

#define CLASSNAME STestCarpCG
#include "phist_gen_s.h"
#include "TestCarpCG_def.hpp"
#undef CLASSNAME

#endif

#define CLASSNAME DTestCarpCG
#include "phist_gen_d.h"
#include "TestCarpCG_def.hpp"
#undef CLASSNAME

#undef TOLA
#undef TOLB

#define TOLA 1.0e-3
#define TOLB 1.0e-4

#ifdef PHIST_HAVE_SP

#define CLASSNAME CTestCarpCG
#include "phist_gen_c.h"
#include "TestCarpCG_def.hpp"
#undef CLASSNAME

#endif

#define CLASSNAME ZTestCarpCG
#include "phist_gen_z.h"
#include "TestCarpCG_def.hpp"


#undef MATNAME
#undef MAXBAS
#undef _N_
#undef _M_
#undef TOLA
#undef TOLB



#undef CLASSNAME
#define MATNAME "symmMat"
#define MATSYMMETRIC
#define _N_ 20
#define _M_ 5

//TODO - tune these settings
#define TOLA 1.0e-3
#define TOLB 1.0e-4

#ifdef PHIST_HAVE_SP

#define CLASSNAME STestCarpCG_symmMat
#include "phist_gen_s.h"
#include "TestCarpCG_def.hpp"
#undef CLASSNAME

#endif

#define CLASSNAME DTestCarpCG_symmMat
#include "phist_gen_d.h"
#include "TestCarpCG_def.hpp"
#undef CLASSNAME

#undef TOLA
#undef TOLB
//TODO - tune these settings
#define TOLA 1.0e-3
#define TOLB 1.0e-4

#ifdef PHIST_HAVE_SP

#define CLASSNAME CTestCarpCG_symmMat
#include "phist_gen_c.h"
#include "TestCarpCG_def.hpp"
#undef CLASSNAME

#endif

#define CLASSNAME ZTestCarpCG_symmMat
#include "phist_gen_z.h"
#include "TestCarpCG_def.hpp"



