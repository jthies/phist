#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"


#include "phist_macros.h"
#include "phist_kernels.h"
#include "phist_enums.h"
#include "phist_blockedgmres.h"
#include "phist_blockedminres.h"
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

#define CLASSNAME STestBlockedGMRES
#include "phist_gen_s.h"
#include "TestBlockedGMRES_def.hpp"
#undef CLASSNAME

#endif

#define CLASSNAME DTestBlockedGMRES
#include "phist_gen_d.h"
#include "TestBlockedGMRES_def.hpp"
#undef CLASSNAME

#undef TOLA
#undef TOLB
//TODO - tune these settings
#define TOLA 2.5e-2
#define TOLB 1.0e-4

#ifdef PHIST_HAVE_SP

#define CLASSNAME CTestBlockedGMRES
#include "phist_gen_c.h"
#include "TestBlockedGMRES_def.hpp"
#undef CLASSNAME

#endif

#define CLASSNAME ZTestBlockedGMRES
#include "phist_gen_z.h"
#include "TestBlockedGMRES_def.hpp"


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

#define MAXBAS 12
//TODO - tune these settings
#define TOLA 2.5e-2
#define TOLB 1.0e-3

#ifdef PHIST_HAVE_SP

#define CLASSNAME STestBlockedGMRES_symmMat
#include "phist_gen_s.h"
#include "TestBlockedGMRES_def.hpp"
#undef CLASSNAME

#endif

#define CLASSNAME DTestBlockedGMRES_symmMat
#include "phist_gen_d.h"
#include "TestBlockedGMRES_def.hpp"
#undef CLASSNAME

#undef TOLA
#undef TOLB
//TODO - tune these settings
#define TOLA 2.5e-2
#define TOLB 1.0e-3

#ifdef PHIST_HAVE_SP

#define CLASSNAME CTestBlockedGMRES_symmMat
#include "phist_gen_c.h"
#include "TestBlockedGMRES_def.hpp"
#undef CLASSNAME

#endif

#define CLASSNAME ZTestBlockedGMRES_symmMat
#include "phist_gen_z.h"
#include "TestBlockedGMRES_def.hpp"



