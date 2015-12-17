#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"


#include "phist_macros.h"
#include "phist_kernels.h"
#include "phist_enums.h"
#include "phist_blockedgmres.h"
#include "phist_jadaOp.hpp"
#include "phist_ScalarTraits.hpp"
#include "phist_orthog.h"
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

#define MATNAME "sprandn"

#define _N_ 25
#define _NV_ 1
#define _NVP_ 2
#define _MAXBAS_ 5

#ifdef PHIST_HAVE_SP

# define CLASSNAME SJadaInnerGmresTest_25_1_2
# include "phist_gen_s.h"
# include "JadaInnerGmresTest_def.hpp"
# undef CLASSNAME

# define CLASSNAME CJadaInnerGmresTest_25_1_2
# include "phist_gen_c.h"
# include "JadaInnerGmresTest_def.hpp"
# undef CLASSNAME

#endif

#define CLASSNAME DJadaInnerGmresTest_25_1_2
#include "phist_gen_d.h"
#include "JadaInnerGmresTest_def.hpp"
#undef CLASSNAME

#define CLASSNAME ZJadaInnerGmresTest_25_1_2
#include "phist_gen_z.h"
#include "JadaInnerGmresTest_def.hpp"
#undef CLASSNAME


#undef MATNAME
#define MATNAME "jadaTestMat"

#undef _N_
#undef _NV_
#undef _NVP_
#undef _MAXBAS_

#define _N_ 512
#define _NV_ 3
#define _NVP_ 7
#define _MAXBAS_ 40

#ifdef PHIST_HAVE_SP

# define CLASSNAME SJadaInnerGmresTest_512_3_7
# include "phist_gen_s.h"
# include "JadaInnerGmresTest_def.hpp"
# undef CLASSNAME

# define CLASSNAME CJadaInnerGmresTest_512_3_7
# include "phist_gen_c.h"
# include "JadaInnerGmresTest_def.hpp"
#undef CLASSNAME

#endif

#define CLASSNAME DJadaInnerGmresTest_512_3_7
#include "phist_gen_d.h"
#include "JadaInnerGmresTest_def.hpp"
#undef CLASSNAME

#define CLASSNAME ZJadaInnerGmresTest_512_3_7
#include "phist_gen_z.h"
#include "JadaInnerGmresTest_def.hpp"
