#include "phist_config.h"

#include "phist_tools.h"
#include "phist_kernels.h"
#include "gtest/phist_gtest.h"

#include "phist_blockedgmres.h"
#include "phist_jadaOp.h"
#include "phist_orthog.h"
#include "../kernels/KernelTestWithSparseMat.h"
#include "../kernels/KernelTestWithSdMats.h"
#include "../kernels/KernelTestWithVectors.h"

using namespace testing;

#ifdef CLASSNAME
#undef CLASSNAME
#endif

#ifdef MATNAME
#undef MATNAME
#endif

#define MATNAME MATNAME_sprandn

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
#define MATNAME MATNAME_jadaTestMat

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
