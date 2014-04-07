#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#include "phist_kernels.h"
#include "KernelTestWithVectors.h"
#include "KernelTestWithSdMats.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

using namespace testing;

#define _N_ 12
#define _NV_ 9
#ifdef CLASSNAME
#undef CLASSNAME
#endif

#ifdef PHIST_HAVE_SP
#define CLASSNAME SGemmTest_12_9
#include "phist_gen_s.h"
#include "GemmTest_def.hpp"
#undef CLASSNAME

#define CLASSNAME CGemmTest_12_9
#include "phist_gen_c.h"
#include "GemmTest_def.hpp"
#undef CLASSNAME

#endif


#define CLASSNAME DGemmTest_12_9
#include "phist_gen_d.h"
#include "GemmTest_def.hpp"
#undef CLASSNAME

#define CLASSNAME ZGemmTest_12_9
#include "phist_gen_z.h"
#include "GemmTest_def.hpp"
#undef _N_

#define _N_ 13
#undef _NV_
#define _NV_ 18
#undef CLASSNAME

#ifdef PHIST_HAVE_SP

#define CLASSNAME SGemmTest_13_18
#include "phist_gen_s.h"
#include "GemmTest_def.hpp"
#undef CLASSNAME

#define CLASSNAME CGemmTest_13_18
#include "phist_gen_c.h"
#include "GemmTest_def.hpp"
#undef CLASSNAME

#endif

#define CLASSNAME DGemmTest_13_18
#include "phist_gen_d.h"
#include "GemmTest_def.hpp"
#undef CLASSNAME

#define CLASSNAME ZGemmTest_13_18
#include "phist_gen_z.h"
#include "GemmTest_def.hpp"
