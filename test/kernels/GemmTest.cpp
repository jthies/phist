
#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#include "phist_kernels.h"
#include "KernelTestWithVectors.h"
#include "KernelTestWithSdMats.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

using namespace testing;

#define _N_ 6
#define _NV_ 4
#ifdef CLASSNAME
#undef CLASSNAME
#endif
#define CLASSNAME SGemmTest_6_4
#include "phist_gen_s.h"
#include "GemmTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME DGemmTest_6_4

#include "phist_gen_d.h"
#include "GemmTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CGemmTest_6_4

#include "phist_gen_c.h"
#include "GemmTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZGemmTest_6_4

#include "phist_gen_z.h"
#include "GemmTest_def.hpp"

#undef _N_
#define _N_ 6
#undef _NV_
#define _NV_ 8

#undef CLASSNAME
#define CLASSNAME SGemmTest_6_8

#include "phist_gen_s.h"
#include "GemmTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME DGemmTest_6_8

#include "phist_gen_d.h"
#include "GemmTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CGemmTest_6_8

#include "phist_gen_c.h"
#include "GemmTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZGemmTest_6_8

#include "phist_gen_z.h"
#include "GemmTest_def.hpp"
