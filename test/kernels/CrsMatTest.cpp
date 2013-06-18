
#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#include "kernels/phist_kernels.h"
#include "KernelTestWithVectors.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

using namespace testing;

#define _N_ 25
#define _NV_ 1
#ifdef CLASSNAME
#undef CLASSNAME
#endif
#define CLASSNAME SCrsMatTest_25_1
#define _TPC_ "s_"

#include "phist_gen_s.h"
#include "CrsMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME DCrsMatTest_25_1
#undef _TPC_
#define _TPC_ "d_"

#include "phist_gen_d.h"
#include "CrsMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CCrsMatTest_25_1
#undef _TPC_
#define _TPC_ "c_"

#include "phist_gen_c.h"
#include "CrsMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZCrsMatTest_25_1
#undef _TPC_
#define _TPC_ "z_"
#include "phist_gen_z.h"
#include "CrsMatTest_def.hpp"

#undef _N_
#undef _NV_
#undef CLASSNAME
#define _N_ 25
#define _NV_ 4

#undef CLASSNAME
#define CLASSNAME SCrsMatTest_25_4
#undef _TPC_
#define _TPC_ "s_"

#include "phist_gen_s.h"
#include "CrsMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME DCrsMatTest_25_4
#undef _TPC_
#define _TPC_ "d_"

#include "phist_gen_d.h"
#include "CrsMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CCrsMatTest_25_4
#undef _TPC_
#define _TPC_ "c_"

#include "phist_gen_c.h"
#include "CrsMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZCrsMatTest_25_4
#undef _TPC_
#define _TPC_ "z_"

#include "phist_gen_z.h"
#include "CrsMatTest_def.hpp"
