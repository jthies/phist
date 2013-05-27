
#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "kernels/phist_kernels.h"
#include "KernelTestWithVectors.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

using namespace testing;

#define _N_ 25
#define _NV_ 1
#define CLASSNAME _TESTNAME2_(CrsMatTest,25,1)

#include "phist_gen_s.h"
#define _TPC_ "s_"
#include "CrsMatTest_def.hpp"

#include "phist_gen_d.h"
#undef _TPC_
#define _TPC_ "d_"
#include "CrsMatTest_def.hpp"

#include "phist_gen_c.h"
#undef _TPC_
#define _TPC_ "c_"
#include "CrsMatTest_def.hpp"

#include "phist_gen_z.h"
#undef _TPC_
#define _TPC_ "z_"
#include "CrsMatTest_def.hpp"

#undef _N_
#undef _NV_
#undef CLASSNAME
#define _N_ 25
#define _NV_ 4
#define CLASSNAME _TESTNAME2_(CrsMatTest,25,4)

#include "phist_gen_s.h"
#undef _TPC_
#define _TPC_ "s_"
#include "CrsMatTest_def.hpp"

#include "phist_gen_d.h"
#undef _TPC_
#define _TPC_ "d_"
#include "CrsMatTest_def.hpp"

#include "phist_gen_c.h"
#undef _TPC_
#define _TPC_ "c_"
#include "CrsMatTest_def.hpp"

#include "phist_gen_z.h"
#undef _TPC_
#define _TPC_ "z_"
#include "CrsMatTest_def.hpp"
