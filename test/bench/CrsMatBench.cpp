
#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#include "kernels/phist_kernels.h"
#include "../kernels/KernelTestWithVectors.h"


#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

using namespace testing;

#define _N_ 1280
#define _NV_ 4

#define CLASSNAME C_BENCH_MHD1280B_4
#define MATNAME "mhd1280b"

#include "phist_gen_c.h"
#include "CrsMatBench_def.hpp"

#undef CLASSNAME
#define CLASSNAME Z_BENCH_MHD1280B_4

#include "phist_gen_z.h"
#include "CrsMatBench_def.hpp"

#undef _N_
#undef _NV_
#define _N_ 4800
#define _NV_ 4

#undef CLASSNAME
#undef MATNAME
#define CLASSNAME S_BENCH_MHD4800B_4
#define MATNAME "mhd4800b"

#include "phist_gen_s.h"
#include "CrsMatBench_def.hpp"

#undef CLASSNAME
#define CLASSNAME D_BENCH_MHD4800B_4

#include "phist_gen_d.h"
#include "CrsMatBench_def.hpp"

