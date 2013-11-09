
#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#include "kernels/phist_kernels.h"
#include "../kernels/KernelTestWithVectors.h"


#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

using namespace testing;

#define MATNAME "mhd1280b"
#define _N_ 1280

#define _NV_ 1
#define CLASSNAME C_BENCH_MHD1280B_1
#include "phist_gen_c.h"
#include "CrsMatBench_def.hpp"
#undef CLASSNAME
#undef _NV_

#define _NV_ 2
#define CLASSNAME C_BENCH_MHD1280B_2
#include "phist_gen_c.h"
#include "CrsMatBench_def.hpp"
#undef CLASSNAME
#undef _NV_

#define _NV_ 4
#define CLASSNAME C_BENCH_MHD1280B_4
#include "phist_gen_c.h"
#include "CrsMatBench_def.hpp"
#undef CLASSNAME
#undef _NV_

#define _NV_ 8
#define CLASSNAME C_BENCH_MHD1280B_8
#include "phist_gen_c.h"
#include "CrsMatBench_def.hpp"
#undef CLASSNAME
#undef _NV_

#define _NV_ 16
#define CLASSNAME C_BENCH_MHD1280B_16
#include "phist_gen_c.h"
#include "CrsMatBench_def.hpp"
#undef CLASSNAME
#undef _NV_


#define _NV_ 1
#define CLASSNAME Z_BENCH_MHD1280B_1
#include "phist_gen_z.h"
#include "CrsMatBench_def.hpp"
#undef CLASSNAME
#undef _NV_

#define _NV_ 2
#define CLASSNAME Z_BENCH_MHD1280B_2
#include "phist_gen_z.h"
#include "CrsMatBench_def.hpp"
#undef CLASSNAME
#undef _NV_

#define _NV_ 4
#define CLASSNAME Z_BENCH_MHD1280B_4
#include "phist_gen_z.h"
#include "CrsMatBench_def.hpp"
#undef CLASSNAME
#undef _NV_

#define _NV_ 8
#define CLASSNAME Z_BENCH_MHD1280B_8
#include "phist_gen_z.h"
#include "CrsMatBench_def.hpp"
#undef CLASSNAME
#undef _NV_

#define _NV_ 16
#define CLASSNAME Z_BENCH_MHD1280B_16
#include "phist_gen_z.h"
#include "CrsMatBench_def.hpp"
#undef CLASSNAME
#undef _NV_


#undef _N_
#undef MATNAME



#define MATNAME "mhd4800b"
#define _N_ 4800


#define _NV_ 1
#define CLASSNAME S_BENCH_MHD4800B_1
#include "phist_gen_s.h"
#include "CrsMatBench_def.hpp"
#undef CLASSNAME
#undef _NV_

#define _NV_ 2
#define CLASSNAME S_BENCH_MHD4800B_2
#include "phist_gen_s.h"
#include "CrsMatBench_def.hpp"
#undef CLASSNAME
#undef _NV_

#define _NV_ 4
#define CLASSNAME S_BENCH_MHD4800B_4
#include "phist_gen_s.h"
#include "CrsMatBench_def.hpp"
#undef CLASSNAME
#undef _NV_

#define _NV_ 8
#define CLASSNAME S_BENCH_MHD4800B_8
#include "phist_gen_s.h"
#include "CrsMatBench_def.hpp"
#undef CLASSNAME
#undef _NV_

#define _NV_ 16
#define CLASSNAME S_BENCH_MHD4800B_16
#include "phist_gen_s.h"
#include "CrsMatBench_def.hpp"
#undef CLASSNAME
#undef _NV_


#define _NV_ 1
#define CLASSNAME D_BENCH_MHD4800B_1
#include "phist_gen_d.h"
#include "CrsMatBench_def.hpp"
#undef CLASSNAME
#undef _NV_

#define _NV_ 2
#define CLASSNAME D_BENCH_MHD4800B_2
#include "phist_gen_d.h"
#include "CrsMatBench_def.hpp"
#undef CLASSNAME
#undef _NV_

#define _NV_ 4
#define CLASSNAME D_BENCH_MHD4800B_4
#include "phist_gen_d.h"
#include "CrsMatBench_def.hpp"
#undef CLASSNAME
#undef _NV_

#define _NV_ 8
#define CLASSNAME D_BENCH_MHD4800B_8
#include "phist_gen_d.h"
#include "CrsMatBench_def.hpp"
#undef CLASSNAME
#undef _NV_

#define _NV_ 16
#define CLASSNAME D_BENCH_MHD4800B_16
#include "phist_gen_d.h"
#include "CrsMatBench_def.hpp"
#undef CLASSNAME
#undef _NV_


#undef _N_
#undef MATNAME
