#include "phist_config.h"

#include "phist_tools.h"
#include "phist_kernels.h"

#include "gtest/phist_gtest.h"

#include "../kernels/KernelTestWithVectors.h"


#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

using namespace testing;

// This file is not currently used for serious testing
// but rather a playground how to integrate automatic
// performance testing. The idea up to now is to simply
// run some kernel (spMVM here) and monitor the time the
// test takes on something like a Jenkins server.

#define MATNAME "mhd1280b"
#define _N_ 1280

#define _NV_ 1


#ifdef PHIST_HAVE_SP
/*
# define CLASSNAME C_BENCH_MHD1280B_1
# include "phist_gen_c.h"
# include "CrsMatBench_def.hpp"
# undef CLASSNAME

# undef _NV_
# define _NV_ 2

# define CLASSNAME C_BENCH_MHD1280B_2
# include "CrsMatBench_def.hpp"
# undef CLASSNAME

# undef _NV_
# define _NV_ 4

# define CLASSNAME C_BENCH_MHD1280B_4
# include "CrsMatBench_def.hpp"
# undef CLASSNAME

# undef _NV_
# define _NV_ 8
# define CLASSNAME C_BENCH_MHD1280B_8
# include "CrsMatBench_def.hpp"
# undef CLASSNAME
*/
#endif
/*
#undef _NV_
#define _NV_ 1
#undef CLASSNAME
#define CLASSNAME Z_BENCH_MHD1280B_1
#include "phist_gen_z.h"
#include "CrsMatBench_def.hpp"
#undef CLASSNAME
#undef _NV_

#define _NV_ 2
#define CLASSNAME Z_BENCH_MHD1280B_2
#include "CrsMatBench_def.hpp"
#undef CLASSNAME
#undef _NV_

#define _NV_ 4
#define CLASSNAME Z_BENCH_MHD1280B_4
#include "CrsMatBench_def.hpp"
#undef CLASSNAME
#undef _NV_

#define _NV_ 8
#define CLASSNAME Z_BENCH_MHD1280B_8
#include "CrsMatBench_def.hpp"
#undef CLASSNAME

#undef _NV_
#undef _N_
#undef MATNAME



#define MATNAME "mhd4800b"
#define _N_ 4800

#ifdef PHIST_HAVE_SP

#define _NV_ 1
#ifdef PHIST_HAVE_SP
# define CLASSNAME S_BENCH_MHD4800B_1
# include "phist_gen_s.h"
# include "CrsMatBench_def.hpp"
# undef CLASSNAME

# undef _NV_
# define _NV_ 2

# define CLASSNAME S_BENCH_MHD4800B_2
# include "CrsMatBench_def.hpp"
# undef CLASSNAME

# undef _NV_
# define _NV_ 4
# define CLASSNAME S_BENCH_MHD4800B_4
# include "CrsMatBench_def.hpp"
# undef CLASSNAME

# undef _NV_
# define _NV_ 8
# define CLASSNAME S_BENCH_MHD4800B_8
# include "CrsMatBench_def.hpp"
# undef CLASSNAME
# undef _NV_
# undef CLASSNAME

#endif



#undef _NV_
#define _NV_ 1

#define CLASSNAME D_BENCH_MHD4800B_1
#include "phist_gen_d.h"
#include "CrsMatBench_def.hpp"
#undef CLASSNAME
#undef _NV_

#define _NV_ 2
#define CLASSNAME D_BENCH_MHD4800B_2
#include "CrsMatBench_def.hpp"
#undef CLASSNAME
#undef _NV_

#define _NV_ 4
#define CLASSNAME D_BENCH_MHD4800B_4
#include "CrsMatBench_def.hpp"
#undef CLASSNAME
#undef _NV_

#define _NV_ 8
#define CLASSNAME D_BENCH_MHD4800B_8
#include "CrsMatBench_def.hpp"
#undef CLASSNAME
#undef _NV_
*/
