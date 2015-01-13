#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#include "kernels/phist_kernels.h"
#include "KernelTestWithVectors.h"
#include "../tools/MatrixIO.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

using namespace ::testing;

#define _N_ 25
#define _NV_ 1

#ifdef PHIST_HAVE_SP
# define CLASSNAME SCarpTest_25_1
# include "phist_gen_s.h"
# include "CarpTest_def.hpp"
# undef CLASSNAME

# define CLASSNAME CCarpTest_25_1
# include "phist_gen_c.h"
# include "CarpTest_def.hpp"
# undef CLASSNAME
#endif

#define CLASSNAME DCarpTest_25_1
#include "phist_gen_d.h"
#include "CarpTest_def.hpp"
#undef CLASSNAME

#define CLASSNAME ZCarpTest_25_1
#include "phist_gen_z.h"
#include "CarpTest_def.hpp"
#undef CLASSNAME

#undef _N_
#undef _NV_
#define _N_ 25
#define _NV_ 4

#ifdef PHIST_HAVE_SP

# define CLASSNAME SCarpTest_25_4
# include "phist_gen_s.h"
# include "CarpTest_def.hpp"
# undef CLASSNAME

# define CLASSNAME CCarpTest_25_4
# include "phist_gen_c.h"
# include "CarpTest_def.hpp"
#undef CLASSNAME

#endif

#define CLASSNAME DCarpTest_25_4
#include "phist_gen_d.h"
#include "CarpTest_def.hpp"
#undef CLASSNAME

#define CLASSNAME ZCarpTest_25_4
#include "phist_gen_z.h"
#include "CarpTest_def.hpp"

#undef CLASSNAME
#undef _N_
#undef _NV_
#define _N_ 25
#define _NV_ 7

#ifdef PHIST_HAVE_SP

# define CLASSNAME SCarpTest_25_7
# include "phist_gen_s.h"
# include "CarpTest_def.hpp"
# undef CLASSNAME

# define CLASSNAME CCarpTest_25_7
# include "phist_gen_c.h"
# include "CarpTest_def.hpp"
#undef CLASSNAME

#endif

#define CLASSNAME DCarpTest_25_7
#include "phist_gen_d.h"
#include "CarpTest_def.hpp"
#undef CLASSNAME

#define CLASSNAME ZCarpTest_25_7
#include "phist_gen_z.h"
#include "CarpTest_def.hpp"
