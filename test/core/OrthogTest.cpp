#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#include "phist_orthog.h"
#include "phist_kernels.h"

#include "../kernels/KernelTest.h"
#include "../kernels/KernelTestWithMap.h"
#include "../kernels/KernelTestWithType.h"
#include "../kernels/KernelTestWithVectors.h"
#include "../kernels/KernelTestWithSdMats.h"


#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

using namespace testing;

#ifdef CLASSNAME
#undef CLASSNAME
#endif

#define _N_ 17
#define _M_ 1
#define _K_ 1

#ifdef PHIST_HAVE_SP

# define CLASSNAME SOrthogTest_17_1_1
# include "phist_gen_s.h"
# include "OrthogTest_def.hpp"
# undef CLASSNAME

# define CLASSNAME COrthogTest_17_1_1
# include "phist_gen_c.h"
# include "OrthogTest_def.hpp"
# undef CLASSNAME
#endif

#define CLASSNAME DOrthogTest_17_1_1
#include "phist_gen_d.h"
#include "OrthogTest_def.hpp"
#undef CLASSNAME

#define CLASSNAME ZOrthogTest_17_1_1
#include "phist_gen_z.h"
#include "OrthogTest_def.hpp"
#undef CLASSNAME


#undef _M_
#define _M_ 6

#ifdef PHIST_HAVE_SP

# define CLASSNAME SOrthogTest_17_6_1
# include "phist_gen_s.h"
# include "OrthogTest_def.hpp"
# undef CLASSNAME

# define CLASSNAME COrthogTest_17_6_1
# include "phist_gen_c.h"
# include "OrthogTest_def.hpp"
# undef CLASSNAME
#endif

#define CLASSNAME DOrthogTest_17_6_1
#include "phist_gen_d.h"
#include "OrthogTest_def.hpp"
#undef CLASSNAME

#define CLASSNAME ZOrthogTest_17_6_1
#include "phist_gen_z.h"
#include "OrthogTest_def.hpp"
#undef CLASSNAME

// larger block size for W
#undef _N_
#define _N_ 64
#undef _M_
#define _M_ 5
#undef _K_
#define _K_ 3

#ifdef PHIST_HAVE_SP

# define CLASSNAME SOrthogTest_64_5_3
# include "phist_gen_s.h"
# include "OrthogTest_def.hpp"
# undef CLASSNAME

# define CLASSNAME COrthogTest_64_5_3
# include "phist_gen_c.h"
# include "OrthogTest_def.hpp"
# undef CLASSNAME

#endif

#define CLASSNAME DOrthogTest_64_5_3
#include "phist_gen_d.h"
#include "OrthogTest_def.hpp"
#undef CLASSNAME

#define CLASSNAME ZOrthogTest_64_5_3
#include "phist_gen_z.h"
#include "OrthogTest_def.hpp"
#undef CLASSNAME

#if 0
// do at least one larger test case
#undef _N_
#define _N_ 10000
#undef _M_
#define _M_ 40
#undef _K_
#define _K_ 8

#ifdef PHIST_HAVE_SP

# define CLASSNAME SOrthogTest_10k_40_8
# include "phist_gen_s.h"
# include "OrthogTest_def.hpp"
# undef CLASSNAME

# define CLASSNAME COrthogTest_10k_40_8
# include "phist_gen_c.h"
# include "OrthogTest_def.hpp"
# undef CLASSNAME

#endif

#define CLASSNAME DOrthogTest_10k_40_8

#include "phist_gen_d.h"
#include "OrthogTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZOrthogTest_10k_40_8

#include "phist_gen_z.h"
#include "OrthogTest_def.hpp"

#endif
