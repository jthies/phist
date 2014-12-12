#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif


#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#include "../kernels/KernelTestWithVectors.h"
#include "../kernels/KernelTestWithSdMats.h"

#include "phist_kernels.h"
#include "phist_svqb.h"

using namespace testing;

#ifdef CLASSNAME
#undef CLASSNAME
#endif

#define _N_ 24
#define _NV_ 1

#ifdef PHIST_HAVE_SP
#define CLASSNAME SSVQB_Test_24_1
#include "phist_gen_s.h"
#include "SVQB_Test_def.hpp"

#undef CLASSNAME
#define CLASSNAME CSVQB_Test_24_1

#include "phist_gen_c.h"
#include "SVQB_Test_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DSVQB_Test_24_1

#include "phist_gen_d.h"
#include "SVQB_Test_def.hpp"


#undef CLASSNAME
#define CLASSNAME ZSVQB_Test_24_1

#include "phist_gen_z.h"
#include "SVQB_Test_def.hpp"

#undef _N_
#define _N_ 59
#undef _NV_
#define _NV_ 5

#ifdef PHIST_HAVE_SP
#undef CLASSNAME
#define CLASSNAME SSVQB_Test_59_5

#include "phist_gen_s.h"
#include "SVQB_Test_def.hpp"

#undef CLASSNAME
#define CLASSNAME CSVQB_Test_59_5

#include "phist_gen_c.h"
#include "SVQB_Test_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DSVQB_Test_59_5

#include "phist_gen_d.h"
#include "SVQB_Test_def.hpp"
#undef CLASSNAME
#define CLASSNAME ZSVQB_Test_59_5

#include "phist_gen_z.h"
#include "SVQB_Test_def.hpp"

#if 1
// let's try something bigger...
#undef _N_
#define _N_ 9999
#undef _NV_
#define _NV_ 65

#ifdef PHIST_HAVE_SP

#undef CLASSNAME
#define CLASSNAME SSVQB_Test_10k_65

#include "phist_gen_s.h"
#include "SVQB_Test_def.hpp"

#undef CLASSNAME
#define CLASSNAME CSVQB_Test_10k_65

#include "phist_gen_c.h"
#include "SVQB_Test_def.hpp"

#endif

#undef CLASSNAME
#define CLASSNAME DSVQB_Test_10k_65

#include "phist_gen_d.h"
#include "SVQB_Test_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZSVQB_Test_10k_65

#include "phist_gen_z.h"
#include "SVQB_Test_def.hpp"

#endif
