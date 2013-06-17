
#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "phist_kernels.h"
#include "KernelTestWithVectors.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

using namespace testing;

#define _N_ 9
#define _NV_ 1
#ifdef CLASSNAME
#undef CLASSNAME
#endif
#define CLASSNAME SMvecTest_9_1
#include "phist_gen_s.h"
#include "MvecTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME DMvecTest_9_1

#include "phist_gen_d.h"
#include "MvecTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecTest_9_1

#include "phist_gen_c.h"
#include "MvecTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecTest_9_1

#include "phist_gen_z.h"
#include "MvecTest_def.hpp"

#undef _N_
#define _N_ 6
#undef _NV_
#define _NV_ 3

#undef CLASSNAME
#define CLASSNAME SMvecTest_6_3

#include "phist_gen_s.h"
#include "MvecTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME DMvecTest_6_3

#include "phist_gen_d.h"
#include "MvecTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CMvecTest_6_3

#include "phist_gen_c.h"
#include "MvecTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZMvecTest_6_3

#include "phist_gen_z.h"
#include "MvecTest_def.hpp"
