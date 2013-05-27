
#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "kernels/phist_kernels.h"
#include "KernelTestWithVectors.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

using namespace testing;

#define _N_ 9
#define _NV_ 1
#define CLASSNAME _TESTNAME2_(MvecTest,9,1)
#include "phist_gen_s.h"
#include "MvecTest_def.hpp"

#include "phist_gen_d.h"
#include "MvecTest_def.hpp"

#include "phist_gen_c.h"
#include "MvecTest_def.hpp"

#include "phist_gen_z.h"
#include "MvecTest_def.hpp"

#undef _N_
#define _N_ 6
#undef _NV_
#define _NV_ 3
#undef CLASSNAME
#define CLASSNAME _TESTNAME2_(MvecTest,6,3)

#include "phist_gen_s.h"
#include "MvecTest_def.hpp"

#include "phist_gen_d.h"
#include "MvecTest_def.hpp"

#include "phist_gen_c.h"
#include "MvecTest_def.hpp"

#include "phist_gen_z.h"
#include "MvecTest_def.hpp"
