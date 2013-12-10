
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#include "phist_kernels.h"
#include "KernelTest.h"

#include "KernelTestWithMap.h"

using namespace testing;

#define _N_ 1
#define CLASSNAME MapTest_1
#include "MapTest_def.hpp"

#undef _N_
#define _N_ 8
#undef CLASSNAME
#define CLASSNAME MapTest_8
#include "MapTest_def.hpp"
