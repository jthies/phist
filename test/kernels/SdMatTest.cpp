
#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "kernels/phist_kernels.h"
#include "KernelTestWithVectors.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

using namespace testing;

#define _NROWS_ 5
#define _NCOLS_ 5
#define CLASSNAME SSdMatTest_5_5
#include "phist_gen_s.h"
#include "SdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME DSdMatTest_5_5
#include "phist_gen_d.h"
#include "SdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CSdMatTest_5_5
#include "phist_gen_c.h"
#include "SdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZSdMatTest_5_5
#include "phist_gen_z.h"
#include "SdMatTest_def.hpp"

#undef _NROWS_
#define _NROWS_ 6
#undef _NCOLS_
#define _NCOLS_ 8

#undef CLASSNAME
#define CLASSNAME SSdMatTest_6_8
#include "phist_gen_s.h"
#include "SdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME DSdMatTest_6_8
#include "phist_gen_d.h"
#include "SdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME CSdMatTest_6_8
#include "phist_gen_c.h"
#include "SdMatTest_def.hpp"

#undef CLASSNAME
#define CLASSNAME ZSdMatTest_6_8
#include "phist_gen_z.h"
#include "SdMatTest_def.hpp"
