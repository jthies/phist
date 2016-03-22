#ifndef KERNEL_TEST_WITHTYPEH
#define KERNEL_TEST_WITHTYPEH

#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/phist_gtest.h"
#include "KernelTest.h"
#include "phist_ScalarTraits.hpp"
#include "phist_random.h"
#include <cstdlib>
#include <limits>

/** 
 */
template<typename _ST>
class TestWithType : public virtual testing::Test
{
public:

//! set to true if the kernel lib supports the type (S/D/C/Z), false otherwise
static bool typeImplemented_;

void SetUpTestCase() {
typeImplemented_=false;
}

};

#include "phist_kernels.h"

#ifdef PHIST_HAVE_SP

# include "phist_gen_s.h"
# include "TestWithType_def.h"

# include "phist_gen_c.h"
# include "TestWithType_def.h"

#endif

#include "phist_gen_d.h"
#include "TestWithType_def.h"

#include "phist_gen_z.h"
#include "TestWithType_def.h"

#endif
