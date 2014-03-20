#ifndef PHIST_KERNEL_TEST_WITH_SDMATS_H
#define PHIST_KERNEL_TEST_WITH_SDMATS_H

#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#include "phist_typedefs.h"
#include "phist_kernels.h"
#include "KernelTestWithType.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

using namespace testing;

/*! Test fixure. */
template<typename T, int _Nrows, int _Ncols>
class KernelTestWithSdMats:
        public virtual KernelTestWithType<T>,
        public virtual KernelTest
  {

public:

  };

// we use a static member from this class for printing the matrix
#include "KernelTestWithVectors.h"

#ifdef PHIST_HAVE_SP

# include "phist_gen_s.h"
# include "KernelTestWithSdMats_def.h"

# include "phist_gen_c.h"
# include "KernelTestWithSdMats_def.h"

#endif

#include "phist_gen_d.h"
#include "KernelTestWithSdMats_def.h"

#include "phist_gen_z.h"
#include "KernelTestWithSdMats_def.h"

#endif
