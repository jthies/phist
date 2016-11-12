#ifndef PHIST_KERNEL_TEST_WITH_SDMATS_H
#define PHIST_KERNEL_TEST_WITH_SDMATS_H

#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/phist_gtest.h"

#include "phist_typedefs.h"
#include "phist_kernels.h"
// for sync'ing sdMats
#include "phist_mpi_kernels.h"
#include "TestWithType.h"

#include "phist_memOwner.hpp"

using namespace testing;

/*! Test fixure. 
 *  _multipleDefinitionCounter is used to enforce multiple definitions of static class variables
 *  if the data is required multiple times!
 */
template<typename T, int _Nrows, int _Ncols, int _useViews=0, int _multipleDefinitionCounter=0>
class KernelTestWithSdMats:
        public virtual TestWithType<T>,
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
