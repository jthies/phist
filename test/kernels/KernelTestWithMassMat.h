#ifndef PHIST_KERNEL_TEST_WITH_MASS_MATRIX_H
#define PHIST_KERNEL_TEST_WITH_MASS_MATRIX_H

#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "phist_kernels.h"
#include "gtest/phist_gtest.h"

#include "phist_typedefs.h"
#include "TestWithType.h"

#include "../tools/MatrixIO.h"

using namespace testing;

 /* for testing functionality that requires two matrices A and B, where B is Hermitian and positive definite,
    One can derive the test class from both KernelTestWithSparseMat and KernelTestWithMassMat. This class adds
    the B_ member, a TYPE(sparseMat_ptr) that represents some H.p.d. matrix of dimension _Nglob.
 */
template<typename T, phist_gidx _Nglob>
class KernelTestWithMassMat: public virtual TestWithType<T>,
                             public virtual KernelTest
  {

public:

  };

#ifdef PHIST_HAVE_SP

#include "phist_gen_s.h"
#include "KernelTestWithMassMat_def.h"

#include "phist_gen_c.h"
#include "KernelTestWithMassMat_def.h"

#endif

#include "phist_gen_d.h"
#include "KernelTestWithMassMat_def.h"

#include "phist_gen_z.h"
#include "KernelTestWithMassMat_def.h"

#endif
