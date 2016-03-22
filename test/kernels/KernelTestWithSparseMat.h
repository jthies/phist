#ifndef PHIST_KERNEL_TEST_WITH_MATRIX_H
#define PHIST_KERNEL_TEST_WITH_MATRIX_H

#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/phist_gtest.h"

#include "phist_typedefs.h"
#include "phist_kernels.h"
#include "KernelTestWithMap.h"
#include "TestWithType.h"

#include "../tools/MatrixIO.h"

using namespace testing;

// available matrices
// need to be preprocessor definitions to allow "#if MATNAME == MATNAME_speye" style comparisons
#define MATNAME_ENUM int
#define MATNAME_spzero 0
#define MATNAME_speye 1
#define MATNAME_sprandn 2
#define MATNAME_sprandn_nodiag 3
#define MATNAME_spshift 4
#define MATNAME_jadaTestMat 5
#define MATNAME_symmMat 6
#define MATNAME_sprandsym 7
#define MATNAME_BENCH3D_8_A1 8
const char* MatNameEnumToStr(MATNAME_ENUM);

/*! Base class for tests using sparse matrices. The class is templated on 
 * the data type T,
 * the global number of rows (_Nglob),
 * the matrix file/name (_MatName),
 * _multipleDefinitionCounter: used to enforce multiple template instantiations of static class variables where needed
   mvec blocks.
 */
template<typename T, phist_gidx _Nglob, MATNAME_ENUM _MatName, int _multipleDefinitionCounter=0>
class KernelTestWithSparseMat:
        public virtual KernelTestWithMap<_Nglob>,
        public virtual TestWithType<T>
  {

public:

  };

#ifdef PHIST_HAVE_SP

#include "phist_gen_s.h"
#include "KernelTestWithSparseMat_def.h"

#include "phist_gen_c.h"
#include "KernelTestWithSparseMat_def.h"

#endif

#include "phist_gen_d.h"
#include "KernelTestWithSparseMat_def.h"

#include "phist_gen_z.h"
#include "KernelTestWithSparseMat_def.h"

#endif
