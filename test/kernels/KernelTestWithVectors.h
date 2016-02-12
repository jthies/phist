#ifndef PHIST_KERNEL_TEST_WITH_VECTORS_H
#define PHIST_KERNEL_TEST_WITH_VECTORS_H

#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"

#include "phist_typedefs.h"
#include "phist_kernels.h"
#include "KernelTestWithMap.h"
#include "TestWithType.h"

using namespace testing;

/*! Base class for tests using mvecs. The class is templated on 
 * the data type T,
 * the global number of rows (_Nglob),
 * the number of colums (_Nvec),
 * _useViews (default false): setup the owned mvecs as views of larger
 * _numberOfVectorsInitialized set to 1-3, makes vec1_, vec2_ and vec3_ available
 * _multipleDefinitionCounter: used to enforce multiple template instantiations of static class variables where needed
   mvec blocks.
 */
template<typename T, gidx_t _Nglob, int _Nvec, int _useViews=0, int _numberOfVectorsInitialized = 1,int _multipleDefinitionCounter=0>
class KernelTestWithVectors:
        public virtual KernelTestWithMap<_Nglob>,
        public virtual TestWithType<T>
  {

public:

  };

#ifdef PHIST_HAVE_SP

#include "phist_gen_s.h"
#include "KernelTestWithVectors_def.h"

#include "phist_gen_c.h"
#include "KernelTestWithVectors_def.h"

#endif

#include "phist_gen_d.h"
#include "KernelTestWithVectors_def.h"

#include "phist_gen_z.h"
#include "KernelTestWithVectors_def.h"

#endif
