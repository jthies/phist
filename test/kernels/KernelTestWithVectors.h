/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef PHIST_KERNEL_TEST_WITH_VECTORS_H
#define PHIST_KERNEL_TEST_WITH_VECTORS_H

#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/phist_gtest.h"

#include "phist_tools.h"
#include "phist_kernels.h"
#include "KernelTestWithMap.h"
#include "TestWithType.h"

#ifdef PHIST_HIGH_PRECISION_KERNELS
#include "phist_prec_helpers.h"
#endif

using namespace testing;

/*! Base class for tests using mvecs. The class is templated on 
 * the data type T,
 * the global number of rows (_Nglob),
 * the number of colums (_Nvec),
 * _useViews (default false): setup the owned mvecs as views of larger
 * _numberOfVectorsInitialized set to 1-4, makes vec1_, vec2_, vec3_ and vec4_ available
 * _multipleDefinitionCounter: used to enforce multiple template instantiations of static class variables where needed
   mvec blocks.
 */
template<typename T, phist_gidx _Nglob, int _Nvec, int _useViews=0, int _numberOfVectorsInitialized = 1,int _multipleDefinitionCounter=0>
class KernelTestWithVectors:
        public virtual KernelTestWithMap<_Nglob>,
        public virtual TestWithType<T>
  {

public:

  };

#ifdef PHIST_HAVE_SP

#include "phist_gen_s.h"
#include "KernelTestWithVectors_def.h"

# ifdef PHIST_HAVE_CMPLX
# include "phist_gen_c.h"
# include "KernelTestWithVectors_def.h"
# endif
#endif

#include "phist_gen_d.h"
#include "KernelTestWithVectors_def.h"

#ifdef PHIST_HAVE_CMPLX
#include "phist_gen_z.h"
#include "KernelTestWithVectors_def.h"
#endif

#endif
