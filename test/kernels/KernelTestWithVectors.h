#ifndef PHIST_KERNEL_TEST_WITH_VECTORS_H
#define PHIST_KERNEL_TEST_WITH_VECTORS_H

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "kernels/phist_kernels.h"
#include "KernelTestWithMap.h"
#include "KernelTestWithType.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

using namespace testing;

/*! Test fixure. */
template<typename T, int _Nglob, int _Nvec>
class KernelTestWithVectors:
        public virtual KernelTestWithMap<_Nglob>,
        public virtual KernelTestWithType<T>
  {

public:
  virtual void SetUp(){}
  virtual void TearDown(){}
  int nvec_, lda_, stride_;

  };

#include "phist_gen_s.h"
#include "KernelTestWithVectors_def.h"

#include "phist_gen_d.h"
#include "KernelTestWithVectors_def.h"

#include "phist_gen_c.h"
#include "KernelTestWithVectors_def.h"

#include "phist_gen_z.h"
#include "KernelTestWithVectors_def.h"

#endif
