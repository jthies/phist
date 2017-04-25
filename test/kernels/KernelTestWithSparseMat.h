/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef PHIST_KERNEL_TEST_WITH_MATRIX_H
#define PHIST_KERNEL_TEST_WITH_MATRIX_H

#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "phist_kernels.h"
#include "gtest/phist_gtest.h"

#include "phist_typedefs.h"
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
#define MATNAME_IDFUNC 9
#define MATNAME_hpd_tridiag 10
#define MATNAME_nhpd_tridiag 11
#define MATNAME_hid_tridiag 12
#define MATNAME_nhid_tridiag 13
#define MATNAME_hpd_tridiag_ainv 14
#define MATNAME_nhpd_tridiag_ainv 15
#define MATNAME_hid_tridiag_ainv 16
#define MATNAME_nhid_tridiag_ainv 17

inline bool MatNameEnumIsMatFunc(MATNAME_ENUM e){return (e==MATNAME_BENCH3D_8_A1)||(e==MATNAME_IDFUNC);}

const char* MatNameEnumToStr(MATNAME_ENUM);

/*! Base class for tests using sparse matrices. The class is templated on 
 * the data type T,
 * the global number of rows (_Nglob),
 * the matrix file/name (_MatName),
 * _multipleDefinitionCounter: used to enforce multiple template instantiations of static class variables where needed
   mvec blocks.
 */
template<typename T, phist_gidx _Nglob, phist_gidx _Mglob, MATNAME_ENUM _MatName, int _multipleDefinitionCounter=0>
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
