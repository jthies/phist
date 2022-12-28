/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_config.h"

#include "phist_tools.h"
#include "phist_kernels.h"
#include "gtest/phist_gtest.h"

#include "phist_precon.h"
#include "phist_krylov.h"
#include "phist_jadaOp.h"
#include "phist_subspacejada.h"
#include "phist_jadaCorrectionSolver.h"
#include "phist_orthog.h"
#include "../kernels/KernelTestWithSparseMat.h"
#include "../kernels/KernelTestWithSdMats.h"
#include "../kernels/KernelTestWithVectors.h"
#include "../jada/JadaTestWithOpts.hpp"

using namespace testing;

#define CLASSFILE_DEF "TestPreconOperators_def.hpp"

#define _WHICH_ phist_SR

#define MATNAME MATNAME_lapl_tridiag
#define PRECNAME MATNAME_lapl_tridiag_ainv
// use GMRES here, too because the left preconditioning would destroy
// the symmetry in the current implementation
//#define _SOLVTYPE_ phist_MINRES
#define _SOLVTYPE_ phist_GMRES
#define _MAXBAS_ 150
// number of inner Krylov iterations inside subspacejada
#define _MAXINNER_ 10

// matrix size
#define _N_ 300

// block size in jada/blockedGMRES _NV_
#define _M_ 1
// number of projection vectors in q (for testing operators)
#define _K_ 1
#define _BASENAME_ TestPreconOperators_lapl_tridiag
#include "../phist_typed_test_gen.h"

#define _N_ 300
#define _M_ 2
#define _K_ 6
#undef _BASENAME_
#define _BASENAME_ TestPreconOperators_lapl_tridiag
#include "../phist_typed_test_gen.h"

#undef MATNAME
#undef PRECNAME
#undef _BASENAME_
#undef _MAXBAS_

#define _N_ 300
#define _M_ 2
#define _K_ 1
#define MATNAME MATNAME_nhpd_tridiag
#define PRECNAME MATNAME_nhpd_tridiag_ainv
#define _SOLVTYPE_ phist_GMRES
#define _MAXBAS_ 300
#define _BASENAME_ TestPreconOperators_nhpd_tridiag
#include "../phist_typed_test_gen.h"

#define _N_ 300
#define _M_ 4
#define _K_ 6
#undef _BASENAME_
#define _BASENAME_ TestPreconOperators_nhpd_tridiag
#include "../phist_typed_test_gen.h"

#undef MATNAME
#undef PRECNAME
#undef _BASENAME_

#define _N_ 300
#define MATNAME MATNAME_hid_tridiag
#define PRECNAME MATNAME_lapl_tridiag_ainv
#define _BASENAME_ TestPreconOperators_hid_tridiag_with_hpd_precon
#define _SOLVTYPE_ phist_GMRES
#define _MAXBAS_ 300

#define _N_ 300
#define _M_ 1
#define _K_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 300
#define _M_ 2
#define _K_ 6
#include "../phist_typed_test_gen.h"


#undef PRECNAME
#undef _BASENAME_

#define PRECNAME MATNAME_hid_tridiag_ainv
#define _BASENAME_ TestPreconOperators_hid_tridiag_with_hid_precon

#define _N_ 300
#define _M_ 1
#define _K_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 300
#define _M_ 2
#define _K_ 6
#include "../phist_typed_test_gen.h"

