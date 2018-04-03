/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_config.h"

#include "phist_tools.h"
#include "phist_kernels.h"

#include "gtest/phist_gtest.h"

#include "../tools/MatrixIO.h"

#include "phist_orthog.h"

#include "../kernels/KernelTestWithMassMat.h"
#include "../kernels/KernelTestWithVectors.h"
#include "../kernels/KernelTestWithSdMats.h"

#ifdef PHIST_HAVE_BELOS
#include "phist_types.hpp"
#include "phist_core.hpp"
#include "phist_BelosMV.hpp"
#include "Belos_PhistAdapter.hpp"
#include "phist_BelosOperatorTraits.hpp"
#include "BelosICGSOrthoManager.hpp"
#include "BelosIMGSOrthoManager.hpp"
#include "BelosDGKSOrthoManager.hpp"

//#define PHIST_CLASSFILE_DEF "phist_BelosICGSOrthoManager.hpp"
//#include "phist_gen_all.h"
#endif
/*
#ifdef PHIST_HAVE_ANASAZI
#include "phist_types.hpp"
#include "phist_core.hpp"
#include "phist_BelosMV.hpp"
#include "Belos_PhistAdapter.hpp"
#include "Anasazi_PhistAdapter.hpp"
#include "phist_AnasaziOperatorTraits.hpp"

#include "AnasaziSVQBOrthoManager.hpp"

//#define PHIST_CLASSFILE_DEF "phist_AnasaziSVQBOrthoManager.hpp"
//#include "phist_gen_all.h"
#endif
*/
using namespace testing;

#define CLASSFILE_DEF "OrthogTest_def.hpp"
#define _BASENAME_ OrthogTestWithB

#define ORTHOG_WITH_HPD_B

#define _N_ 17
#define _M_ 1
#define _K_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 17
#define _M_ 6
#define _K_ 1
#include "../phist_typed_test_gen.h"

// larger block size for W
#define _N_ 64
#define _M_ 5
#define _K_ 3
#include "../phist_typed_test_gen.h"

