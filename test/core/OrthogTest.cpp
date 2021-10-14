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

#include "phist_orthog.h"

#include "../kernels/KernelTestWithVectors.h"
#include "../kernels/KernelTestWithSdMats.h"

#if defined(PHIST_HAVE_BELOS)||defined(PHIST_HAVE_ANASAZI)
#include "phist_types.hpp"
#include "phist_core.hpp"
#include "phist_BelosMV.hpp"
#endif

#ifdef PHIST_HAVE_BELOS
#include "Belos_PhistAdapter.hpp"
#include "phist_BelosOperatorTraits.hpp"
#include "BelosICGSOrthoManager.hpp"
// reference ortho managers
#include "BelosIMGSOrthoManager.hpp"
#include "BelosDGKSOrthoManager.hpp"

// our own orthog-based adaptation:
#include "BelosICGSOrthoManager.hpp"
#define PHIST_CLASSFILE_DEF "phist_BelosICGSOrthoManager.hpp"
#include "phist_gen_all.h"
#endif

#ifdef PHIST_HAVE_ANASAZI
#include "Anasazi_PhistAdapter.hpp"
#include "phist_AnasaziOperatorTraits.hpp"

// we use this as reference ortho manager:
#include "AnasaziICGSOrthoManager.hpp"
// and this is our orthog-based adaptation of the SVQB ortho manager
#include "AnasaziSVQBOrthoManager.hpp"
#define PHIST_CLASSFILE_DEF "phist_AnasaziSVQBOrthoManager.hpp"
#include "phist_gen_all.h"
#endif

using namespace testing;

#define _BASENAME_ OrthogTest

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

