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

#include "phist_mpi_kernels.h"
#include "phist_schur_decomp.h"
#include "../kernels/KernelTest.h"
#include "../kernels/TestWithType.h"
#include "../kernels/KernelTestWithSdMats.h"

using namespace testing;

#ifdef CLASSNAME
#undef CLASSNAME
#endif

#define _N_ 10

#ifdef PHIST_HAVE_SP

#define CLASSNAME STestSchurDecomp_10
#include "phist_gen_s.h"
#include "TestSchurDecomp_def.hpp"
#undef CLASSNAME

#ifdef PHIST_HAVE_CMPLX
#define CLASSNAME CTestSchurDecomp_10
#include "phist_gen_c.h"
#include "TestSchurDecomp_def.hpp"
#undef CLASSNAME
#endif
#endif

#define CLASSNAME DTestSchurDecomp_10
#include "phist_gen_d.h"
#include "TestSchurDecomp_def.hpp"
#undef CLASSNAME

#ifdef PHIST_HAVE_CMPLX
#define CLASSNAME ZTestSchurDecomp_10
#include "phist_gen_z.h"
#include "TestSchurDecomp_def.hpp"
#undef CLASSNAME
#endif

#undef _N_

#define _N_ 50

#ifdef PHIST_HAVE_SP
#define CLASSNAME STestSchurDecomp_50
#include "phist_gen_s.h"
#include "TestSchurDecomp_def.hpp"
#undef CLASSNAME

#ifdef PHIST_HAVE_CMPLX
#define CLASSNAME CTestSchurDecomp_50
#include "phist_gen_c.h"
#include "TestSchurDecomp_def.hpp"
#undef CLASSNAME
#endif

#endif

#define CLASSNAME DTestSchurDecomp_50
#include "phist_gen_d.h"
#include "TestSchurDecomp_def.hpp"
#undef CLASSNAME

#ifdef PHIST_HAVE_CMPLX
#define CLASSNAME ZTestSchurDecomp_50
#include "phist_gen_z.h"
#include "TestSchurDecomp_def.hpp"
#undef CLASSNAME
#endif
#undef _N_
