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

#include "KernelTestWithMap.h"

using namespace testing;

#define _N_ 1
#define CLASSNAME XMapTest_1
#include "MapTest_def.hpp"

#undef _N_
#define _N_ 8
#undef CLASSNAME
#define CLASSNAME XMapTest_8
#include "MapTest_def.hpp"
