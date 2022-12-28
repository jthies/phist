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
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>

#include "phist_blockedidrs.h"
#include "phist_orthog.h"
#include "phist_MemOwner.hpp"

#define PHIST_CLASSFILE_DEF "phist_blockedidrs_def.hpp"
#include "phist_gen_all.h"
