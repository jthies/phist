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
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "phist_simple_lanczos_ft.h"

#ifndef GHOST_CP
#define GHOST_CP
#endif
#ifndef PHIST_CP
#define PHIST_CP
#endif

#include "checkpoint.hpp"
#include "cpTypes/cpPhistMvec/cpPhistMvec.h"
#include "cpOptions.h"



#define PHIST_CLASSFILE_DEF "phist_simple_lanczos_ft_def.hpp"
#include "phist_gen_all.h"
