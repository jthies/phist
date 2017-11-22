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
#include <vector>
#ifdef PHIST_HAVE_TEUCHOS
#include <Teuchos_RCP.hpp>
#endif

#include "phist_blockedgmres.h"
#include "phist_orthog.h"

#define PHIST_CLASSFILE_DEF "phist_blockedgmres_def.hpp"
#include "phist_gen_all.h"
