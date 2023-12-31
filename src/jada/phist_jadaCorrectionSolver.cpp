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
#include "phist_core.h"
#include "phist_krylov.h"
#include "phist_ScalarTraits.hpp"
#include "phist_jadaCorrectionSolver.h"
#include "phist_jadaOp.h"
#include "phist_precon.h"

#include "phist_MemOwner.hpp"

#include <cstdlib>
#include <vector>

#define PHIST_CLASSFILE_DEF "phist_jadaCorrectionSolver_def.hpp"
#include "phist_gen_all.h"
