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
#include "phist_core.h"
#include "phist_ScalarTraits.hpp"
#include "phist_MemOwner.hpp"

#include "phist_subspacejada.h"

#include "phist_schur_decomp.h"
#include "phist_simple_arnoldi.h"
#include "phist_transform_searchspace.h"
#include "phist_jadaCorrectionSolver.h"

#include "phist_gen_d.h"
#include "phist_subspacejada_def.hpp"

