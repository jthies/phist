/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
/*! \file builtin/prec_kernels.cpp
 * wraps implementation of builtin kernel routines
 * \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
 * \author "Jonas Thies <j.thies@tudelft.nl>
 *
*/

#include "phist_config.h"

#include "phist_tools.h"
#include "phist_kernel_perfmodels.hpp"

#include "phist_typedefs.h"
#include "phist_MemOwner.hpp"
#include "phist_sdFact.h"

#include "phist_ScalarTraits.hpp"

#define PHIST_CLASSFILE_DEF "phist_sdFact_kernels_def.hpp"
#include "phist_gen_all.h"
#undef PHIST_CLASSFILE_DEF

#define PHIST_CLASSFILE_DEF "phist_sdFact_prec_kernels_def.hpp"
#include "phist_gen_all.h"
#undef PHIST_CLASSFILE_DEF

#define PHIST_CLASSFILE_DEF "phist_sdFact_def.hpp"
#include "phist_gen_all.h"
#undef PHIST_CLASSFILE_DEF


#ifdef PHIST_HIGH_PRECISION_KERNELS
#include "phist_gen_d.h"
extern "C" {
// high precision variants only available in "D" case up to now
#include "DsdFact_prec_kernels.c"
}
#include "phist_gen_clean.h"
#endif


