/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
/*! \file builtin/prec_kernels.cpp
 * wraps implementation of builtin kernel routines
 * \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
 * \author "Jonas Thies <Jonas.Thies@DLR.de>
 *
*/

#include "phist_config.h"

#include "phist_tools.h"
#include "phist_kernel_perfmodels.hpp"

#include "phist_typedefs.h"
#include "phist_MemOwner.hpp"
#include "phist_sdFact.h"

#include "phist_ScalarTraits.hpp"

#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_sdFact_kernels_def.hpp"
#include "phist_sdFact_prec_kernels_def.hpp"
#include "phist_sdFact_def.hpp"

# ifdef PHIST_HAVE_CMPLX
# include "phist_gen_c.h"
# include "phist_sdFact_kernels_def.hpp"
# include "phist_sdFact_prec_kernels_def.hpp"
# include "phist_sdFact_def.hpp"
# endif
#endif

#ifdef PHIST_HAVE_CMPLX
#include "phist_gen_z.h"
#include "phist_sdFact_kernels_def.hpp"
#include "phist_sdFact_prec_kernels_def.hpp"
#include "phist_sdFact_def.hpp"
#endif

#include "phist_gen_d.h"
#include "phist_sdFact_kernels_def.hpp"
#ifdef PHIST_HIGH_PRECISION_KERNELS
extern "C" {
// high precision variants only available in "D" case up to now
#include "DsdFact_prec_kernels.c"
}
#endif
#include "phist_sdFact_def.hpp"


