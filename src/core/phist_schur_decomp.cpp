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
#include "phist_schur_decomp.h"
#include "jada_helpers.hpp"
#include "phist_kernels.h"
#include "phist_lapack.h"
#include "phist_MemOwner.hpp"

#ifdef __cplusplus
extern "C" {
#endif

#define PHIST_CLASSFILE_DEF "phist_schur_decomp_def.hpp"
#include "phist_gen_all.h"

#ifdef __cplusplus
}
#endif
