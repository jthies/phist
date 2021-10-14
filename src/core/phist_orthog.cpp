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
#include "phist_MemOwner.hpp"
#include "phist_kernels.h"
#include "phist_orthog.h"

// fallback routine: if the kernel library does not implement mvec_QR,
// we simply switch to our own implementation of Cholesky-QR.
#include "phist_chol_QR.h"
// under the hood we use this implementation, which exploits fused kernels and
// high precision operations if available.
#include "phist_orthogrrfused.h"
#include "phist_orthogrr.h"

#define PHIST_CLASSFILE_DEF "phist_orthog_def.hpp"
#include "phist_gen_all.h"
