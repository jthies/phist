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
#include "phist_orthogrr.h"
#include "phist_normF.h"
#include "phist_sdFact.h"
#include "phist_MemOwner.hpp"

// there are in principle two options for the block normalization: SVQB and Cholesky-QR.
// The relation between the input and output arguments V,W,Q,R1,R2 that we require, however,
// can't be fulfilled with SVQB (I think) in the singular case, so we use Cholesky-QR for now.
//#define ORTHOGRR_USE_SVQB 1
#define PHIST_CLASSFILE_DEF "phist_orthogrr_def.hpp"
#include "phist_gen_all.h"
