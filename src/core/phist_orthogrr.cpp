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
#include "phist_orthogrr.h"
#include "phist_normF.h"
#include "phist_sdFact.h"
#include "phist_MemOwner.hpp"

#ifndef PHIST_HIGH_PRECISION_KERNELS
/*#define ORTHOGRR_USE_SVQB 1*/
#endif
#define PHIST_CLASSFILE_DEF "phist_orthogrr_def.hpp"
#include "phist_gen_all.h"
