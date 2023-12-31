/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef PHIST_SDFACT_H
#define PHIST_SDFACT_H

#include "phist_config.h"

#ifndef DOXYGEN

#include "phist_void_aliases.h"

#endif //DOXYGEN

#ifdef __cplusplus
extern "C" {
#endif

//! \defgroup sdFact sdFact: various factorization methods for small dense matrices
//! \ingroup core
//!
//! These kernels are highly accurate but not necessarily efficient, If the
//! kernel lib supports PHIST_HIGH_PRECISION_KERNELS (e.g. the builtin kernels
//! with AVX2), the input *iflag=PHIST_ROBUST_REDUCTIONS enables the high precision
//! variants. Otherwise, standard precision is used.
//! 
//! \note you should always make sure that on hybrid CPU/GPU systems the sdMat values on the
//! host and device are synchronized by calling sdMat_from(to)_device before (after) using 
//! these functions. They all assume that the data obtained from sdMat_extract_view/error is
//! up-to-date with the device, and do not call to_device afterwards.
//!@{

#define PHIST_CLASSFILE_DEF "phist_sdFact_decl.h"
#include "phist_gen_all.h"

//!@}
#ifdef __cplusplus
} // extern "C"
#endif

#endif
