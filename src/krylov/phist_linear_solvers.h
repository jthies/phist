/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef PHIST_LINEAR_SOLVERS_H
#define PHIST_LINEAR_SOLVERS_H

#include "phist_config.h"

#ifndef DOXYGEN

#include "phist_operator.h"
#include "phist_enums.h"
#include "phist_typedefs.h"

#endif //DOXYGEN

#ifdef __cplusplus
extern "C" {
#endif

//! \defgroup linear_solvers linear_solvers: Iterative methods for linear systems
//! \ingroup krylov
//!@{

#define PHIST_CLASSFILE_DEF "phist_linear_solvers_decl.h"
#include "phist_gen_all.h"
//!@}
#ifdef __cplusplus
}
#endif

#endif
