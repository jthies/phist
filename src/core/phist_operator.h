/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef PHIST_OPERATOR_H
#define PHIST_OPERATOR_H

#include "phist_config.h"

/*! \file phist_operator.h \brief linear operator interface */

#ifndef DOXYGEN

#include "phist_void_aliases.h"
#include "phist_typedefs.h"

#endif //DOXYGEN

#ifdef __cplusplus
extern "C" {
#endif

//! \defgroup linearOp linearOp: operator interface
//! \ingroup core
//!@{

#define PHIST_CLASSFILE_DEF "phist_operator_decl.h"
#include "phist_gen_all.h"

//!@}

#ifdef __cplusplus
}
#endif

#endif
