/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef PHIST_BLOCKEDGMRES_H
#define PHIST_BLOCKEDGMRES_H

#include "phist_config.h"

#ifndef DOXYGEN

/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "phist_operator.h"
#include "phist_enums.h"
#include "phist_typedefs.h"

#endif //DOXYGEN

#ifdef __cplusplus
extern "C" {
#endif

//! \defgroup blockedGMRES blocked GMRES solver for general linear systems
//! \ingroup linear_solvers
//!@{

#define PHIST_CLASSFILE_DEF "phist_blockedgmres_decl.h"
#include "phist_gen_all.h"

//!@}
#ifdef __cplusplus
}
#endif

#endif
