/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef PHIST_SUBSPACEJADA_H
#define PHIST_SUBSPACEJADA_H

#ifndef DOXYGEN

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "phist_operator.h"
#include "phist_enums.h"
#include "phist_jadaOpts.h"

#endif /* DOXYGEN */

#ifdef __cplusplus
extern "C" {
#endif

//! \defgroup bjdqr subspacejada: Jacobi-Davidson QR method with blocking, locking and restart
//! \ingroup jada
//!@{

#define PHIST_CLASSFILE_DEF "phist_subspacejada_decl.h"
#include "phist_gen_all.h"

//!@}

#ifdef __cplusplus
}
#endif

#endif
