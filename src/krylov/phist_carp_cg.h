/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef PHIST_CARP_CG_H
#define PHIST_CARP_CG_H

#include "phist_config.h"

#ifndef DOXYGEN

/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "phist_enums.h"
#include "phist_typedefs.h"

#endif //DOXYGEN

#ifdef __cplusplus
extern "C" {
#endif

//! \defgroup carp_cg CARP-CG row projection method for general (shifted) linear systems
//! \ingroup linear_solvers
//!
//! Blocked CARP-CG solver for general shifted matrices A-sigma[j]I.
//!
//! Where a different shift is allowed for each RHS (as is required in our
//! block JDQR method). A special feature of this implementation is that it
//! can handle complex shifts even if the kernel library doesn't offer complex
//! kernels. To implement this, we use some wrapper classes defined in phist_carp_cg_kernels_decl.hpp.
//!
//! Another feature introduced this way is using additional projections to `precondition' the iteration.
//! The linear system solved is then in fact
//!
//!     |A-sigma[j]I    V||x+i*xi  |   |b|
//!     | V'            0||x'+i*xi'| = |0|
//!
//! The algorithm is CGMN (CG on the minimum norm problem AA'x=b with SSOR pre-
//! conditioning, implemented following the work of Bjoerck and Elfving (1979)).
//! The parallelization of the Kaczmarz sweeps is left to the kernel library
//! (functions carp_setup, carp_sweep). 
//!
//!@{

#define PHIST_CLASSFILE_DEF "phist_carp_cg_decl.h"
#include "phist_gen_all.h"

//!@}
#ifdef __cplusplus
}
#endif

#endif
