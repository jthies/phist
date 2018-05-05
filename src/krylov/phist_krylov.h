#ifndef PHIST_KRYLOV_H
#define PHIST_KRYLOV_H

#include "phist_config.h"

//! \defgroup krylov krylov: Krylov subspace methods
//!@{

/* blocked Krylov solvers used as correction solvers in jada */
#include "phist_blockedbicgstab.h"
#include "phist_blockedgmres.h"
#include "phist_blockedminres.h"
#include "phist_blockedqmr.h"

/* Krylov methods for linear systems - high level single-vector interface */
#include "phist_linear_solvers.h"

#include "phist_carp_cg.h"

#ifdef PHIST_HAVE_ANASAZI
#include "phist_anasazi.h"
#endif

#ifdef PHIST_HAVE_BELOS
#include "phist_belos.h"
#endif

/* simple iteration schemes to start the Jacobi-Davidson method */
#include "phist_simple_arnoldi.h"
#include "phist_simple_lanczos.h"

#include "phist_iter_op.h"
//!@}
#endif
