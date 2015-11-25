#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <cstdio>
#include <cstdlib>

#include <iostream>
#include "phist_macros.h"
#include "../phist_kernels.h"
#include "phist_kernel_perfmodels.hpp"

#include "phist_typedefs.h"
#include "typedefs.hpp"
#include "phist_ScalarTraits.hpp"

#include "phist_ghost_internal.h"

// these are from Trilinos, we need them to interface
// the TSQR library for orthogonalizing tall skinny matrices.
#ifdef PHIST_HAVE_BELOS
#include "phist_trilinos_macros.h"
#include "Ghost_TsqrAdaptor.hpp"
#include "Belos_GhostAdapter.hpp"
#include "BelosTsqrOrthoManager.hpp"
#endif

#ifdef PHIST_HAVE_ANASAZI
#include "phist_trilinos_macros.h"
#include "Anasazi_GhostAdapter.hpp"
#include "phist_AnasaziOperatorTraits.hpp"
#include "AnasaziSVQBOrthoManager.hpp"
#endif

#include "phist_GhostMV.hpp"

#ifdef PHIST_HAVE_LIKWID
#include <likwid.h>
#endif

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <ghost.h>
#include <ghost/machine.h>
#include <ghost/thpool.h>
#include <ghost/pumap.h>
#include <ghost/locality.h>
#include <ghost/timing.h>
#include <limits>
#include <map>

#if defined(PHIST_HAVE_BELOS)||defined(PHIST_HAVE_KOKKOS)
# if defined(GHOST_HAVE_LONGIDX_LOCAL)
# warning "The interfaces between GHOST and Belos/TSQR cause problems unless you compile GHOST with LONGIDX_GLOBAL but *without* LONGIDX_LOCAL"
# endif
#endif

#include "phist_gen_d.h"
#include "kernels_def.hpp"
#include "carp_def.hpp"
#include "../common/kernels_no_fused.cpp"

