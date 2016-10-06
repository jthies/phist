#include "phist_config.h"
#ifdef PHIST_HAVE_SP
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "phist_macros.h"
#include "phist_subspacejada.h"
#include "phist_kernels.h"
#include "phist_orthog.h"

#include "phist_ScalarTraits.hpp"
#include "phist_schur_decomp.h"
#include "phist_simple_arnoldi.h"
#include "phist_transform_searchspace.h"
#include "phist_jadaCorrectionSolver.h"
#include "phist_core_flags.h"


#include "phist_gen_c.h"
#include "phist_subspacejada_def.hpp"

#endif
