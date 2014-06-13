#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "phist_macros.h"
#include "phist_kernels.h"
#include "phist_enums.h"
#include "phist_feastCorrectionSolver.h"
#include "phist_ScalarTraits.hpp"

#include <cstdlib>
#include <vector>

#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_feastCorrectionSolver_def.hpp"
#endif
#include "phist_gen_d.h"
#include "phist_feastCorrectionSolver_def.hpp"
