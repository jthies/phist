#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "phist_macros.h"
#include "phist_simple_lanczos_ft.h"
#include "phist_kernels.h"
#include "phist_lapack.h"

#include "phist_ScalarTraits.hpp"

#ifndef GHOST_CP
#define GHOST_CP
#endif
#ifndef PHIST_CP
#define PHIST_CP
#endif

#include "Checkpoint.hpp"
#include "cpTypes/cpPhistMvec/CpPhistMvec.h"
#include "cp_options.h"



#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_simple_lanczos_def_ft.hpp"
#include "phist_gen_c.h"
#include "phist_simple_lanczos_def_ft.hpp"
#endif

#include "phist_gen_d.h"
#include "phist_simple_lanczos_def_ft.hpp"
#include "phist_gen_z.h"
#include "phist_simple_lanczos_def_ft.hpp"
