#include "phist_config.h"

#ifdef PHIST_HAVE_SP

/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "phist_operator.h"
#include "phist_precon.h"
#include "phist_macros.h"
#include "phist_ScalarTraits.hpp"
#include <cstdlib>
#include <strings.h>

#include "phist_PreconTraits.hpp"

#include "phist_gen_s.h"
#include "phist_precon_def.hpp"
#endif
