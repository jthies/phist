#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "phist_jadaOp.hpp"
#include "phist_macros.h"
#include "phist_ScalarTraits.hpp"

#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_jadaOp_def.hpp"
#include "phist_gen_c.h"
#include "phist_jadaOp_def.hpp"
#endif
#include "phist_gen_d.h"
#include "phist_jadaOp_def.hpp"
#include "phist_gen_z.h"
#include "phist_jadaOp_def.hpp"
