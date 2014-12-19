#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "phist_svqb.h"

#include "phist_macros.h"
#include "phist_kernels.h"
#include "phist_lapack.h"

#include "phist_ScalarTraits.hpp"

#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_svqb_def.hpp"
#include "phist_gen_c.h"
#include "phist_svqb_def.hpp"
#endif
#include "phist_gen_d.h"
#include "phist_svqb_def.hpp"
#include "phist_gen_z.h"
#include "phist_svqb_def.hpp"