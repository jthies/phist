#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#ifdef PHIST_HAVE_TEUCHOS
#include <Teuchos_RCP.hpp>
#endif

#include "phist_macros.h"
#include "phist_blockedgmres.h"
#include "phist_kernels.h"
#include "phist_lapack.h"
#include "phist_orthog.h"

#include "phist_ScalarTraits.hpp"

#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_blockedgmres_def.hpp"
#include "phist_gen_c.h"
#include "phist_blockedgmres_def.hpp"
#endif
#include "phist_gen_d.h"
#include "phist_blockedgmres_def.hpp"
#include "phist_gen_z.h"
#include "phist_blockedgmres_def.hpp"

