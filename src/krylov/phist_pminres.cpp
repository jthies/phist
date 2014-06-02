#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>

#include "phist_config.h"
#include "phist_macros.h"
#include "phist_pminres.h"
#include "phist_jadaOp.hpp"
#include "phist_kernels.h"
#include "phist_lapack.h"
#include "phist_orthog.h"
#include "phist_pgmres.h"

#include "phist_ScalarTraits.hpp"

#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_pminres_def.hpp"
#include "phist_gen_c.h"
#include "phist_pminres_def.hpp"
#endif
#include "phist_gen_d.h"
#include "phist_pminres_def.hpp"
#include "phist_gen_z.h"
#include "phist_pminres_def.hpp"

