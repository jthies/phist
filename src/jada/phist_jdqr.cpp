#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <iomanip>

#include "phist_macros.h"
#include "phist_jdqr.h"
#include "phist_kernels.h"
#include "phist_lapack.h"
#include "phist_orthog.h"

#include "phist_ScalarTraits.hpp"
#include "phist_schur_decomp.h"
#include "phist_simple_arnoldi.h"
#include "phist_jadaOp.hpp"
#include "phist_jadaOpts.h"

#include "phist_jadaCorrectionSolver.h"
// for testing CARP-CG, we use Belos BlockCG for now.
// If it turns out to work well, we could include CG
// as an option in the jadaCorrectionSolver.
#include "phist_belos.h"

#include "phist_gen_s.h"
#include "phist_jdqr_def.hpp"
#include "phist_gen_d.h"
#include "phist_jdqr_def.hpp"
#include "phist_gen_c.h"
#include "phist_jdqr_def.hpp"
#include "phist_gen_z.h"
#include "phist_jdqr_def.hpp"

