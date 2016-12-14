#include "phist_config.h"

#include "phist_tools.h"
#include "phist_kernels.h"
#include "phist_core.h"
#include "phist_jdqr.h"

#include "phist_schur_decomp.h"
#include "phist_simple_arnoldi.h"
#include "phist_jadaOp.h"
#include "phist_jadaOpts.h"
#include "phist_core_flags.h"

#include "phist_jadaCorrectionSolver.h"

#include <iostream>
#include <iomanip>

#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_jdqr_def.hpp"
#include "phist_gen_c.h"
#include "phist_jdqr_def.hpp"
#endif
#include "phist_gen_d.h"
#include "phist_jdqr_def.hpp"
#include "phist_gen_z.h"
#include "phist_jdqr_def.hpp"
