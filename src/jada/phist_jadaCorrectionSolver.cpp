#include "phist_config.h"

#include "phist_tools.h"
#include "phist_kernels.h"
#include "phist_core.h"
#include "phist_ScalarTraits.hpp"
#include "phist_jadaCorrectionSolver.h"
#include "phist_jadaOp.h"
#include "phist_precon.h"

#include "phist_MemOwner.hpp"

#include <cstdlib>
#include <vector>

#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_jadaCorrectionSolver_def.hpp"
#include "phist_gen_c.h"
#include "phist_jadaCorrectionSolver_def.hpp"
#endif
#include "phist_gen_d.h"
#include "phist_jadaCorrectionSolver_def.hpp"
#include "phist_gen_z.h"
#include "phist_jadaCorrectionSolver_def.hpp"
