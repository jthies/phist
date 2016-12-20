#include "phist_config.h"

#include "phist_tools.h"
#include "phist_kernels.h"
#include "phist_iter_op.h"

// todo: rename feastCorrectionSolver,
// it is in fact a kind of general wrapper
// for our linear solvers.
#include "phist_feastCorrectionSolver.h"

#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_iter_op_def.hpp"
#include "phist_gen_c.h"
#include "phist_iter_op_def.hpp"
#endif

#include "phist_gen_d.h"
#include "phist_iter_op_def.hpp"
#include "phist_gen_z.h"
#include "phist_iter_op_def.hpp"

