#include "phist_feastCorrectionSolver.h"
#include "phist_macros.h"
#include "phist_enums.h"
#include "phist_ScalarTraits.hpp"

#include <cstdlib>
#include <vector>

#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_jadaCorrectionSolver_def.hpp"
#endif
#include "phist_gen_d.h"
#include "phist_jadaCorrectionSolver_def.hpp"
