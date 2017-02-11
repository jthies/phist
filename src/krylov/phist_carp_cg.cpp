#include "phist_config.h"

#include "phist_tools.h"
#include "phist_kernels.h"
#include "phist_core.h"
#include "phist_ScalarTraits.hpp"

#include <vector>

#include "phist_carp_cg.h"

#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_carp_cg_def.hpp"
#include "phist_gen_c.h"
#include "phist_carp_cg_def.hpp"
#endif
#include "phist_gen_d.h"
#include "phist_carp_cg_def.hpp"
#include "phist_gen_z.h"
#include "phist_carp_cg_def.hpp"
