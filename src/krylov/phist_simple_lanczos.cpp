#include "phist_config.h"

#include "phist_tools.h"
#include "phist_kernels.h"
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "phist_simple_lanczos.h"

#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_simple_lanczos_def.hpp"
#include "phist_gen_c.h"
#include "phist_simple_lanczos_def.hpp"
#endif
#include "phist_gen_d.h"
#include "phist_simple_lanczos_def.hpp"
#include "phist_gen_z.h"
#include "phist_simple_lanczos_def.hpp"

