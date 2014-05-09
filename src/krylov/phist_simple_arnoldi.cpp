#include <cstdlib>
#include <cstring>
#include <cmath>

#include "phist_config.h"
#include "phist_macros.h"
#include "phist_simple_arnoldi.h"
#include "phist_kernels.h"
#include "phist_lapack.h"
#include "phist_orthog.h"

#include "phist_ScalarTraits.hpp"
#include "jada_helpers.hpp"

#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_simple_arnoldi_def.hpp"
#include "phist_gen_c.h"
#include "phist_simple_arnoldi_def.hpp"
#endif
#include "phist_gen_d.h"
#include "phist_simple_arnoldi_def.hpp"
#include "phist_gen_z.h"
#include "phist_simple_arnoldi_def.hpp"

