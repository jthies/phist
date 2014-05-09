#include "phist_config.h"
#include "phist_transform_searchspace.h"
#include "phist_macros.h"
#include "phist_lapack.h"
#include "phist_ScalarTraits.hpp"

#include <stdlib.h>

#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_transform_searchspace_def.hpp"
#include "phist_gen_c.h"
#include "phist_transform_searchspace_def.hpp"
#endif
#include "phist_gen_d.h"
#include "phist_transform_searchspace_def.hpp"
#include "phist_gen_z.h"
#include "phist_transform_searchspace_def.hpp"
