#include "phist_config.h"
#include "phist_operator.h"
#include "phist_enums.h"
#include "phist_kernels.h"

#include <stdio.h>

#include "phist_config.h"
#include "phist_ScalarTraits.hpp"

#include "cpBase.hpp"


#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "cpPhistSdMat_def.hpp"

#include "phist_gen_c.h"
#include "cpPhistSdMat_def.hpp"
#endif

#include "phist_gen_d.h"
#include "cpPhistSdMat_def.hpp"

#include "phist_gen_z.h"
#include "cpPhistSdMat_def.hpp"
