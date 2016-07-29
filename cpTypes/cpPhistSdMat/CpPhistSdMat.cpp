#include "phist_config.h"
#include "phist_operator.h"
#include "phist_enums.h"
#include "phist_kernels.h"

#include <stdio.h>

#include "phist_config.h"
#include "phist_ScalarTraits.hpp"

#include "../../CpBase.hpp"


#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "CpPhistSdMat_def.hpp"

#include "phist_gen_c.h"
#include "CpPhistSdMat_def.hpp"
#endif

#include "phist_gen_d.h"
#include "CpPhistSdMat_def.hpp"

#include "phist_gen_z.h"
#include "CpPhistSdMat_def.hpp"
