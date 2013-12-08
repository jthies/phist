#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <iomanip>

#include "phist_macros.h"
#include "phist_subspacejada.h"
#include "phist_kernels.h"
#include "phist_lapack.h"
#include "phist_orthog.h"

#include "phist_ScalarTraits.hpp"
#include "jada_helpers.hpp"
#include "phist_schur_decomp.h"
#include "phist_jadaOp.h"
#include "phist_simple_arnoldi.h"
#include "phist_mvec_times_sdMat_inplace.h"

#include "phist_bgmres.h"

#include "phist_gen_s.h"
#include "phist_subspacejada_def.hpp"
#include "phist_gen_d.h"
#include "phist_subspacejada_def.hpp"
#include "phist_gen_c.h"
#include "phist_subspacejada_def.hpp"
#include "phist_gen_z.h"
#include "phist_subspacejada_def.hpp"

