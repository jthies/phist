#include <cstdlib>
#include <cstring>
#include <cmath>

#include "phist_macros.h"
#include "phist_schur_decomp.h"
#include "phist_kernels.h"
#include "phist_lapack.h"

#include "phist_ScalarTraits.hpp"
#include "jada_helpers.hpp"

#ifdef __cplusplus
extern "C" {
#endif

#include "phist_gen_s.h"
#include "phist_schur_decomp_def.hpp"
#include "phist_gen_d.h"
#include "phist_schur_decomp_def.hpp"
#include "phist_gen_c.h"
#include "phist_schur_decomp_def.hpp"
#include "phist_gen_z.h"
#include "phist_schur_decomp_def.hpp"

#ifdef __cplusplus
}
#endif