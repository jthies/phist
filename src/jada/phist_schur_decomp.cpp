#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <cstdlib>
#include <cstring>
#include <cmath>

#include "phist_macros.h"
#include "phist_schur_decomp.h"
#include "jada_helpers.hpp"
#include "phist_kernels.h"
#include "phist_lapack.h"

#include "phist_ScalarTraits.hpp"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_schur_decomp_def.hpp"
#include "phist_gen_c.h"
#include "phist_schur_decomp_def.hpp"
#endif
#include "phist_gen_d.h"
#include "phist_schur_decomp_def.hpp"
#include "phist_gen_z.h"
#include "phist_schur_decomp_def.hpp"

#ifdef __cplusplus
}
#endif
