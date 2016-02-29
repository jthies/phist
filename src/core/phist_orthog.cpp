#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "phist_orthog.h"
#include "phist_kernels.h"

// fallback routine: if the kernel library does not implement mvec_QR,
// we simply switch to our own implementation of Cholesky-QR.
#include "phist_chol_QR.h"

#include "phist_ScalarTraits.hpp"

#ifdef PHIST_HAVE_SP
# include "phist_gen_s.h"
# include "phist_orthog_def.hpp"
# include "phist_gen_c.h"
# include "phist_orthog_def.hpp"
#endif

#include "phist_gen_d.h"
# include "phist_orthog_def.hpp"
#include "phist_gen_z.h"
#include "phist_orthog_def.hpp"
