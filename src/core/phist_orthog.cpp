#include "phist_config.h"
#include "phist_tools.h"
#include "phist_kernels.h"
#include "phist_orthog.h"

// fallback routine: if the kernel library does not implement mvec_QR,
// we simply switch to our own implementation of Cholesky-QR.
#include "phist_chol_QR.h"
// under the hood we use this implementation, which exploits fused kernels and
// high precision operations if available.
#include "phist_orthogrrfused.h"
#include "phist_orthogrr.h"

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
