#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <stdlib.h>

#include "phist_kernels.h"
#include "phist_random.h"

typedef struct {
  int lda;
  double* data;
} dwrap;


#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "kernels_common_impl_def.h"

#include "phist_gen_c.h"
#include "kernels_common_impl_def.h"
#endif

#include "phist_gen_d.h"
#include "kernels_common_impl_def.h"

#include "phist_gen_z.h"
#include "kernels_common_impl_def.h"

