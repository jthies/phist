#ifdef __cplusplus
#include <cstring>
#else
#include <string.h>
#endif

#include "phist_config.h"
#include "phist_kernels.h"

int phist_sizeof_lidx()
{
  return sizeof(lidx_t);
}

int phist_sizeof_gidx()
{
  return sizeof(gidx_t);
}

#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_driver_utils_def.h"
#include "phist_gen_c.h"
#include "phist_driver_utils_def.h"
#endif
#include "phist_gen_d.h"
#include "phist_driver_utils_def.h"
#include "phist_gen_z.h"
#include "phist_driver_utils_def.h"

