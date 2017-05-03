#ifndef PHIST_F_MACROS_H
#define PHIST_F_MACROS_H

#include "phist_config.h"
#ifdef PHIST_HAVE_GHOST
#include <ghost/config.h>
#endif

#ifndef PHIST_HAVE_GHOST
#define G_LIDX_T C_INT32_T
# ifdef PHIST_FORCE_INT_GIDX
# define G_GIDX_T C_INT32_T
# else
# define G_GIDX_T C_INT64_T
# endif
#else
# ifdef GHOST_IDX64_LOCAL
# define G_LIDX_T C_INT64_T
# else
# define G_LIDX_T C_INT32_T
# endif
# ifdef GHOST_IDX64_GLOBAL
# define G_GIDX_T C_INT64_T
# else
# define G_GIDX_T C_INT32_T
# endif
#endif

#include "fdebug.h"

#endif
