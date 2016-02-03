#ifndef PHIST_PRECON_H
#define PHIST_PRECON_H

#include "phist_operator.h"

#include "phist_config.h"

#ifndef DOXYGEN

/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "phist_kernels.h"

#endif //DOXYGEN

typedef enum {
  NONE,
#ifdef PHIST_HAVE_IFPACK
  IFPACK,
#endif
#ifdef PHIST_HAVE_ML
  ML,
#endif
#ifdef PHIST_HAVE_IFPACK2
  IFPACK2,
#endif
#ifdef PHIST_HAVE_MUELU
  MUELU,
#endif
  INVALID_PRECON_TYPE
} precon_t;

#ifdef __cplusplus
extern "C" {
#endif
#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_precon_decl.h"
#include "phist_gen_c.h"
#include "phist_precon_decl.h"
#endif
#include "phist_gen_d.h"
#include "phist_precon_decl.h"
#include "phist_gen_z.h"
#include "phist_precon_decl.h"

#ifdef __cplusplus
}
#endif


#endif
