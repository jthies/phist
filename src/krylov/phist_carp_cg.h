#ifndef PHIST_PCG_H
#define PHIST_PCG_H

#include "phist_config.h"

#ifndef DOXYGEN

/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "phist_enums.h"
#include "phist_typedefs.h"

#endif //DOXYGEN

#ifdef __cplusplus
extern "C" {
#endif
#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_carp_cg_decl.h"
#include "phist_gen_c.h"
#include "phist_carp_cg_decl.h"
#endif
#include "phist_gen_d.h"
#include "phist_carp_cg_decl.h"
#include "phist_gen_z.h"
#include "phist_carp_cg_decl.h"
#ifdef __cplusplus
}
#endif

#endif
