#ifndef PHIST_JADA_OP_H
#define PHIST_JADA_OP_H

#include "phist_config.h"

#ifndef DOXYGEN

/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "phist_operator.h"

#endif //DOXYGEN

#ifdef __cplusplus
extern "C" {
#endif

#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_jadaOp_decl.h"
#include "phist_gen_c.h"
#include "phist_jadaOp_decl.h"
#endif
#include "phist_gen_d.h"
#include "phist_jadaOp_decl.h"
#include "phist_gen_z.h"
#include "phist_jadaOp_decl.h"
#include "phist_gen_clean.h"

#ifdef __cplusplus
} // extern "C"
#endif

#endif