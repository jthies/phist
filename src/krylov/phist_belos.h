#ifndef PHIST_BELOS_H
#define PHIST_BELOS_H

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "phist_operator.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_belos_decl.h"
#include "phist_gen_c.h"
#include "phist_belos_decl.h"
#endif
#include "phist_gen_d.h"
#include "phist_belos_decl.h"
#include "phist_gen_z.h"
#include "phist_belos_decl.h"

#ifdef __cplusplus
}
#endif

#endif
