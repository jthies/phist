#ifndef PHIST_JADACORRECTIONSOLVER_H
#define PHIST_JADACORRECTIONSOLVER_H

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "phist_kernels.h"
#include "phist_pgmres.h"
#include "phist_pminres.h"

#ifdef __cplusplus
extern "C" {
#endif
#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_jadaCorrectionSolver_decl.h"
#include "phist_gen_c.h"
#include "phist_jadaCorrectionSolver_decl.h"
#endif
#include "phist_gen_d.h"
#include "phist_jadaCorrectionSolver_decl.h"
#include "phist_gen_z.h"
#include "phist_jadaCorrectionSolver_decl.h"
#ifdef __cplusplus
}
#endif
#endif
