#ifndef PHIST_JADACORRECTIONSOLVER_H
#define PHIST_JADACORRECTIONSOLVER_H

#include "phist_kernels.h"

#ifdef __cplusplus
extern "C" {
#endif
// TODO - for now we only provide the real version of
// this solver interface, but do allow complex shifts
#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_feastCorrectionSolver_decl.h"
#endif
#include "phist_gen_d.h"
#include "phist_feastCorrectionSolver_decl.h"
#ifdef __cplusplus
}
#endif
#endif
