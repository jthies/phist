#ifndef PHIST_BGMRES_H
#define PHIST_BGMRES_H

#include "phist_operator.h"

#ifdef __cplusplus
extern "C" {
#endif

#include "phist_gen_d.h"
#include "phist_bgmres_decl.h"

// disallow using Belos with other than double type for Epetra,
// that makes no sense because Epetra has only double as type.
#ifndef PHIST_KERNEL_LIB_EPETRA
#include "phist_gen_s.h"
#include "phist_bgmres_decl.h"
#include "phist_gen_c.h"
#include "phist_bgmres_decl.h"
#include "phist_gen_z.h"
#include "phist_bgmres_decl.h"
#endif

#ifdef __cplusplus
}
#endif

#endif
