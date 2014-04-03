#ifndef PHIST_OPERATOR_H
#define PHIST_OPERATOR_H

#include "phist_config.h"
#include "phist_kernels.h"

#ifdef __cplusplus
extern "C" {
#endif
#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_operator_decl.h"
#include "phist_gen_c.h"
#include "phist_operator_decl.h"
#endif
#include "phist_gen_d.h"
#include "phist_operator_decl.h"
#include "phist_gen_z.h"
#include "phist_operator_decl.h"

#ifdef __cplusplus
}
#endif

#endif
