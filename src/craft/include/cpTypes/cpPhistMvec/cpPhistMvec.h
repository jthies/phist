#ifndef __CPPHISTMVEC_H__
#define __CPPHISTMVEC_H__


#include "phist_config.h"
#include "phist_operator.h"
#include "phist_enums.h"
#include "phist_kernels.h"

#ifdef __cplusplus
extern "C" {
#endif
#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "cpPhistMvec_decl.h"
#include "phist_gen_c.h"
#include "cpPhistMvec_decl.h"
#endif

#include "phist_gen_d.h"
#include "cpPhistMvec_decl.h"

#include "phist_gen_z.h"
#include "cpPhistMvec_decl.h"

#include "phist_gen_clean.h"

#ifdef __cplusplus
}
#endif


#endif
