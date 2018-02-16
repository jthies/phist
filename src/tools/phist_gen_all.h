#ifndef PHIST_CLASSFILE_DEF
#error "file included incorrectly"
#endif

#include "phist_config.h"

#include "phist_gen_d.h"
#include PHIST_CLASSFILE_DEF

#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include PHIST_CLASSFILE_DEF
#endif

#ifdef PHIST_HAVE_CMPLX

#include "phist_gen_z.h"
#include PHIST_CLASSFILE_DEF

#ifdef PHIST_HAVE_SP
#include "phist_gen_c.h"
#include PHIST_CLASSFILE_DEF
#endif

#endif

/* cleanup macros */
#include "phist_gen_clean.h"
#undef PHIST_CLASSFILE_DEF
