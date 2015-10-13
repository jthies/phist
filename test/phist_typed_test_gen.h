/*! \file phist_typed_test_gen.hpp
 * small helper file to shorten definition of tests for different dimensions and types
*/

// expects the following defined variables:
// class base name _BASENAME_
// dimensions _N_, _M_ and possible _K_

// construct CLASSNAME from basename, type and dimensions
#ifdef _K_

#define CN_EVAL_CONCATENATE(a,b,c,d,e) CN_CONCATENATE(a,b,c,d,e)
#define CN_CONCATENATE(a,b,c,d,e) a ## b ## _ ## c ## _ ## d ## _ ## e
#define CLASSNAME_FROM_TYPE(t) CN_EVAL_CONCATENATE(t,_BASENAME_,_N_,_M_,_K_)

#else

#define CN_EVAL_CONCATENATE(a,b,c,d) CN_CONCATENATE(a,b,c,d)
#define CN_CONCATENATE(a,b,c,d) a ## b ## _ ## c ## _ ## d
#define CLASSNAME_FROM_TYPE(t) CN_EVAL_CONCATENATE(t,_BASENAME_,_N_,_M_)

#endif

// some tests use other dimension specifiers
#define _NV_ _M_
#define _NROWS_ _N_
#define _NCOLS_ _M_

// construct test def file name. The user may overrule this by pre-defining CLASSFILE_DEF
#ifndef CLASSFILE_DEF
#define CLASSFILE_DEF_WAS_UNDEFINED
#define CLASSFILE_DEF CF_EVAL_CONCATENATE_AND_QUOTE(_BASENAME_)
#define CF_EVAL_QUOTE(y) CF_QUOTE(y)
#define CF_QUOTE(x) #x
#define CF_EVAL_CONCATENATE_AND_QUOTE(x) CF_CONCATENATE_AND_QUOTE(x)
#define CF_CONCATENATE_AND_QUOTE(x) CF_EVAL_QUOTE(x##_def.hpp)
#endif

// different types
#ifdef PHIST_HAVE_SP
#define CLASSNAME CLASSNAME_FROM_TYPE(S)
#include "phist_gen_s.h"
#include CLASSFILE_DEF
#undef CLASSNAME

#ifdef PHIST_HAVE_CMPLX
#define CLASSNAME CLASSNAME_FROM_TYPE(C)
#include "phist_gen_c.h"
#include CLASSFILE_DEF
#undef CLASSNAME
#endif
#endif

#define CLASSNAME CLASSNAME_FROM_TYPE(D)
#include "phist_gen_d.h"
#include CLASSFILE_DEF
#undef CLASSNAME

#ifdef PHIST_HAVE_CMPLX
#define CLASSNAME CLASSNAME_FROM_TYPE(Z)
#include "phist_gen_z.h"
#include CLASSFILE_DEF
#undef CLASSNAME
#endif

// remove definitions
#undef CF_EVAL_QUOTE
#undef CF_QUOTE
#undef CF_EVAL_CONCATENATE_AND_QUOTE
#undef CF_CONCATENATE_AND_QUOTE

#undef CLASSNAME_FROM_TYPE
#undef CN_CONCATENATE
#undef CN_EVAL_CONCATENATE

#ifdef CLASSFILE_DEF_WAS_UNDEFINED
#undef CLASSFILE_DEF
#undef CLASSFILE_DEF_WAS_UNDEFINED
#endif
#undef BASENAME_DIM

#undef _N_
#undef _M_
#undef _NV_
#undef _NROWS_
#undef _NCOLS_
#ifdef _K_
#undef _K_
#endif
