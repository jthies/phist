/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
/*! \file phist_typed_test_gen.hpp
 * small helper file to shorten definition of tests for different dimensions and types
*/

// expects the following defined variables:
// class base name _BASENAME_
// dimensions _N_, _M_ and possible _K_
// and possible DISABLE_TESTCASE to prepend DISABLED_

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
#ifdef _K_
#define _NVP_ _K_
#endif

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

// disable single precision kernels for larger test cases because
// they often fail the tolerances suitable for double, e.g.
// 10*mt::eps()
#if defined(PHIST_HAVE_SP)&&((_N_)<=250)
#ifdef DISABLE_TESTCASE
#define CLASSNAME CLASSNAME_FROM_TYPE(DISABLED_S)
#else
#define CLASSNAME CLASSNAME_FROM_TYPE(S)
#endif
#include "phist_gen_s.h"
#include CLASSFILE_DEF
#undef CLASSNAME

#ifdef PHIST_HAVE_CMPLX
#ifdef DISABLE_TESTCASE
#define CLASSNAME CLASSNAME_FROM_TYPE(DISABLED_C)
#else
#define CLASSNAME CLASSNAME_FROM_TYPE(C)
#endif
#include "phist_gen_c.h"
#include CLASSFILE_DEF
#undef CLASSNAME
#endif
#endif

#ifdef DISABLE_TESTCASE
#define CLASSNAME CLASSNAME_FROM_TYPE(DISABLED_D)
#else
#define CLASSNAME CLASSNAME_FROM_TYPE(D)
#endif
#include "phist_gen_d.h"
#include CLASSFILE_DEF
#undef CLASSNAME

#ifdef PHIST_HAVE_CMPLX
#ifdef DISABLE_TESTCASE
#define CLASSNAME CLASSNAME_FROM_TYPE(DISABLED_Z)
#else
#define CLASSNAME CLASSNAME_FROM_TYPE(Z)
#endif
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
#undef _NVP_
#endif
