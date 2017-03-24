/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_gen_clean.h"

// scalar type
#define _ST_ float

#define PHIST_TG_PREFIX(name) S ## name
#define SPHIST_TG_PREFIX(name) s ## name
#define SUBR(name) phist_S ## name
#define TYPE(name) phist_S ## name

// how to call a lapack routine via the C interface, e.g. PHIST_LAPACKE(getrf)
#define PHIST_LAPACKE(name) LAPACKE_s ## name

// C++ users should use the class phist::ScalarTraits and
// phist_std_typedefs.hpp for all of this:
#ifndef __cplusplus

#define SQRT(X) sqrtf(X)
#define ABS(X) fabs(X)

// type specifier
#define _TP_ 'S'

#define ZERO 0.0f

#define ONE 1.0f
#endif

#include "phist_gen_common.h"
