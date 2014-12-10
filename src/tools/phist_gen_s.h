#include "phist_gen_clean.h"

// scalar type
#define _ST_ float

#define PREFIX(name) S ## name
#define SPREFIX(name) s ## name
#define SUBR(name) phist_S ## name
#define TYPE(name) S ## name ## _t

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
