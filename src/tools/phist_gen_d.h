#include "phist_gen_clean.h"

#define IS_DOUBLE

// scalar type
#define _ST_ double

// C++ users should use the class phist::ScalarTraits and
// phist_std_typedefs.hpp for all of this:
#ifndef __cplusplus

// type specifier
#define _TP_ 'D'

#define ZERO 0.0

#define ONE 1.0

// squareroot
#define SQRT(X) sqrt(X)

// absolute value
#define ABS(X) abs(X)

#endif

// adds type prefix to a specifier, e.g. PHIST_TG_PREFIX(gemm) -> Dgemm
#define PHIST_TG_PREFIX(name) D ## name

// adds lower case type prefix to a specifier, e.g. SPHIST_TG_PREFIX(gemm) -> dgemm
#define SPHIST_TG_PREFIX(name) d ## name

// how to build up the name of a subroutine (void function)
#define SUBR(name) phist_D ## name

// how to build up the name of a type
#define TYPE(name) phist_D ## name

// how to call a lapack routine via the C interface, e.g. PHIST_LAPACKE(getrf)
#define PHIST_LAPACKE(name) LAPACKE_d ## name

#include "phist_gen_common.h"
