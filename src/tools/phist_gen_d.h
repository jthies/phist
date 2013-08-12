#include "phist_gen_pre.h"

#define _IS_DOUBLE_

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

// adds type prefix to a specifier, e.g. PREFIX(gemm) -> Dgemm
#define PREFIX(name) D ## name

// how to build up the name of a subroutine (void function)
#define SUBR(name) phist_D ## name

// how to build up the name of a type
#define TYPE(name) D ## name ## _t

#include "phist_gen_post.h"
