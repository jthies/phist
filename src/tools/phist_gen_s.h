#include "phist_gen_pre.h"

// scalar type
#define _ST_ float

#define PREFIX(name) S ## name
#define SUBR(name) phist_S ## name
#define TYPE(name) S ## name ## _t

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

#include "phist_gen_post.h"
