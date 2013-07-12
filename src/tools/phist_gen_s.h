#include "phist_gen_pre.h"

// scalar type
#define _ST_ float

#define _PREF_(name) S ## name
#define _SUBR_(name) phist_S ## name
#define _TYPE_(name) S ## name ## _t

// C++ users should use the class phist::ScalarTraits and
// phist_std_typedefs.hpp for all of this:
#ifndef __cplusplus

#define _SQRT_(X) sqrtf(X)
#define _ABS_(X) fabs(X)

// type specifier
#define _TP_ 'S'

#define _ZERO_ 0.0f

#define _ONE_ 1.0f
#endif

#include "phist_gen_post.h"
