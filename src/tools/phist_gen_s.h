#include "phist_gen_pre.h"

// type specifier
#define _TP_ S

// scalar type
#define _ST_ float

#define _ZERO_ 0.0f

#define _ONE_ 1.0f

#ifdef __cplusplus
#define _SQRT_(X) std::sqrt(X)
#define _ABS_(X) std::abs(X)
#else
#define _SQRT_(X) sqrtf(X)
#define _ABS_(X) fabs(X)
#endif
#define _PREF_(name) S ## name
#define _SUBR_(name) phist_S ## name
#define _TYPE_(name) S ## name ## _t
#include "phist_gen_post.h"
