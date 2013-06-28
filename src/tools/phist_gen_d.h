#include "phist_gen_pre.h"

// type specifier
#define _TP_ D

#define _IS_DOUBLE_

// scalar type
#define _ST_ double

#define _ZERO_ 0.0

#define _ONE_ 1.0

// squareroot
#define _SQRT_(X) sqrt(X)

// absolute value
#define _ABS_(X) abs(X)

// adds type prefix to a specifier, e.g. _PREF_(gemm) -> Dgemm
#define _PREF_(name) D ## name

// how to build up the name of a subroutine (void function)
#define _SUBR_(name) phist_D ## name

// how to build up the name of a type
#define _TYPE_(name) D ## name ## _t

#include "phist_gen_post.h"
