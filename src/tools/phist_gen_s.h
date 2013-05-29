#include "phist_gen_pre.h"

// type specifier
#define _TP_ S

// scalar type
#define _ST_ float

#define _ZERO_ 0.0f

#define _ONE_ 1.0f

// adds type prefix to a specifier, e.g. _PREF_(gemm) -> Dgemm
#define _PREF_(name) S ## name

// how to build up the name of a subroutine (void function)
#define _SUBR_(name) phist_S ## name

// how to build up the name of a type
#define _TYPE_(name) S ## name ## _t

// how to build up the name of a test with two 'template' params
#define _TESTNAME2_(name,p1,p2) S ## name ## _ ## p1 ## _ ## p2

#include "phist_gen_post.h"
