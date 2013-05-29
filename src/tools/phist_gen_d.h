#include "phist_gen_pre.h"

// type specifier
#define _TP_ D

#define _IS_DOUBLE_

// scalar type
#define _ST_ double

#define _ZERO_ 0.0

#define _ONE_ 1.0

// adds type prefix to a specifier, e.g. _PREF_(gemm) -> Dgemm
#define _PREF_(name) D ## name

// how to build up the name of a subroutine (void function)
#define _SUBR_(name) phist_D ## name

// how to build up the name of a type
#define _TYPE_(name) D ## name ## _t

// how to build up the name of a test with two 'template' params
#define _TESTNAME2_(name,p1,p2) D ## name ## _ ## p1 ## _ ## p2

#include "phist_gen_post.h"
