/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_gen_clean.h"

#define IS_COMPLEX
#define _ST_ phist_s_complex

// C++ users should use the class phist::ScalarTraits and
// phist_std_typedefs.hpp for all of this:
#ifndef __cplusplus

// type specifier
#define _TP_ 'C'

#define _CMPLX_I_ (0.0f+1.0f*I)
#define ZERO (0.0f+0.0f*I)
#define ONE (1.0f+0.0f*I)
#define SQRT(X) csqrtf(X)
#define MSQRT(X) sqrtf(X)
#define CONJ(X) conjf(X)
#define ABS(X) cabsf(X)
#define MABS(X) absf(X)
#define REAL(X) crealf(X)
#define IMAG(X) cimagf(X)

#endif

// adds type prefix to a specifier, e.g. PHIST_TG_PREFIX(gemm) -> Dgemm
#define PHIST_TG_PREFIX(name) C ## name

// adds lower case type prefix to a specifier, e.g. SPHIST_TG_PREFIX(gemm) -> dgemm
#define SPHIST_TG_PREFIX(name) c ## name

// how to build up the name of a subroutine (void function)
#define SUBR(name) phist_C ## name

// how to build up the name of a type
#define TYPE(name) phist_C ## name

// how to call a lapack routine via the C interface, e.g. PHIST_LAPACKE(getrf)
#define PHIST_LAPACKE(name) LAPACKE_c ## name

#include "phist_gen_common.h"
