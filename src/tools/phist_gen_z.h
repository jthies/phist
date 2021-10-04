/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_gen_clean.h"

// type specifier
#define _TP_ 'Z'
#define _ST_ phist_d_complex

#define IS_DOUBLE
#define IS_COMPLEX

// C++ users should use the class phist::ScalarTraits and
// phist_std_typedefs.hpp for all of this:
#ifndef __cplusplus
#define _CMPLX_I_ (0.0+1.0*I)
#define ZERO (0.0+0.0*I)
#define ONE (1.0+0.0*I)
#define SQRT(X) csqrt(X)
#define MSQRT(X) sqrt(X)
#define ABS(X) cabs(X)
#define MABS(X) abs(X)
#define CONJ(X) conj(X)
#define REAL(X) creal(X)
#define IMAG(X) cimag(X)
#endif

// adds type prefix to a specifier, e.g. PHIST_TG_PREFIX(gemm) -> Dgemm
#define PHIST_TG_PREFIX(name) Z ## name

// adds lower case type prefix to a specifier, e.g. SPHIST_TG_PREFIX(gemm) -> dgemm
#define SPHIST_TG_PREFIX(name) z ## name

// how to build up the name of a subroutine (void function)
#define SUBR(name) phist_Z ## name

// how to build up the name of a type
#define TYPE(name) phist_Z ## name

// how to call a lapack routine via the C interface, e.g. PHIST_LAPACKE(getrf)
#define PHIST_LAPACKE(name) LAPACKE_z ## name

#include "phist_gen_common.h"
