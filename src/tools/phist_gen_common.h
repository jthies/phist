/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
// some common macros for all types:
#ifdef IS_DOUBLE
#define _MT_ double
#define _CT_ phist_d_complex
#define ASSERT_REAL_EQ(expected,actual) ASSERT_DOUBLE_EQ(expected,actual)
#define EXPECT_REAL_EQ(expected,actual) EXPECT_DOUBLE_EQ(expected,actual)
#else
#define _MT_ float
#define _CT_ phist_s_complex
#define ASSERT_REAL_EQ(expected,actual) ASSERT_FLOAT_EQ(expected,actual)
#define EXPECT_REAL_EQ(expected,actual) EXPECT_FLOAT_EQ(expected,actual)
#endif

// check if the standard things are defined in a previously included gen_x.h file
#ifndef _ST_
#error "macro _ST_ not defined"
#endif

#ifndef PHIST_TG_PREFIX
#error "macro PHIST_TG_PREFIX not defined"
#endif

#ifndef SPHIST_TG_PREFIX
#error "macro SPHIST_TG_PREFIX not defined"
#endif

#ifndef PHIST_LAPACKE
#error "macro PHIST_LAPACKE not defined"
#endif

#ifndef SUBR
#error "macro SUBR not defined"
#endif

#ifndef TYPE
#error "macro TYPE not defined"
#endif

// additional things for the C interfaces
#ifndef __cplusplus

#ifndef IS_COMPLEX
#define CONJ(x) x
#define REAL(x) x
#define IMAG(x) ZERO
#define _CMPLX_I_ ZERO
#define MSQRT(x) SQRT(x)
#define MABS(x) ABS(x)
#endif

// check wether all necessary macros are there:
#ifndef _TP_
#error "macro _TP_ not defined"
#endif

#ifndef _CMPLX_I_
#warning "macro _CMPLX_I_ not defined"
#endif

#ifndef CONJ
#error "macro CONJ not defined"
#endif

#ifndef SQRT
#error "macro SQRT not defined"
#endif

#ifndef ABS
#error "macro ABS not defined"
#endif

#ifndef REAL
#error "macro REAL not defined"
#endif

#ifndef IMAG
#error "macro IMAG not defined"
#endif

#ifndef ZERO
#error "macro ZERO not defined"
#endif

#ifndef ONE
#error "macro ONE not defined"
#endif
#endif
