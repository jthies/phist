/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/

#ifdef COMPLEX_T
#undef COMPLEX_T
#endif

// this file undefs all the macros defined by the phist_gen_* 
// headers
#ifdef _TP_
#undef _TP_
#endif

#ifdef _ST_
#undef _ST_
#endif

#ifdef _CMPLX_I_
#undef _CMPLX_I_
#endif

#ifdef _MT_
#undef _MT_
#endif

#ifdef _CT_
#undef _CT_
#endif

#ifdef CONJ
#undef CONJ
#endif

#ifdef SQRT
#undef SQRT
#endif

#ifdef MSQRT
#undef MSQRT
#endif

#ifdef ABS
#undef ABS
#endif

#ifdef MABS
#undef MABS
#endif

#ifdef REAL
#undef REAL
#endif

#ifdef IMAG
#undef IMAG
#endif

#ifdef IS_COMPLEX
#undef IS_COMPLEX
#endif

#ifdef IS_DOUBLE
#undef IS_DOUBLE
#endif

#ifdef ZERO
#undef ZERO
#endif

#ifdef ONE
#undef ONE
#endif

#ifdef PHIST_TG_PREFIX
#undef PHIST_TG_PREFIX
#endif

#ifdef SPHIST_TG_PREFIX
#undef SPHIST_TG_PREFIX
#endif

#ifdef SUBR
#undef SUBR
#endif

#ifdef TYPE
#undef TYPE
#endif

#ifdef PHIST_LAPACKE
#undef PHIST_LAPACKE
#endif

#ifdef ASSERT_REAL_EQ
#undef ASSERT_REAL_EQ
#endif

#ifdef EXPECT_REAL_EQ
#undef EXPECT_REAL_EQ
#endif
