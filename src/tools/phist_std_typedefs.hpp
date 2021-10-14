/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/

//! \file phist_std_typedefs.hpp
//! you can include this file in a scope to get C++ typedefs for
//! the  macros defined in a previous include of phist_gen_x.h. 
//! this allows easier access to traits classes and saves the C++
//! programmer from writing down all the macros like TYPE(...)
//! all the time. The programmer must take care that this file
//! is included once per scope only, otherwise the compiler will 
//! complain.

#ifndef DOXYGEN

#ifndef __cplusplus
#error this file is only intended for C++ code. C/Fortran users should \
stick to the macros defined in 'phist_gen_x.h'
#endif

#ifndef PHIST_SCALAR_TRAITS_HPP
#error you have to include 'phist_ScalarTraits.hpp' before \
'phist_std_typedefs.hpp'.
#endif

typedef ::phist::ScalarTraits< _ST_ >::scalar_t ST;
typedef ::phist::ScalarTraits<ST> st;

typedef st::magn_t MT;
typedef ::phist::ScalarTraits<MT> mt;

typedef std::complex<MT> CT;
typedef ::phist::ScalarTraits<CT> ct;

typedef st::mvec_t* mvec_ptr;
typedef const st::mvec_t* const_mvec_ptr;
typedef st::sparseMat_t* sparseMat_ptr;
typedef const st::sparseMat_t* const_sparseMat_ptr;
typedef st::sdMat_t* sdMat_ptr;
typedef const st::sdMat_t* const_sdMat_ptr;
typedef st::linearOp_t* linearOp_ptr;
typedef const st::linearOp_t* const_linearOp_ptr;

typedef mt::blas_cmplx_t blas_cmplx;

#endif /* DOXYGEN */
