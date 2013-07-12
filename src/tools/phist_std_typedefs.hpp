#ifndef __cplusplus
#error this file is only intended for C++ code. C/Fortran users should \
stick to the macros defined in 'phist_gen_x.h'
#endif

#ifndef PHIST_SCALAR_TRAITS_HPP
#error you have to include 'phist_ScalarTraits.hpp' before \
'phist_std_typedefs.hpp'.
#endif

//! you can include this file in a scope to get C++ typedefs for
//! the  macros defined in a previous include of phist_gen_x.h. 
//! this allows easier access to traits classes and saves the C++
//! programmer from writing down all the macros like _TYPE_(...)
//! all the time. The programmer must take care that this file
//! is included once per scope only, otherwise the compiler will 
//! complain.
typedef ::phist::ScalarTraits< _ST_ > st;
typedef st::scalar_t ST;
typedef st::magn_t MT;
typedef ::phist::ScalarTraits<MT> mt;
typedef _TYPE_(mvec_ptr) mvec_ptr_t;
typedef _TYPE_(const_mvec_ptr) const_mvec_ptr_t;
typedef _TYPE_(sdMat_ptr) sdMat_ptr_t;
typedef _TYPE_(const_sdMat_ptr) const_sdMat_ptr_t;
typedef _TYPE_(crsMat_ptr) crsMat_ptr_t;
typedef _TYPE_(const_crsMat_ptr) const_crsMat_ptr_t;
#ifdef PHIST_OPERATOR_H
typedef _TYPE_(op_ptr) op_ptr_t;
typedef _TYPE_(const_op_ptr) const_op_ptr_t;
#endif
