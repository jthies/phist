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
//! programmer from writing down all the macros like TYPE(...)
//! all the time. The programmer must take care that this file
//! is included once per scope only, otherwise the compiler will 
//! complain.
typedef ::phist::ScalarTraits< _ST_ > st;
typedef st::scalar_t ST;
typedef st::magn_t MT;
typedef st::mvec_t* mvec_ptr_t;
typedef const st::mvec_t* const_mvec_ptr_t;
typedef st::crsMat_t* crsMat_ptr_t;
typedef const st::crsMat_t* const_crsMat_ptr_t;
typedef st::sdMat_t* sdMat_ptr_t;
typedef const st::sdMat_t* const_sdMat_ptr_t;
typedef st::op_t* op_ptr_t;
typedef const st::op_t* const_op_ptr_t;

typedef ::phist::ScalarTraits<MT> mt;
