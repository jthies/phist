/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
//! \file phist_AnasaziOperatorTraits.hpp

#ifndef PHIST_ANASAZI_OPERATOR_TRAITS_HPP
#define PHIST_ANASAZI_OPERATOR_TRAITS_HPP

#include "phist_config.h"

#if defined(PHIST_HAVE_ANASAZI) && defined(PHIST_HAVE_BELOS)

#ifndef DOXYGEN

#include "AnasaziOperatorTraits.hpp"
#include "phist_BelosOperatorTraits.hpp"

#endif /* DOXYGEN */
//! namespace for Anasazi
namespace Anasazi 
{

//! class to allow using phist operators in Anasazi 
//! (block Eigensolvers from Trilinos). Simply inherit
//! the traits from Belos (block iterative linear solver package),
//! as they are essentially the same.
template <class Scalar, class MV>
class OperatorTraits <Scalar, MV, typename phist::ScalarTraits<Scalar>::linearOp_t >
  : public ::Belos::OperatorTraits<Scalar, MV, typename phist::ScalarTraits<Scalar>::linearOp_t >
{
};

}// namespace
#endif
#endif
