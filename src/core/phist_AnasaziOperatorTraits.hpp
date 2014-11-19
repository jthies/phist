#ifndef PHIST_ANASAZI_OPERATOR_TRAITS_HPP
#define PHIST_ANASAZI_OPERATOR_TRAITS_HPP

#include "phist_config.h"

/* needs to be included before system headers for some intel compilers+mpi */
#ifndef DOXYGEN

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#endif

#if defined(PHIST_HAVE_ANASAZI) && defined(PHIST_HAVE_BELOS)

#include "AnasaziOperatorTraits.hpp"
#include "phist_BelosOperatorTraits.hpp"

namespace Anasazi 
{

//! class to allow using phist operators in Anasazi 
//! (block Eigensolvers from Trilinos). Simply inherit
//! the traits from Belos (block iterative linear solver package),
//! as they are essentially the same.
template <class Scalar, class MV>
class OperatorTraits <Scalar, MV, typename phist::ScalarTraits<Scalar>::op_t >
  : public ::Belos::OperatorTraits<Scalar, MV, typename phist::ScalarTraits<Scalar>::op_t >
{
};

}// namespace
#endif
#endif
