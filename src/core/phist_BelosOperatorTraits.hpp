/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
//! \file phist_BelosOperatorTraits.hpp
//! \brief Implementation of the Belos::OperatorTraits for phist_linearOp_t.

#ifndef PHIST_BELOS_OPERATOR_TRAITS_HPP
#define PHIST_BELOS_OPERATOR_TRAITS_HPP

#include "phist_config.h"

#ifdef PHIST_HAVE_BELOS

#ifndef DOXYGEN

#include "phist_operator.h"
#include "phist_ScalarTraits.hpp"
#include "phist_BelosMV.hpp"
#include "BelosTypes.hpp"
#include <BelosOperatorTraits.hpp>

#endif //DOXYGEN

// this file is mostly copied from the Belos Tpetra adapter implementation in Trilinos 11.2.4
//! namespace for Belos
namespace Belos {

  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Belos::OperatorTraits for phist_linearOp_t.
  //
  ////////////////////////////////////////////////////////////////////

/// \brief Partial specialization of OperatorTraits for phist_linearOp_t.
template <typename ST>
class OperatorTraits <ST, ::phist::BelosMV<ST>, typename phist::ScalarTraits<ST>::linearOp_t >
{
  public:

    typedef typename phist::ScalarTraits<ST>::linearOp_t phist_linearOp_t;
    typedef typename phist::ScalarTraits<ST>::mvec_t phist_mvec_t;
    typedef phist::BelosMV<ST> MV;

    static void 
    Apply (const phist_linearOp_t& Op, 
     const MV& X,
     MV& Y,
     ETrans trans=NOTRANS)
    {
      TEUCHOS_TEST_FOR_EXCEPTION(trans!=NOTRANS,std::invalid_argument,
          "Belos::OperatorTraits<ST,MV,phist_linearOp_t>:: Apply: only implemented for trans=NOTRANS up to now.");
      int iflag=0;
      ST alpha = phist::ScalarTraits<ST>::one();
      ST beta = phist::ScalarTraits<ST>::zero();
      Op.apply(alpha,Op.A,X.get(), beta, Y.get(),&iflag);
    }

    static bool
    HasApplyTranspose (const phist_linearOp_t& Op)
    {
      return false;
    }
};

} // end of Belos namespace 
#endif /* PHIST_HAVE_BELOS */
#endif // end of file BELOS_PHIST_ADAPTER_HPP
