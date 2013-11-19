#ifndef PHIST_BELOS_OPERATOR_TRAITS_HPP
#define PHIST_BELOS_OPERATOR_TRAITS_HPP

#include "phist_rcp_helpers.hpp"
#include "phist_operator.h"
#include "phist_ScalarTraits.hpp"
#include "BelosTypes.hpp"
#include <BelosOperatorTraits.hpp>

// this file is mostly copied from the Belos Tpetra adapter implementation in Trilinos 11.2.4

namespace Belos {

  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Belos::OperatorTraits for phist_op_t.
  //
  ////////////////////////////////////////////////////////////////////

  /// \brief Partial specialization of OperatorTraits for phist_op_t.
  /// Note: this is only going to work if MV is the class used by the kernel
  /// library, for instance, you can't compile phist with PHIST_KERNEL_LIB=
  /// tpetra and then use ghost vectors in Belos.
  template <class Scalar, class MV> 
  class OperatorTraits <Scalar, MV, typename phist::ScalarTraits<Scalar>::op_t >
  {
  public:

    typedef typename phist::ScalarTraits<Scalar>::op_t phist_op_t;
    typedef typename phist::ScalarTraits<Scalar>::mvec_t phist_mvec_t;

    static void 
    Apply (const phist_op_t& Op, 
	   const MV& X,
	   MV& Y,
	   ETrans trans=NOTRANS)
    {
    TEUCHOS_TEST_FOR_EXCEPTION(trans!=NOTRANS,std::invalid_argument,
          "Belos::OperatorTraits<Scalar,MV,phist_op_t>:: Apply: only implemented for trans=NOTRANS up to now.");
    int ierr;
    Scalar alpha = phist::ScalarTraits<Scalar>::one();
    Scalar beta = phist::ScalarTraits<Scalar>::zero();
    Op.apply(alpha,Op.A,phist::ref2ptr(X), beta, phist::ref2ptr(Y),&ierr);
    }

    static bool
    HasApplyTranspose (const phist_op_t& Op)
    {
      return false;
    }
  };

} // end of Belos namespace 

#endif 
// end of file BELOS_GHOST_ADAPTER_HPP
