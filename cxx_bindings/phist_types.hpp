/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/

#include "phist_typedefs.h"
#include "phist_ScalarTraits.hpp"

namespace phist {

template<typename ST>
class types
{
  public:
  
    typedef ScalarTraits<ST> st;
    typedef typename ScalarTraits<ST>::magn_t MT;
    typedef ScalarTraits<MT> mt;
    
    typedef typename st::mvec_t* mvec_ptr;
    typedef typename st::mvec_t const* const_mvec_ptr;

    typedef typename st::sdMat_t* sdMat_ptr;
    typedef typename st::sdMat_t const* const_sdMat_ptr;

    typedef typename st::sparseMat_t* sparseMat_ptr;
    typedef typename st::sparseMat_t const* const_sparseMat_ptr;

    typedef typename st::linearOp_t* linearOp_ptr;
    typedef typename st::linearOp_t const* const_linearOp_ptr;
};

}
