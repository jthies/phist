/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef ANASAZI_PHIST_ADAPTER_HPP
#define ANASAZI_PHIST_ADAPTER_HPP

#include "phist_config.h"

#if defined(PHIST_HAVE_ANASAZI)&&defined(PHIST_HAVE_BELOS)

#include "AnasaziMultiVecTraits.hpp"
#include "phist_BelosMV.hpp"
#include "Belos_PhistAdapter.hpp"

using namespace phist;

namespace Anasazi 
{

  template<typename ST>
  class MultiVecTraits<ST,BelosMV<ST> >: public Belos::MultiVecTraits<ST,BelosMV<ST> >
  {
  };
}// namespace Anasazi


#endif
#endif
