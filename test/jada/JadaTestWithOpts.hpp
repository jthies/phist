/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef PHIST_JADA_TEST_WITH_OPTS_HPP
#define PHIST_JADA_TEST_WITH_OPTS_HPP

#include "phist_jadaOpts.h"

//! base class for tests that require a jadaOpts_t,
//! it just provides such an object and initializes
//! it in its SetUp routine.
class JadaTestWithOpts
{

  public:
  
  phist_jadaOpts jadaOpts_;
  
  virtual void SetUp()
  {
    phist_jadaOpts_setDefaults(&jadaOpts_);
  }

};


#endif
