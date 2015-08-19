#ifndef PHIST_JADA_TEST_WITH_OPTS_HPP
#define PHIST_JADA_TEST_WITH_OPTS_HPP

#include "phist_jadaOpts.h"

//! base class for tests that require a jadaOpts_t,
//! it just provides such an object and initializes
//! it in its SetUp routine.
class JadaTestWithOpts
{

  public:
  
  phist_jadaOpts_t jadaOpts_;
  
  virtual void SetUp()
  {
    phist_jadaOpts_setDefaults(&jadaOpts_);
  }

};


#endif
