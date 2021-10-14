/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/

#ifndef PHIST_EXCEPTION_HPP
#define PHIST_EXCEPTION_HPP

#ifndef DOXYGEN
#include <exception>
#include "phist_tools.h"
#endif

namespace phist
{

//! \ingroup cxx_bindings
//!@{

  //! class for reporting an error (negative) or warning (positive) iflag on return
  class Exception : public std::exception
  {
    public:
    
    //! constructor
    Exception(int iflag) : std::exception(), iflag_(iflag) {}
 
   //! destructor
   ~Exception(){}
   
   //!
   const char* what() const noexcept {return phist_retcode2str(iflag_);}
   
   //!
   inline int iflag() const noexcept {return iflag_;}

  protected:
  
    //!
    int iflag_;
  
  };
//!@}
}

#endif
