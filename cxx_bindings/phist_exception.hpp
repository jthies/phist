/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
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

//! \addtogroup cxx_interface
//!@{

  //! class for reporting an error (negative) or warning (positive) iflag on return
  class Exception : public std::exception
  {
    public:
    
    //!
    Exception(int iflag) : std::exception(), iflag_(iflag) {}
 
   //!
   ~Exception(){}
   
   //!
   const char* what() const noexcept {return phist_retcode2str(iflag_);}
   
   //!
   inline int iflag() const noexcept {return iflag_;}

  protected:
  
    //!
    int iflag_;
  
  };
//@!}
}

#endif
