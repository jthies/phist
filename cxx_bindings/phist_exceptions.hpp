/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/

#include <exception>

namespace phist
{
  //! class for reporting an error (negative iflag on return)
  class Erro: public std::exception<const char*>
  {
    public:
    
    //!
    Error(const char* message) : std::exception(message) {}
  };

  //! class for reporting a warning (positive iflag on return)
  class Erro: public std::exception<const char*>
  {
    public:
    
    //!
    Warning(const char* message) : std::exception(message) {}
  };
}
