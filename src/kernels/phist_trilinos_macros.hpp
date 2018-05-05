/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef DOXYGEN
#include "Teuchos_Utils.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#ifndef PHIST_TRY_CATCH
#define PHIST_TRY_CATCH(func,_FLAG_) \
{bool _success=true; \
try {func;} TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,_success); \
if (!_success) {_FLAG_=-77; return;}}
#endif
#endif /* DOXYGEN */
