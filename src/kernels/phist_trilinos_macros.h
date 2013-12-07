#include "Teuchos_Utils.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#ifndef TRY_CATCH
#define TRY_CATCH(func,_FLAG_) \
{bool _success=true; \
try {func;} TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,_success); \
if (!_success) {_FLAG_=-77; return;}}
#endif
