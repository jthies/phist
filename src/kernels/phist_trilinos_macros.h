#include "Teuchos_Utils.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#ifndef CAST_PTR_FROM_VOID
#define CAST_PTR_FROM_VOID(_TYPE_,_PTR_,_VPTR_,_FLAG_) \
_TYPE_ *_PTR_ = NULL; \
_PTR_=(_TYPE_*)(_VPTR_); \
if (_PTR_==NULL) {_FLAG_=-88; PHIST_OUT(PHIST_ERROR,"bad cast in file %s, line %d.",\
__FILE__,__LINE__); return;}
#endif

#ifndef TRY_CATCH
#define TRY_CATCH(func,_FLAG_) \
{bool _success=true; \
try {func;} TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,_success); \
if (!_success) {_FLAG_=-77; return;}}
#endif
