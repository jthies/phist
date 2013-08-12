#include "Teuchos_Utils.hpp"
#include "Teuchos_StandardCatchMacros.hpp"


#ifndef _PRINT_MSG_
#define _PRINT_MSG_(s) {std::cerr << (s) << std::endl;}
#endif

#ifndef _DEBUG_
#define _DEBUG_(s) {std::cerr<<#s<<" = "<<(s)<<std::endl;}
#endif

#ifndef _CAST_PTR_FROM_VOID_
#define _CAST_PTR_FROM_VOID_(TYPE,_PTR_,_VPTR_,_FLAG_) \
TYPE *_PTR_ = NULL; \
_PTR_=(TYPE*)(_VPTR_); \
if (_PTR_==NULL) {_FLAG_=-88; return;}
#endif

#ifndef _CHECKZERO
#define _CHECKZERO(func,_FLAG_) \
{*ierr=func; if (_FLAG_) { \
std::string msg="Error code "+Teuchos::toString(_FLAG_)+" returned from call "\
+"'"+#func+"'\n"+\
"(file "+__FILE__+", line "+Teuchos::toString(__LINE__)+")"; \
_PRINT_MSG_(msg); return;}}
#endif

#ifndef _TRY_CATCH_
#define _TRY_CATCH_(func,_FLAG_) \
{bool _success=true; \
try {func;} TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,_success); \
if (!_success) {_FLAG_=-77; return;}}
#endif
