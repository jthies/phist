#include "phist_tools.h"

#ifndef _PHIST_RETURN_TYPES
#define _PHIST_RETURN_TYPES

#define _PHIST_SUCCESS_ 0
#define _PHIST_FUNCTIONAL_ERROR_ -1
#define _PHIST_CAUGHT_EXCEPTION_ -77
#define _PHIST_BAD_CAST_ -88
#define _PHIST_NOT_IMPLEMENTED_ -99
#endif

#ifndef _PHIST_ERROR_HANDLER_
#define _PHIST_ERROR_HANDLER_(func,_FLAG_) \
{func; if (_FLAG_!=_PHIST_SUCCESS_) { \
fprintf(stderr,"Error code %d (%s) returned from call %s\n(file %s, line %d)",\
(_FLAG_),(phist_retcode2str(_FLAG_)),(#func),(__FILE__),(__LINE__)); return;}}

#endif

