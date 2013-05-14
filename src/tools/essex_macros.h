#ifndef ESSEX_MACROS_H
#define ESSEX_MACROS_H

#include "essex_tools.h"

#define _ESSEX_SUCCESS_ 0
#define _ESSEX_FUNCTIONAL_ERROR_ -1
#define _ESSEX_CAUGHT_EXCEPTION_ -77
#define _ESSEX_BAD_CAST_ -88
#define _ESSEX_NOT_IMPLEMENTED_ -99

#ifndef _ESSEX_ERROR_HANDLER_
#define _ESSEX_ERROR_HANDLER_(func,_FLAG_) \
{func; if (_FLAG_!=_ESSEX_SUCCESS_) { \
fprintf(stderr,"Error code %d (%s) returned from call %s\n(file %s, line %d)",\
_FLAG_,essex_retcode2str(_FLAG_),#func,__FILE__,__LINE__); return;}}
#endif


#ifndef _ESSEX_TEST_HANDLER_
#define _ESSEX_TEST_HANDLER_(func,_FLAG_,_FLAG_EXPECT_) \
{func; if (_FLAG_!=_FLAG_EXPECT_) { \
fprintf(stderr,"Unexpected error code %d (%s) returned from call %s\n(file %s, line %d)",\
_FLAG_,essex_retcode2str(_FLAG_),#func,__FILE__,__LINE__);}}
#endif
#endif
