#ifndef PHIST_MACROS_H
#define PHIST_MACROS_H

#include "phist_tools.h"

#define _PHIST_SUCCESS_ 0
#define _PHIST_FUNCTIONAL_ERROR_ -1
#define _PHIST_CAUGHT_EXCEPTION_ -77
#define _PHIST_BAD_CAST_ -88
#define _PHIST_NOT_IMPLEMENTED_ -99

#ifndef _PHIST_ERROR_HANDLER_
#define _PHIST_ERROR_HANDLER_(func,_FLAG_) \
{func; if (_FLAG_!=_PHIST_SUCCESS_) { \
fprintf(stderr,"Error code %d (%s) returned from call %s\n(file %s, line %d)",\
(_FLAG_),(phist_retcode2str(_FLAG_)),(#func),(__FILE__),(__LINE__)); return;}}
#endif


#ifndef _PHIST_TEST_HANDLER_
#define _PHIST_TEST_HANDLER_(func,_FLAG_,_FLAG_EXPECT_) \
{func; if (_FLAG_!=_FLAG_EXPECT_) { \
fprintf(stderr,"Unexpected error code %d (%s) returned from call %s\n(file %s, line %d)",\
(_FLAG_),(phist_retcode2str(_FLAG_)),(#func),(__FILE__),(__LINE__));}}
#endif

#endif

// if a phist_gen_*.h file has been included, define
// macros to build function names and such.
#ifdef _TP_

// _MT_ (magnitude type) for real numbers
#ifndef _IS_COMPLEX_
#ifdef _MT_
#undef _MT_
#endif

#define _MT_ _ST_
#endif

#ifndef _IS_COMPLEX_
#ifdef _Complex_I
#undef _Complex_I
#endif
#define _Complex_I (_ST_)0.0
#endif

#endif
