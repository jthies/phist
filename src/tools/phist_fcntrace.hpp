#ifndef NO_INCLUDES_IN_HEADERS
#include "phist_tools.h"
#include <stdio.h>
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#endif

#ifndef DEBUG
#define DEBUG 0
#endif

#ifndef _PHIST_RETURN_TYPES
#define _PHIST_RETURN_TYPES

#define _PHIST_SUCCESS_ 0
#define _PHIST_FUNCTIONAL_ERROR_ -1
#define _PHIST_CAUGHT_EXCEPTION_ -77
#define _PHIST_BAD_CAST_ -88
#define _PHIST_NOT_IMPLEMENTED_ -99
#endif

//! checks an ierr flag passed to a void function for non-zero value, assigns it to _FLAG_,
//! prints an error message and returns if non-zero (to be used in void functions)
#ifndef PHIST_CHK_IERR
#define PHIST_CHK_IERR(func,_FLAG_) \
{func; if (_FLAG_!=_PHIST_SUCCESS_) { \
PHIST_OUT(0,"Error code %d (%s) returned from call %s\n(file %s, line %d)",\
(_FLAG_),(phist_retcode2str(_FLAG_)),(#func),(__FILE__),(__LINE__)); return;}}
#endif

//! this macro is deprecated, PHIST_CHK_IERR should now be used.
#ifndef _PHIST_ERROR_HANDLER_
#define _PHIST_ERROR_HANDLER_(func,_FLAG_) PHIST_CHK_IERR(func,_FLAG_)
#endif


//! like PHIST_CHK_IERR, but returns ierr (to be used in int functions returning an error code)
#ifndef PHIST_ICHK_IERR
#define PHIST_ICHK_IERR(func,_FLAG_) \
{func; if (_FLAG_!=_PHIST_SUCCESS_) { \
PHIST_OUT(0,"Error code %d (%s) returned from call %s\n(file %s, line %d)",\
(_FLAG_),(phist_retcode2str(_FLAG_)),(#func),(__FILE__),(__LINE__)); return _FLAG_;}}
#endif

//! checks return value of an int function for non-zero value, assigns it to _FLAG_,
//! prints an error message and returns if non-zero (to be used in void functions)
#ifndef PHIST_CHK_IRET
#define PHIST_CHK_IRET(func,_FLAG_) \
{FLAG=func; if (_FLAG_!=_PHIST_SUCCESS_) { \
PHIST_OUT(0,"Error code %d (%s) returned from call %s\n(file %s, line %d)",\
(_FLAG_),(phist_retcode2str(_FLAG_)),(#func),(__FILE__),(__LINE__)); return;}}
#endif

//! like PHIST_CHK_IERR, but returns ierr (to be used in int functions returning an error code)
#ifndef PHIST_ICHK_IRET
#define PHIST_ICHK_IRET(func,_FLAG_) \
{FLAG=func; if (_FLAG_!=_PHIST_SUCCESS_) { \
PHIST_OUT(0,"Error code %d (%s) returned from call %s\n(file %s, line %d)",\
(_FLAG_),(phist_retcode2str(_FLAG_)),(#func),(__FILE__),(__LINE__)); return _FLAG_;}}
#endif

#ifdef PHIST_HAVE_MPI
#define PHIST_OUT(level,msg, ...) {\
        if(DEBUG >= level) {\
                int __me,__ini,__fini;\
                MPI_Initialized(&__ini); \
                if (__ini) MPI_Finalized(&__fini); \
                if (__ini && (!__fini)) { \
                MPI_Comm_rank(MPI_COMM_WORLD,&__me);\
                } else {__me=0;}\
                fprintf(stderr,"PE%d: "msg"\n",__me,##__VA_ARGS__);\
                fflush(stderr);\
        }\
}
#else
#define PHIST_OUT(level,msg, ...) {\
        if(DEBUG >= level) {\
                fprintf(stderr,msg"\n",##__VA_ARGS__);\
                fflush(stderr);\
        }\
}
#endif

#ifdef MAX
#undef MAX
#endif

#define MAX(a,b) (a)<(b)?(b):(a);

#ifdef MIN
#undef MIN
#endif

#define MIN(a,b) (b)<(a)?(b):(a);
extern "C" {

#ifndef PHIST_FCN_TRACER
#define PHIST_FCN_TRACER

class FcnTracer
  {
  public:
  
  FcnTracer(const char* fcn) : fcn_(fcn)
    {
    PHIST_OUT(0,"PHIST ENTER %s\n",fcn_.c_str());
    }

  ~FcnTracer()
    { 
    PHIST_OUT(0,"PHIST LEAVE %s\n",fcn_.c_str()); 
    }

  std::string fcn_;
  };

#endif

#if (DEBUG>=9)
#define ENTER_FCN(s) FcnTracer FnT(s);
#else
#define ENTER_FCN(s)
#endif
