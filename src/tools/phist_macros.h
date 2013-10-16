#ifndef NO_INCLUDES_IN_HEADERS
#include "phist_tools.h"
#include <stdio.h>
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#endif

// these can be passed to PHIST_OUT(FLAG,...)
// to get a certain amount of coherence in the 
// way screen output is handled. The PHIST_OUTLEV
// macro defines which messages are printed and
// which aren't.
#define PHIST_ERROR   0
#define PHIST_WARNING 1
#define PHIST_INFO 2
#define PHIST_VERBOSE 3
#define PHIST_DEBUG 4
#define PHIST_TRACE 5

#ifndef PHIST_OUTLEV
#define PHIST_OUTLEV PHIST_WARNIUNG
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
#ifdef __cplusplus
#define PHIST_CHK_IERR(func,_FLAG_) \
try {func; if (_FLAG_!=_PHIST_SUCCESS_) { \
PHIST_OUT(PHIST_ERROR,"Error code %d (%s) returned from call %s\n(file %s, line %d)",\
(_FLAG_),(phist_retcode2str(_FLAG_)),(#func),(__FILE__),(__LINE__)); return;}} \
catch (std::exception e) {PHIST_OUT(PHIST_ERROR,"Exception caught in call %s (%s)\n(file %s, line %d)",\
(#func),e.what(),(__FILE__),(__LINE__)); (_FLAG_)=-77; return;} \
catch (std::string s) {PHIST_OUT(PHIST_ERROR,"Exception caught in call %s (%s)\n(file %s, line %d)",\
(#func),s.c_str(),(__FILE__),(__LINE__)); (_FLAG_)=-77; return;} \
catch (int iexc) {PHIST_OUT(PHIST_ERROR,"int Exception caught in call %s (value %d)\n(file %s, line %d)",\
(#func),iexc,(__FILE__),(__LINE__)); (_FLAG_)=-77; return;} \
catch (...) {PHIST_OUT(PHIST_ERROR,"unknown Exception caught in call %s\n(file %s, line %d)",\
(#func),(__FILE__),(__LINE__)); (_FLAG_)=-77; return;}
#else
#define PHIST_CHK_IERR(func,_FLAG_) \
{func; if (_FLAG_!=_PHIST_SUCCESS_) { \
PHIST_OUT(PHIST_ERROR,"Error code %d (%s) returned from call %s\n(file %s, line %d)",\
(_FLAG_),(phist_retcode2str(_FLAG_)),(#func),(__FILE__),(__LINE__)); return;}}
#endif
#endif

//! this macro is deprecated, PHIST_CHK_IERR should now be used.
#ifndef _PHIST_ERROR_HANDLER_
#define _PHIST_ERROR_HANDLER_(func,_FLAG_) PHIST_CHK_IERR(func,_FLAG_)
#endif


//! like PHIST_CHK_IERR, but returns ierr (to be used in int functions returning an error code)
#ifndef PHIST_ICHK_IERR
#define PHIST_ICHK_IERR(func,_FLAG_) \
{func; if (_FLAG_!=_PHIST_SUCCESS_) { \
PHIST_OUT(PHIST_ERROR,"Error code %d (%s) returned from call %s\n(file %s, line %d)",\
(_FLAG_),(phist_retcode2str(_FLAG_)),(#func),(__FILE__),(__LINE__)); return _FLAG_;}}
#endif

//! checks return value of an int function for non-zero value, assigns it to _FLAG_,
//! prints an error message and returns if non-zero (to be used in void functions)
#ifndef PHIST_CHK_IRET
#define PHIST_CHK_IRET(func,_FLAG_) \
{FLAG=func; if (_FLAG_!=_PHIST_SUCCESS_) { \
PHIST_OUT(PHIST_ERROR,"Error code %d (%s) returned from call %s\n(file %s, line %d)",\
(_FLAG_),(phist_retcode2str(_FLAG_)),(#func),(__FILE__),(__LINE__)); return;}}
#endif

//! like PHIST_CHK_IERR, but returns ierr (to be used in int functions returning an error code)
#ifndef PHIST_ICHK_IRET
#define PHIST_ICHK_IRET(func,_FLAG_) \
{FLAG=func; if (_FLAG_!=_PHIST_SUCCESS_) { \
PHIST_OUT(PHIST_ERROR,"Error code %d (%s) returned from call %s\n(file %s, line %d)",\
(_FLAG_),(phist_retcode2str(_FLAG_)),(#func),(__FILE__),(__LINE__)); return _FLAG_;}}
#endif

#ifdef PHIST_HAVE_MPI
#define PHIST_OUT(level,msg, ...) {\
        if(PHIST_OUTLEV >= level) {\
                FILE* out= (level<=PHIST_WARNING)? stderr:stdout;\
                int me,ini,fini;\
                MPI_Initialized(&ini); \
                if (ini) MPI_Finalized(&fini); \
                if (ini && (!fini)) { \
                MPI_Comm_rank(MPI_COMM_WORLD,&me);\
                } else {me=0;}\
                fprintf(out,"PE%d: "msg"\n",me,##__VA_ARGS__);\
                fflush(out);\
        }\
}
#else
#define PHIST_OUT(level,msg, ...) {\
        if(PHIST_OUTLEV >= level) {\
                FILE* out= (level<=PHIST_WARNING)? stderr:stdout;
                fprintf(out,msg"\n",##__VA_ARGS__);\
                fflush(out);\
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

// function tracer can only be used in C++ code
#ifdef __cplusplus

#ifndef PHIST_FCN_TRACER
#define PHIST_FCN_TRACER

class FcnTracer
  {
  public:
  
  FcnTracer(const char* fcn) : fcn_(fcn)
    {
    PHIST_OUT(PHIST_TRACE,"PHIST ENTER %s\n",fcn_.c_str());
    }

  ~FcnTracer()
    { 
    PHIST_OUT(PHIST_TRACE,"PHIST LEAVE %s\n",fcn_.c_str()); 
    }

  std::string fcn_;
  };

#endif

#if (PHIST_OUTLEV>=PHIST_TRACE)
#define ENTER_FCN(s) FcnTracer FnT(s);
#else
#define ENTER_FCN(s)
#endif

#endif // cplusplus
