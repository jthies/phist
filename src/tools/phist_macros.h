#ifndef PHIST_MACROS_H
#define PHIST_MACROS_H

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#ifndef NO_INCLUDES_IN_HEADERS
#include "phist_tools.h"
#ifdef __cplusplus
#include <cstdio>
#else
#include <stdio.h>
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
#define PHIST_OUTLEV PHIST_INFO
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
                fflush(out);\ }\ }
#endif

#ifdef PHIST_HAVE_MPI
#define PHIST_SOUT(level,msg, ...) {\
        if(PHIST_OUTLEV >= level) {\
                FILE* out= (level<=PHIST_WARNING)? stderr:stdout;\
                int me,ini,fini;\
                MPI_Initialized(&ini); \
                if (ini) MPI_Finalized(&fini); \
                if (ini && (!fini)) { \
                MPI_Comm_rank(MPI_COMM_WORLD,&me);\
                } else {me=0;}\
                if(me==0){\
                fprintf(out,msg"\n",##__VA_ARGS__);\
                fflush(out);}\
        }\
}
#else
#define PHIST_SOUT(level,msg, ...) {\
        if(PHIST_OUTLEV >= level) {\
                FILE* out= (level<=PHIST_WARNING)? stderr:stdout;
                fprintf(out,msg"\n",##__VA_ARGS__);\
                fflush(out);\
        }\
}
#endif


// return types
#define PHIST_SUCCESS 0
#define PHIST_FUNCTIONAL_ERROR -1
#define PHIST_CAUGHT_EXCEPTION -77
#define PHIST_BAD_CAST -88
#define PHIST_NOT_IMPLEMENTED -99

//! checks an ierr flag passed to a void function for non-zero value, assigns it to FLAG,
//! prints an error message and returns if non-zero (to be used in void functions)
#ifdef __cplusplus
#define PHIST_CHK_IERR(func,FLAG) {\
try {func; if (FLAG!=PHIST_SUCCESS) { \
PHIST_OUT(PHIST_ERROR,"Error code %d (%s) returned from call %s\n(file %s, line %d)",\
(FLAG),(phist_retcode2str(FLAG)),(#func),(__FILE__),(__LINE__)); return;}} \
catch (std::exception e) {PHIST_OUT(PHIST_ERROR,"Exception caught in call %s (%s)\n(file %s, line %d)",\
(#func),e.what(),(__FILE__),(__LINE__)); (FLAG)=-77; return;} \
catch (std::string s) {PHIST_OUT(PHIST_ERROR,"Exception caught in call %s (%s)\n(file %s, line %d)",\
(#func),s.c_str(),(__FILE__),(__LINE__)); (FLAG)=-77; return;} \
catch (int iexc) {PHIST_OUT(PHIST_ERROR,"int Exception caught in call %s (value %d)\n(file %s, line %d)",\
(#func),iexc,(__FILE__),(__LINE__)); (FLAG)=-77; return;} \
catch (...) {PHIST_OUT(PHIST_ERROR,"unknown Exception caught in call %s\n(file %s, line %d)",\
(#func),(__FILE__),(__LINE__)); (FLAG)=-77; return;}}
#else
#define PHIST_CHK_IERR(func,FLAG) {\
{func; if (FLAG!=PHIST_SUCCESS) { \
PHIST_OUT(PHIST_ERROR,"Error code %d (%s) returned from call %s\n(file %s, line %d)",\
(FLAG),(phist_retcode2str(FLAG)),(#func),(__FILE__),(__LINE__)); return;}}}
#endif

//! checks an ierr flag passed to a void function for negative value, assigns it to FLAG,
//! prints an error message and returns if non-zero (to be used in void functions)
#ifndef PHIST_CHK_NEG_IERR
#ifdef __cplusplus
#define PHIST_CHK_NEG_IERR(func,FLAG) {\
try {func; if (FLAG < PHIST_SUCCESS) { \
PHIST_OUT(PHIST_ERROR,"Error code %d (%s) returned from call %s\n(file %s, line %d)",\
(FLAG),(phist_retcode2str(FLAG)),(#func),(__FILE__),(__LINE__)); return;}} \
catch (std::exception e) {PHIST_OUT(PHIST_ERROR,"Exception caught in call %s (%s)\n(file %s, line %d)",\
(#func),e.what(),(__FILE__),(__LINE__)); (FLAG)=-77; return;} \
catch (std::string s) {PHIST_OUT(PHIST_ERROR,"Exception caught in call %s (%s)\n(file %s, line %d)",\
(#func),s.c_str(),(__FILE__),(__LINE__)); (FLAG)=-77; return;} \
catch (int iexc) {PHIST_OUT(PHIST_ERROR,"int Exception caught in call %s (value %d)\n(file %s, line %d)",\
(#func),iexc,(__FILE__),(__LINE__)); (FLAG)=-77; return;} \
catch (...) {PHIST_OUT(PHIST_ERROR,"unknown Exception caught in call %s\n(file %s, line %d)",\
(#func),(__FILE__),(__LINE__)); (FLAG)=-77; return;}}
#else
#define PHIST_CHK_NEG_IERR(func,FLAG) {\
{func; if (FLAG<PHIST_SUCCESS) { \
PHIST_OUT(PHIST_ERROR,"Error code %d (%s) returned from call %s\n(file %s, line %d)",\
(FLAG),(phist_retcode2str(FLAG)),(#func),(__FILE__),(__LINE__)); return;}}}
#endif
#endif


//! like PHIST_CHK_IERR, but returns ierr (to be used in int functions returning an error code)
#define PHIST_ICHK_IERR(func,FLAG) {\
{func; if (FLAG!=PHIST_SUCCESS) { \
PHIST_OUT(PHIST_ERROR,"Error code %d (%s) returned from call %s\n(file %s, line %d)",\
(FLAG),(phist_retcode2str(FLAG)),(#func),(__FILE__),(__LINE__)); return FLAG;}}}

#ifdef MAX
#undef MAX
#endif
#define MAX(a,b) ((a)<(b)?(b):(a))

#ifdef MIN
#undef MIN
#endif
#define MIN(a,b) ((b)<(a)?(b):(a))

#if PHIST_OUTLEV>=PHIST_DEBUG
#define PHIST_DEB(msg, ...) PHIST_OUT(PHIST_DEBUG,msg,##__VA_ARGS__);
#else
#define PHIST_DEB(msg, ...)
#endif

#ifdef __cplusplus
# include "phist_fcntrace.hpp"
# ifdef PHIST_TIMEMONITOR
#   include <Teuchos_TimeMonitor.hpp>
#   if (PHIST_OUTLEV>=PHIST_TRACE) || (LIKWID_PERFMON)
#     define ENTER_FCN(s) FcnTracer YouCantHaveMultiple_ENTER_FCN_StatementsInOneScope(s);\
                          static Teuchos::RCP<Teuchos::Time> TimeFrom_ENTER_FCN; \
                          if( TimeFrom_ENTER_FCN.is_null() ) \
                              TimeFrom_ENTER_FCN = Teuchos::TimeMonitor::getNewTimer(s); \
                          Teuchos::TimeMonitor TimeMonFrom_ENTER_FCN( *TimeFrom_ENTER_FCN );
#   else
#     define ENTER_FCN(s) static Teuchos::RCP<Teuchos::Time> TimeFrom_ENTER_FCN; \
                          if( TimeFrom_ENTER_FCN.is_null() ) \
                              TimeFrom_ENTER_FCN = Teuchos::TimeMonitor::getNewTimer(s); \
                          Teuchos::TimeMonitor TimeMonFrom_ENTER_FCN( *TimeFrom_ENTER_FCN );
#   endif
# else
#   if (PHIST_OUTLEV>=PHIST_TRACE) || (LIKWID_PERFMON)
#     define ENTER_FCN(s) FcnTracer YouCantHaveMultiple_ENTER_FCN_StatementsInOneScope(s);
#   else
#     define ENTER_FCN(s)
#   endif
# endif
#else
# define ENTER_FCN(s)
#endif

#ifndef TOUCH
#define TOUCH(x) (void)(x);
#endif

#ifndef CAST_PTR_FROM_VOID
#define CAST_PTR_FROM_VOID(_TYPE_,_PTR_,_VPTR_,_FLAG_) \
_TYPE_ *_PTR_ = NULL; \
_PTR_=(_TYPE_*)(_VPTR_); \
if (_PTR_==NULL) {_FLAG_=-88; PHIST_OUT(PHIST_ERROR,"bad cast in file %s, line %d.",\
__FILE__,__LINE__); return;}
#endif

#endif
