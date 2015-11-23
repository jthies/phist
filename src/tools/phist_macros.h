#ifndef PHIST_MACROS_H
#define PHIST_MACROS_H

#include "phist_config.h"

#ifndef DOXYGEN

/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "phist_tools.h"
#ifdef __cplusplus
#include <cstdio>
#else
#include <stdio.h>
#endif
#endif

#include "phist_defs.h"

#ifdef PHIST_HAVE_MPI
#define PHIST_OUT(level,msg, ...) {\
        if(PHIST_OUTLEV >= level) {\
                FILE* PHIST_OUT_out= (level<=PHIST_WARNING)? stderr:stdout;\
                int PHIST_OUT_me,PHIST_OUT_np,PHIST_OUT_ini,PHIST_OUT_fini;\
                MPI_Initialized(&PHIST_OUT_ini); \
                if (PHIST_OUT_ini) MPI_Finalized(&PHIST_OUT_fini); \
                if (PHIST_OUT_ini && (!PHIST_OUT_fini)) { \
                MPI_Comm_rank(MPI_COMM_WORLD,&PHIST_OUT_me);\
                MPI_Comm_size(MPI_COMM_WORLD,&PHIST_OUT_np);\
                } else {PHIST_OUT_me=0; PHIST_OUT_np=1;}\
                if (PHIST_OUT_np>1) \
                {fprintf(PHIST_OUT_out,"PE%d: " msg,PHIST_OUT_me,##__VA_ARGS__);}\
                else \
                {fprintf(PHIST_OUT_out,msg,##__VA_ARGS__);}\
                fflush(PHIST_OUT_out);\
        }\
}
#else
#define PHIST_OUT(level,msg, ...) {\
        if(PHIST_OUTLEV >= level) {\
                FILE* PHIST_OUT_out= (level<=PHIST_WARNING)? stderr:stdout;\
                fprintf(PHIST_OUT_out,msg,##__VA_ARGS__);\
                fflush(PHIST_OUT_out);\
                }\
        }
#endif

#ifdef PHIST_HAVE_MPI
#define PHIST_SOUT(level,msg, ...) {\
        if(PHIST_OUTLEV >= level) {\
                FILE* PHIST_OUT_out= (level<=PHIST_WARNING)? stderr:stdout;\
                int PHIST_OUT_me,PHIST_OUT_ini,PHIST_OUT_fini;\
                MPI_Initialized(&PHIST_OUT_ini); \
                if (PHIST_OUT_ini) MPI_Finalized(&PHIST_OUT_fini); \
                if (PHIST_OUT_ini && (!PHIST_OUT_fini)) { \
                MPI_Comm_rank(MPI_COMM_WORLD,&PHIST_OUT_me);\
                } else {PHIST_OUT_me=0;}\
                if(PHIST_OUT_me==0){\
                fprintf(PHIST_OUT_out,msg,##__VA_ARGS__);\
                fflush(PHIST_OUT_out);}\
        }\
}
#else
#define PHIST_SOUT(level,msg, ...) {\
        if(PHIST_OUTLEV >= level) {\
                FILE* PHIST_OUT_out= (level<=PHIST_WARNING)? stderr:stdout;\
                fprintf(PHIST_OUT_out,msg,##__VA_ARGS__);\
                fflush(PHIST_OUT_out);\
        }\
}
#endif


// "line-level" timings using PHIST_CHK macros
#if defined(__cplusplus) && defined(PHIST_TIMEMONITOR_PERLINE)
#include "phist_timemonitor.hpp"
#include <string>
#define PHIST_STRINGIFY_MACRO(l) #l
#define PHIST_FILE_LINE_REMOVE_PATH(f) (strrchr(f, '/') ? strrchr(f, '/') + 1 : f)
#define PHIST_FILE_LINE_MACRO(f,l) PHIST_FILE_LINE_REMOVE_PATH(f ":" PHIST_STRINGIFY_MACRO(l))
#define PHIST_TIMEMONITOR_PERLINE_MACRO PHIST_CXX_TIMER(PHIST_FILE_LINE_MACRO(__FILE__,__LINE__));
#else
#define PHIST_TIMEMONITOR_PERLINE_MACRO
#endif

//! checks an iflag flag passed to a void function for non-zero value, assigns it to FLAG,
//! prints an error message and returns if non-zero (to be used in void functions)
#ifdef __cplusplus
#define PHIST_CHK_IERR(func,FLAG) { PHIST_TIMEMONITOR_PERLINE_MACRO \
try {func; if (FLAG==PHIST_DEPRECATED) { \
PHIST_OUT(PHIST_WARNING,"Warning: function %s is DEPRECATED!\n (file %s, line %d)\n",(#func),(__FILE__),(__LINE__)); \
FLAG=PHIST_SUCCESS; \
} else if (FLAG!=PHIST_SUCCESS) { \
PHIST_OUT(PHIST_ERROR,"Error code %d (%s) returned from call %s\n(file %s, line %d)\n",\
(FLAG),(phist_retcode2str(FLAG)),(#func),(__FILE__),(__LINE__)); return;}} \
catch (const std::exception &e) {PHIST_OUT(PHIST_ERROR,"Exception caught in call %s (%s)\n(file %s, line %d)\n",\
(#func),e.what(),(__FILE__),(__LINE__)); (FLAG)=-77; return;} \
catch (const std::string &s) {PHIST_OUT(PHIST_ERROR,"Exception caught in call %s (%s)\n(file %s, line %d)\n",\
(#func),s.c_str(),(__FILE__),(__LINE__)); (FLAG)=-77; return;} \
catch (const int &iexc) {PHIST_OUT(PHIST_ERROR,"int Exception caught in call %s (value %d)\n(file %s, line %d)\n",\
(#func),iexc,(__FILE__),(__LINE__)); (FLAG)=-77; return;} \
catch (...) {PHIST_OUT(PHIST_ERROR,"unknown Exception caught in call %s\n(file %s, line %d)\n",\
(#func),(__FILE__),(__LINE__)); (FLAG)=-77; return;}}
#else
#define PHIST_CHK_IERR(func,FLAG) { PHIST_TIMEMONITOR_PERLINE_MACRO \
{func; if (FLAG==PHIST_DEPRECATED) { \
PHIST_OUT(PHIST_WARNING,"Warning: function %s is DEPRECATED!\n (file %s, line %d)\n",(#func),(__FILE__),(__LINE__)); \
FLAG=PHIST_SUCCESS; \
} else if (FLAG!=PHIST_SUCCESS) { \
PHIST_OUT(PHIST_ERROR,"Error code %d (%s) returned from call %s\n(file %s, line %d)\n",\
(FLAG),(phist_retcode2str(FLAG)),(#func),(__FILE__),(__LINE__)); return;}}}
#endif

#ifdef PHIST_HAVE_GHOST
#include "ghost/config.h"
#include "ghost/types.h"
// check return value from GHOST
#define PHIST_CHK_GERR(func,FLAG) { PHIST_TIMEMONITOR_PERLINE_MACRO \
ghost_error_t gerr=func; FLAG=(int)gerr; if (gerr!=GHOST_SUCCESS) { \
PHIST_OUT(PHIST_ERROR,"Error code %d (%s) returned from call %s\n(file %s, line %d)\n",\
(FLAG),(phist_ghost_error2str(gerr)),(#func),(__FILE__),(__LINE__)); return;}\
}
#endif
//! checks an iflag flag passed to a void function for negative value, assigns it to FLAG,
//! prints an error message and returns if non-zero (to be used in void functions)
#ifndef PHIST_CHK_NEG_IERR
#ifdef __cplusplus
#define PHIST_CHK_NEG_IERR(func,FLAG) { PHIST_TIMEMONITOR_PERLINE_MACRO \
try {func; if (FLAG < PHIST_SUCCESS) { \
PHIST_OUT(PHIST_ERROR,"Error code %d (%s) returned from call %s\n(file %s, line %d)\n",\
(FLAG),(phist_retcode2str(FLAG)),(#func),(__FILE__),(__LINE__)); return;}} \
catch (const std::exception &e) {PHIST_OUT(PHIST_ERROR,"Exception caught in call %s (%s)\n(file %s, line %d)\n",\
(#func),e.what(),(__FILE__),(__LINE__)); (FLAG)=-77; return;} \
catch (const std::string &s) {PHIST_OUT(PHIST_ERROR,"Exception caught in call %s (%s)\n(file %s, line %d)\n",\
(#func),s.c_str(),(__FILE__),(__LINE__)); (FLAG)=-77; return;} \
catch (const int &iexc) {PHIST_OUT(PHIST_ERROR,"int Exception caught in call %s (value %d)\n(file %s, line %d)\n",\
(#func),iexc,(__FILE__),(__LINE__)); (FLAG)=-77; return;} \
catch (...) {PHIST_OUT(PHIST_ERROR,"unknown Exception caught in call %s\n(file %s, line %d)\n",\
(#func),(__FILE__),(__LINE__)); (FLAG)=-77; return;}}
#else
#define PHIST_CHK_NEG_IERR(func,FLAG) {\
{func; if (FLAG<PHIST_SUCCESS) { \
PHIST_OUT(PHIST_ERROR,"Error code %d (%s) returned from call %s\n(file %s, line %d)\n",\
(FLAG),(phist_retcode2str(FLAG)),(#func),(__FILE__),(__LINE__)); return;}}}
#endif
#endif


//! like PHIST_CHK_IERR, but returns iflag (to be used in int functions returning an error code)
#define PHIST_ICHK_IERR(func,FLAG) { PHIST_TIMEMONITOR_PERLINE_MACRO \
{func; if ((FLAG)!=PHIST_SUCCESS && (FLAG)!=PHIST_DEPRECATED) { \
PHIST_OUT(PHIST_ERROR,"Error code %d (%s) returned from call %s\n(file %s, line %d)\n",\
(FLAG),(phist_retcode2str(FLAG)),(#func),(__FILE__),(__LINE__)); return FLAG;} \
else if (FLAG==PHIST_DEPRECATED) { FLAG=PHIST_SUCCESS; \
PHIST_OUT(PHIST_WARNING,"Warning, function %s is DEPRECATED.\n(file %s, line %d)\n",(#func),(__FILE__),(__LINE__));}}}

#define PHIST_ICHK_NEG_IERR(func,FLAG) { PHIST_TIMEMONITOR_PERLINE_MACRO \
{func; if (FLAG<PHIST_SUCCESS) { \
PHIST_OUT(PHIST_ERROR,"Error code %d (%s) returned from call %s\n(file %s, line %d)\n",\
(FLAG),(phist_retcode2str(FLAG)),(#func),(__FILE__),(__LINE__)); return FLAG;}}}

#ifdef __cplusplus
//! like PHIST_CHK_IERR, but throw an exception on non-zero iflag
#define PHIST_TCHK_IERR(func,FLAG) { PHIST_TIMEMONITOR_PERLINE_MACRO \
{func; if (FLAG!=PHIST_SUCCESS && FLAG!=PHIST_DEPRECATED) { \
PHIST_OUT(PHIST_ERROR,"Error code %d (%s) returned from call %s\n(file %s, line %d)\n",\
(FLAG),(phist_retcode2str(FLAG)),(#func),(__FILE__),(__LINE__)); throw FLAG;} \
else if (FLAG==PHIST_DEPRECATED) { FLAG=PHIST_SUCCESS; \
PHIST_OUT(PHIST_WARNING,"Warning, function %s is DEPRECATED.\n(file %s, line %d)\n",(#func),__FILE__,__LINE__);}}}

#endif

#if PHIST_OUTLEV>=PHIST_DEBUG
#define PHIST_DEB(msg, ...) PHIST_OUT(PHIST_DEBUG,msg,##__VA_ARGS__);
#else
#define PHIST_DEB(msg, ...)
#endif

/* PHIST_ENTER_FCN definition */
#ifdef __cplusplus
# ifdef SCOREP_USER_ENABLE
#   include <scorep/SCOREP_User.h>
# else
#   define SCOREP_USER_REGION(a, b)
# endif
# ifdef PHIST_TIMEMONITOR
#   include "phist_timemonitor.hpp"
# endif
# ifdef PHIST_KERNEL_LIB_GHOST
#   include "ghost/task.h"
#   define PHIST_GHOST_CHK_IN_TASK(s, iflag) { \
      ghost_task_t *curtask = NULL; \
      PHIST_CHK_GERR(ghost_task_cur(&curtask), iflag); \
      if( curtask == NULL ) { \
        PHIST_SOUT(PHIST_DEBUG, "Called %s outside a ghost task!\n", s); }\
      }
# else
#   define PHIST_GHOST_CHK_IN_TASK(s, iflag)
# endif
# include "phist_fcntrace.hpp"
# if defined(PHIST_TIMEMONITOR)
#     define PHIST_ENTER_FCN(s) phist_FcnTrace YouCantHaveMultiple_PHIST_ENTER_FCN_StatementsInOneScope(s);\
                                PHIST_CXX_TIMER(YouCantHaveMultiple_PHIST_ENTER_FCN_StatementsInOneScope.fcn().c_str());\
                                SCOREP_USER_REGION(s, SCOREP_USER_REGION_TYPE_FUNCTION);
# elif (PHIST_OUTLEV>=PHIST_TRACE) || defined(LIKWID_PERFMON)
#     define PHIST_ENTER_FCN(s) phist_FcnTrace YouCantHaveMultiple_PHIST_ENTER_FCN_StatementsInOneScope(s);\
                                SCOREP_USER_REGION(s, SCOREP_USER_REGION_TYPE_FUNCTION);
# else
#     define PHIST_ENTER_FCN(s) SCOREP_USER_REGION(s, SCOREP_USER_REGION_TYPE_FUNCTION);
# endif
#else
# define PHIST_ENTER_FCN(s)
#endif

/* PHIST_ENTER_KERNEL_FCN definition (used to prevent measuring nested kernel calls) */
#if defined(__cplusplus) || !defined(PHIST_TIMINGS_FULL_TRACE)
# define PHIST_ENTER_KERNEL_FCN(s) phist_CheckKernelFcnNesting s_(s); PHIST_ENTER_FCN(s_.str());
#else
# define PHIST_ENTER_KERNEL_FCN(s)
#endif

/* print a warning that an untested / experimental function is called */
#define PHIST_MARK_AS_EXPERIMENTAL(s) PHIST_SOUT(PHIST_WARNING, "Called experimental (untested) %s\n", s);

/* this macro can be used to avoid compiler warnings about unused variables */
#ifndef PHIST_TOUCH
#define PHIST_TOUCH(x) (void)(x);
#endif

#ifndef PHIST_CAST_PTR_FROM_VOID
#define PHIST_CAST_PTR_FROM_VOID(_TYPE_,_PTR_,_VPTR_,_FLAG_) \
_TYPE_ *_PTR_ = NULL; \
_PTR_=(_TYPE_*)(_VPTR_); \
if (_PTR_==NULL) {_FLAG_=-88; PHIST_OUT(PHIST_ERROR,"bad cast in file %s, line %d.\n",\
__FILE__,__LINE__); return;}
#endif

#ifdef PHIST_TIMEMONITOR
// this macro is used in the kernel lib adaptors so that the 
// timemonitor output contains information on the total number
// fo matrix-vector products (#spMVMs) instead of just #spMMVMs
#define PHIST_COUNT_MATVECS(vec) \
{\
  int nvec;\
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(vx, &nvec, iflag), *iflag);\
  for(int i = 0; i < nvec; i++)\
  {\
    PHIST_ENTER_FCN("phist_totalMatVecCount");\
  }\
}
#else
#define PHIST_COUNT_MATVECS(vec)
#endif

#endif
