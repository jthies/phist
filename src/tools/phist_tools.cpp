#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "phist_tools.h"
#include "phist_macros.h"
#include "phist_enums.h"


extern "C" const char* phist_retcode2str(int code)
  {
  if (code==PHIST_SUCCESS) return "success";
  else if (code==PHIST_FUNCTIONAL_ERROR) return "functional error";
  else if (code==PHIST_MEM_ALLOC_FAILED) return "memory allocation failed";
  else if (code==PHIST_INVALID_INPUT) return "invalid input";
  else if (code==PHIST_INTEGER_OVERFLOW) return "int overflow";
  else if (code==PHIST_CAUGHT_EXCEPTION) return "caught exception";
  else if (code==PHIST_BAD_CAST) return "bad cast or NULL pointer";
  else if (code==PHIST_NOT_IMPLEMENTED) return "not implemented";
  return "unknown error";
  }

#ifdef PHIST_KERNEL_LIB_GHOST
extern "C" const char* phist_ghost_error2str(ghost_error_t code)
  {
      return ghost_error_string(code);
  }
#endif

extern "C" const char* eigSort2str(eigSort_t s)
  {
  return s==LM?"LM":s==SM?"SM":s==LR?"LR":s==SR?"SR":s==NONE?"none":s==TARGET?"target":"unknown eigsort_t";
  }
      

extern "C" const char* linSolv2str(linSolv_t s)
  {
  return s==GMRES?"GMRES":s==CARP_CG?"CARP-CG":s==DO_NOTHING?"do nothing":"unknown linSolv_t";
  }

#ifdef PHIST_TIMEMONITOR
# ifndef PHIST_HAVE_TEUCHOS
# include "phist_TimeMonitor.hpp"
namespace phist_TimeMonitor
{
  Timer::TimeDataMap Timer::_timingResults;
}
# endif  
#endif  
