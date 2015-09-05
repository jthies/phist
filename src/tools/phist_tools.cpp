#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "phist_tools.h"
#include "phist_macros.h"
#include "phist_enums.h"
#include "phist_fcntrace.hpp"


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
  else if (code<=-1 && code>-10) return "function-specific error";
  else if (code==PHIST_DEPRECATED) return "deprecated";
  else if (code>0) return "warning";
  return "unknown error";
}

#ifdef PHIST_HAVE_GHOST
extern "C" const char* phist_ghost_error2str(ghost_error_t code)
{
      return ghost_error_string(code);
}
#endif

extern "C" const char* eigSort2str(eigSort_t s)
{
  return s==LM?"LM":
         s==SM?"SM":
         s==LR?"LR":
         s==SR?"SR":
         s==NONE?"none":
         s==TARGET?"target":
                   "invalid";
}

extern "C" const char* linSolv2str(linSolv_t s)
{
  return   s==GMRES?"GMRES":
           s==MINRES?"MINRES":
           s==CARP_CG?"CARP_CG":
           s==DO_NOTHING?"do_nothing":
           s==USER_DEFINED?"user_defined":
                         "invalid";
}

extern "C" eigSort_t str2eigSort(const char* c_str)
{
  std::string str(c_str);
  eigSort_t s=INVALID_EIGSORT_T;
  if (str=="LM") s=LM;
  else if (str=="SM") s=SM;
  else if (str=="LR") s=LR;
  else if (str=="SR") s=SR;
  else if (str=="none"||str=="NONE") s=NONE;
  else if (str=="target"||str=="TARGET") s=TARGET;
  return s;
}
      

extern "C" linSolv_t str2linSolv(const char* c_str)
{
  std::string str(c_str);
  linSolv_t s=INVALID_LINSOLV_T;
  if (str=="gmres"||str=="GMRES") s=GMRES;
  else if (str=="minres"||str=="MINRES") s=MINRES;
  else if (str=="carp_cg"||str=="CARP_CG") s=CARP_CG;
  else if (str=="do_nothing"||str=="DO_NOTHING") s=DO_NOTHING;
  else if (str=="user_defined"||str=="USER_DEFINED") s=USER_DEFINED;
  return s;
}

std::istream& operator>>(std::istream& is, eigSort_t& s)
{
  std::string tmp;
  is>>tmp;
  s=str2eigSort(tmp.c_str());
  return is;
}

std::istream& operator>>(std::istream& is, linSolv_t& s)
{
  std::string tmp;
  is>>tmp;
  s=str2linSolv(tmp.c_str());
  return is;
}

#ifdef PHIST_TIMINGS_FULL_TRACE
std::vector<const char*> phist_FcnTrace::fcnTrace_;
#endif

extern "C" const char* phist_kernel_lib()
{
#ifdef PHIST_KERNEL_LIB_BUILTIN
  return "builtin";
#elif defined(PHIST_KERNEL_LIB_GHOST)
  return "ghost";
#elif defined(PHIST_KERNEL_LIB_EPETRA)
  return "epetra";
#elif defined(PHIST_KERNEL_LIB_TPETRA)
  return "tpetra";
#else
#error "No appropriate PHIST_KERNEL_LIB_... defined!"
#endif
}

#ifdef PHIST_HAVE_CXX11_THREADLOCAL
thread_local bool phist_CheckKernelFcnNesting::nestedKernelCall_ = false;
#else
bool phist_CheckKernelFcnNesting::nestedKernelCall_ = false;
#endif
