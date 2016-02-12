#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "phist_tools.h"
#include "phist_macros.h"
#include "phist_enums.h"
#include "phist_fcntrace.hpp"

#include <algorithm>
#include <string>

// little helper utiliity so that we can recognize strings regardless of case,
// e.g. carp_cg, CARP_CG, carp_CG => CARP_CG. Note that "carp-cg" won't work
// because the dash is not transformed, though.
std::string phist_str2upper(const std::string& s)
{
  std::string S = s;
  std::transform(S.begin(), S.end(), S.begin(), ::toupper);
  return S;
}

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
extern "C" const char* phist_ghost_error2str(ghost_error code)
{
      return ghost_error_string(code);
}
#endif

extern "C" const char* eigSort2str(eigSort_t s)
{
  return s==phist_LM?"LM":
         s==phist_SM?"SM":
         s==phist_LR?"LR":
         s==phist_SR?"SR":
         s==phist_NO_EIGSORT?"NONE":
         s==phist_TARGET?"TARGET":
                   "INVALID";
}

extern "C" const char* eigExtr2str(eigExtr_t s)
{
  return s==phist_STANDARD?"STANDARD":
         s==phist_HARMONIC?"HARMONIC":
                   "INVALID";
}

extern "C" const char* linSolv2str(linSolv_t s)
{
  return   s==phist_GMRES?"GMRES":
           s==phist_MINRES?"MINRES":
           s==phist_CARP_CG?"CARP_CG":
           s==phist_NO_LINSOLV?"NONE":
           s==phist_USER_DEFINED?"USER_DEFINED":
                         "INVALID";
}

extern "C" const char* precon2str(precon_t s)
{
  return   s==phist_NO_PRECON?"NONE":
#ifdef PHIST_HAVE_IFPACK
           s==phist_IFPACK?"IFPACK":
#elif defined(PHIST_HAVE_IFPACK2)
           s==phist_IFPACK?"IFPACK2":
#endif
#ifdef PHIST_HAVE_ML
           s==phist_ML?"ML":
#endif
#ifdef PHIST_HAVE_MUELU
           s==phist_MUELU?"MUELU":
#endif
#ifdef PHIST_HAVE_AMESOS2
           s==phist_AMESOS2?"AMESOS2":
#endif
                "INVALID";
}

extern "C" eigSort_t str2eigSort(const char* c_str)
{
  std::string str(c_str);
  str=phist_str2upper(str);
  eigSort_t s=phist_INVALID_EIGSORT_T;
  if (str=="LM") s=phist_LM;
  else if (str=="SM") s=phist_SM;
  else if (str=="LR") s=phist_LR;
  else if (str=="SR") s=phist_SR;
  else if (str=="NONE") s=phist_NO_EIGSORT;
  else if (str=="TARGET") s=phist_TARGET;
  return s;
}

extern "C" eigExtr_t str2eigExtr(const char* c_str)
{
  std::string str(c_str);
  str=phist_str2upper(str);
  eigExtr_t s=phist_INVALID_EIGEXTR_T;
  if (str=="STANDARD") s=phist_STANDARD;
  else if (str=="HARMONIC") s=phist_HARMONIC;
  return s;
}
      

extern "C" linSolv_t str2linSolv(const char* c_str)
{
  std::string str(c_str);
  str=phist_str2upper(str);
  linSolv_t s=phist_INVALID_LINSOLV_T;
  if (str=="GMRES") s=phist_GMRES;
  else if (str=="MINRES") s=phist_MINRES;
  else if (str=="CARP_CG") s=phist_CARP_CG;
  else if (str=="NONE") s=phist_NO_LINSOLV;
  else if (str=="user_defined"||str=="USER_DEFINED") s=phist_USER_DEFINED;
  return s;
}

extern "C" precon_t str2precon(const char* c_str)
{
  std::string str(c_str);
  str=phist_str2upper(str);
  precon_t s=phist_INVALID_PRECON_T;
  if (str=="NONE") s=phist_NO_PRECON;
#ifdef PHIST_HAVE_IFPACK
  else if (str=="IFPACK") s=phist_IFPACK;
#endif
#ifdef PHIST_HAVE_IFPACK2
  else if (str=="IFPACK"||str=="IFPACK2") s=phist_IFPACK;
#endif
#ifdef PHIST_HAVE_ML
  else if (str=="ML") s=phist_ML;
#endif
#ifdef PHIST_HAVE_MUELU
  else if (str=="MUELU") s=phist_MUELU;
#endif
#ifdef PHIST_HAVE_AMESOS2
  else if (str=="AMESOS2") s=phist_AMESOS2;
#endif
  return s;
}

std::istream& operator>>(std::istream& is, eigSort_t& s)
{
  std::string tmp;
  is>>tmp;
  PHIST_SOUT(PHIST_DEBUG,"try to parse eigSort_t '%s'\n",tmp.c_str());
  s=str2eigSort(tmp.c_str());
  return is;
}

std::istream& operator>>(std::istream& is, eigExtr_t& s)
{
  std::string tmp;
  is>>tmp;
  PHIST_SOUT(PHIST_DEBUG,"try to parse eigExtr_t '%s'\n",tmp.c_str());
  s=str2eigExtr(tmp.c_str());
  return is;
}

std::istream& operator>>(std::istream& is, linSolv_t& s)
{
  std::string tmp;
  is>>tmp;
  PHIST_SOUT(PHIST_DEBUG,"try to parse linSolv_t '%s'\n",tmp.c_str());
  s=str2linSolv(tmp.c_str());
  return is;
}

std::istream& operator>>(std::istream& is, precon_t& s)
{
  std::string tmp;
  is>>tmp;
  PHIST_SOUT(PHIST_DEBUG,"try to parse precon_t '%s'\n",tmp.c_str());
  s=str2precon(tmp.c_str());
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
#elif defined(PHIST_KERNEL_LIB_MAGMA)
  return "magma";
#else
#error "No appropriate PHIST_KERNEL_LIB_... defined!"
#endif
}

#ifdef PHIST_HAVE_CXX11_THREADLOCAL
thread_local bool phist_CheckKernelFcnNesting::nestedKernelCall_ = false;
#else
bool phist_CheckKernelFcnNesting::nestedKernelCall_ = false;
#endif
