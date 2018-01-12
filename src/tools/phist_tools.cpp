/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
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

#include <cstring>
#include <cstdarg>

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

extern "C" const char* eigSort2str(phist_EeigSort s)
{
  return s==phist_LM?"LM":
         s==phist_SM?"SM":
         s==phist_LR?"LR":
         s==phist_SR?"SR":
         s==phist_NO_EIGSORT?"NONE":
         s==phist_TARGET?"TARGET":
                   "INVALID";
}

extern "C" const char* matSym2str(phist_EmatSym s)
{
  return s==phist_GENERAL?"GENERAL":
         s==phist_HERMITIAN?"HERMITIAN":
         s==phist_COMPLEX_SYMMETRIC?"COMPLEX_SYMMETRIC":
         s==phist_PATTERN_SYMMETRIC?"PATTERN_SYMMETRIC":
                   "INVALID";
}

extern "C" const char* eigExtr2str(phist_EeigExtr s)
{
  return s==phist_STANDARD?"STANDARD":
         s==phist_HARMONIC?"HARMONIC":
                   "INVALID";
}

extern "C" const char* linSolv2str(phist_ElinSolv s)
{
  return   s==phist_GMRES?"GMRES":
           s==phist_MINRES?"MINRES":
           s==phist_QMR?"QMR":
           s==phist_BICGSTAB?"BICGSTAB":
           s==phist_CARP_CG?"CARP_CG":
           s==phist_NO_LINSOLV?"NONE":
           s==phist_USER_LINSOLV?"USER_DEFINED":
                         "INVALID";
}

extern "C" const char* precon2str(phist_Eprecon s)
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
           s==phist_USER_PRECON?"USER_DEFINED":
                "INVALID";
}

extern "C" const char* projection2str(phist_Eprojection s)
{
  return   s==phist_PROJ_NONE?"NONE":
           s==phist_PROJ_PRE?"PRE":
           s==phist_PROJ_POST?"POST":
           s==phist_PROJ_PRE_POST?"PRE_POST":
           s==phist_PROJ_SKEW?"SKEW":
           s==phist_PROJ_ALL?"ALL":
                "INVALID";
}

extern "C" phist_EeigSort str2eigSort(const char* c_str)
{
  std::string str(c_str);
  str=phist_str2upper(str);
  phist_EeigSort s=phist_INVALID_EIGSORT;
  if (str=="LM") s=phist_LM;
  else if (str=="SM") s=phist_SM;
  else if (str=="LR") s=phist_LR;
  else if (str=="SR") s=phist_SR;
  else if (str=="NONE") s=phist_NO_EIGSORT;
  else if (str=="TARGET") s=phist_TARGET;
  return s;
}

extern "C" phist_EmatSym str2matSym(const char* c_str)
{
  std::string str(c_str);
  str=phist_str2upper(str);
  phist_EmatSym s=phist_INVALID_MATSYM;
  if (str=="GENERAL") s=phist_GENERAL;
  else if (str=="HERMITIAN") s=phist_HERMITIAN;
  else if (str=="COMPLEX_SYMMETRIC") s=phist_COMPLEX_SYMMETRIC;
  else if (str=="PATTERN_SYMMETRIC") s=phist_PATTERN_SYMMETRIC;
  else if (str=="NONE") s=phist_GENERAL;
  return s;
}

extern "C" phist_EeigExtr str2eigExtr(const char* c_str)
{
  std::string str(c_str);
  str=phist_str2upper(str);
  phist_EeigExtr s=phist_INVALID_EIGEXTR;
  if (str=="STANDARD") s=phist_STANDARD;
  else if (str=="HARMONIC") s=phist_HARMONIC;
  return s;
}
      

extern "C" phist_ElinSolv str2linSolv(const char* c_str)
{
  std::string str(c_str);
  str=phist_str2upper(str);
  phist_ElinSolv s=phist_INVALID_LINSOLV;
  if (str=="GMRES") s=phist_GMRES;
  else if (str=="MINRES") s=phist_MINRES;
  else if (str=="QMR") s=phist_QMR;
  else if (str=="BICGSTAB") s=phist_BICGSTAB;
  else if (str=="CARP_CG") s=phist_CARP_CG;
  else if (str=="NONE") s=phist_NO_LINSOLV;
  else if (str=="USER_DEFINED") s=phist_USER_LINSOLV;
  return s;
}

extern "C" phist_Eprecon str2precon(const char* c_str)
{
  std::string str(c_str);
  str=phist_str2upper(str);
  phist_Eprecon s=phist_INVALID_PRECON;
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
  else if (str=="USER_DEFINED") s=phist_USER_PRECON;
  return s;
}

extern "C" phist_Eprojection str2projection(const char* c_str)
{
  std::string str(c_str);
  str=phist_str2upper(str);
  phist_Eprojection s=phist_INVALID_PROJ;
  if (str=="NONE") s=phist_PROJ_NONE;
  else if (str=="PRE") s=phist_PROJ_PRE;
  else if (str=="POST") s=phist_PROJ_POST;
  else if (str=="PRE_POST") s=phist_PROJ_PRE_POST;
  else if (str=="SKEW") s=phist_PROJ_SKEW;
  else if (str=="ALL") s=phist_PROJ_ALL;
  return s;
}

#define INST_IO_OP(_EWHAT,_WHAT2STR,_STR2WHAT) \
std::istream& operator>>(std::istream& is, _EWHAT& e) \
{ \
  std::string tmp; \
  is>>tmp; \
  PHIST_SOUT(PHIST_DEBUG,"try to parse %s '%s'\n",#_EWHAT,tmp.c_str()); \
  e=_STR2WHAT(tmp.c_str()); \
  return is; \
} \
\
std::ostream& operator<< (std::ostream& os, const _EWHAT& e) \
{ \
  os << _WHAT2STR(e); \
  return os; \
}

INST_IO_OP(phist_EeigSort,eigSort2str,str2eigSort)
INST_IO_OP(phist_EeigExtr,eigExtr2str,str2eigExtr)
INST_IO_OP(phist_ElinSolv,linSolv2str,str2linSolv)
INST_IO_OP(phist_EmatSym,matSym2str,str2matSym)
INST_IO_OP(phist_Eprecon,precon2str,str2precon)
INST_IO_OP(phist_Eprojection,projection2str,str2projection)


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
#elif defined(PHIST_KERNEL_LIB_PETSC)
  return "petsc";
#elif defined(PHIST_KERNEL_LIB_EIGEN)
  return "Eigen";
#else
#error "No appropriate PHIST_KERNEL_LIB_... defined!"
#endif
}

#ifdef PHIST_HAVE_CXX11_THREADLOCAL
thread_local bool phist_CheckKernelFcnNesting::nestedKernelCall_ = false;
#else
bool phist_CheckKernelFcnNesting::nestedKernelCall_ = false;
#endif

#ifdef PHIST_HAVE_MPI
//! pretty-print process-local strings in order. This function should
//! not be used directly but via the wrapper macro PHIST_ORDERED_OUT(...)
//! the return value is the number of characters contributed by this process.
extern "C" int phist_ordered_fprintf(FILE* stream, MPI_Comm comm, const char* fmt, ...)
{
  int rank, size;
  char *local_string=NULL;
  int local_length, global_length;
  int *char_counts=NULL, *char_disps=NULL;
  char* global_string=NULL;
  va_list args;
  va_start(args, fmt);
  local_length=vasprintf(&local_string,fmt,args);
  
  // use MPI to gather the global string
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);
  if (rank==0)
  {
    char_counts=new int[size];
    char_disps=new int[size+1];
  }
  // gather remote string lengths
  MPI_Gather(&local_length,1,MPI_INT,
              char_counts, 1,MPI_INT,
              0, comm);
  
  if (rank==0)
  {
    char_disps[0]=0;
    for (int i=0; i<size;i++) char_disps[i+1]=char_disps[i]+char_counts[i];
    global_length=char_disps[size];
    if (global_length>1e8)
    {
      fprintf(stderr,"WARNING: you're gathering a very large string, PHIST_ORDERED_OUT is intended for short messages\n");
    }
    global_string=new char[global_length+1];
  }
  MPI_Gatherv(local_string,local_length,MPI_CHAR,
              global_string,char_counts,char_disps,MPI_CHAR,
              0, comm);
  // print the global string
  if (rank==0)
  {
    global_string[global_length]='\0';
    fprintf(stream,global_string);
  }
  
  // clean up the mess
  // we use free() here because we got the string from vasprintf and the GNU address
  // sanitizer gets angry if we delete [] it.
  free(local_string);
  if (rank==0)
  {
    delete [] char_counts;
    delete [] char_disps;
    delete [] global_string;
  }
  // return the number of characters contributed by this process
  return local_length;
}

#endif
