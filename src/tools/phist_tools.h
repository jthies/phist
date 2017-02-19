#ifndef PHIST_TOOLS_H
#define PHIST_TOOLS_H

#include "phist_config.h"

#ifndef DOXYGEN

/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#ifdef __cplusplus
#include "phist_ScalarTraits.hpp"
#endif

#ifdef __cplusplus
# include <cstdio>
#else
# include <stdio.h>
#endif

#ifdef PHIST_HAVE_LIKWID
#include <likwid.h>
#endif

#endif //DOXYGEN

# include "phist_macros.h"

#ifdef PHIST_HAVE_GHOST
#include <ghost.h>
#endif

#ifdef __cplusplus

// get the size of a cache line on the current MPI process,
// or 1 if the architecture doesn't have a cache. The unit is
// sizeof(T), so on a typical CPU with 64 bit cache lines
// phist_cacheline_size<double>()==8 and phist_cacheline_size<char>()==64.
template<typename T>
unsigned int phist_cacheline_size()
{
  static unsigned int clsize=-1;
  if (clsize>0) return clsize;
#ifdef PHIST_HAVE_GHOST
  ghost_machine_cacheline_size(&clsize);
# ifdef PHIST_KERNEL_LIB_GHOST
  ghost_type type;
  ghost_type_get(&type);
  if (type==GHOST_TYPE_CUDA) clsize=sizeof(T);
# endif  
  
#else
  PHIST_SOUT(PHIST_WARNING,"guessing cacheline length (64bit) because ghost is not available.\n");
  clsize=64;
#endif
  clsize/=sizeof(T);
  return clsize;
}

extern "C" {
#endif
const char* phist_retcode2str(int code);
#ifdef PHIST_HAVE_GHOST
#include <ghost/config.h>
#include <ghost/types.h>
const char* phist_ghost_error2str(ghost_error code);
#endif
const char* phist_kernel_lib();

# ifdef PHIST_HAVE_MPI

//! pretty-print process-local strings in order. This function should
//! not be used directly but via the wrapper macro PHIST_ORDERED_OUT(...)
int phist_ordered_fprintf(FILE* stream, MPI_Comm comm, const char* fmt, ...);

# endif

#ifdef __cplusplus
}
#endif

#include "phist_tasks.h"

#ifndef DOXYGEN
#include "phist_lapack.h"
#endif

#endif
