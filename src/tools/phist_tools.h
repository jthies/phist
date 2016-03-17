#ifndef PHIST_TOOLS_H
#define PHIST_TOOLS_H

#include "phist_config.h"

#ifndef DOXYGEN

/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#ifdef __cplusplus
# include <cstdio>
#else
# include <stdio.h>
#endif
#endif //DOXYGEN

#ifdef __cplusplus
extern "C" {
#endif
const char* phist_retcode2str(int code);
#ifdef PHIST_HAVE_GHOST
#include "ghost/config.h"
#include "ghost/types.h"
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

#endif
