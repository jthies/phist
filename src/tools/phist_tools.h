#ifndef PHIST_TOOLS_H
#define PHIST_TOOLS_H

#include "phist_config.h"

#ifndef DOXYGEN

/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#endif //DOXYGEN

#ifdef __cplusplus
extern "C" {
#endif
const char* phist_retcode2str(int code);
#ifdef PHIST_HAVE_GHOST
#include "ghost/config.h"
#include "ghost/types.h"
const char* phist_ghost_error2str(ghost_error_t code);
#endif
const char* phist_kernel_lib();
#ifdef __cplusplus
}
#endif
#endif
