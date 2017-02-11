#ifndef PHIST_MPI_KERNELS_H
#define PHIST_MPI_KERNELS_H

#include "phist_config.h"

#ifndef DOXYGEN

/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#else
typedef int MPI_Comm;
typedef int MPI_Request;
typedef int MPI_Status;
#endif

/* note: the phist_typedefs.h file is provided in the subdirectory
 where the interface is implemented (e.g. ghost/, tpetra/).
 */
#include "phist_typedefs.h"

#endif /* doxygen */

#ifdef __cplusplus
extern "C" {
#endif

/*! \defgroup mpi_kernels additional operations implemented for all libs at once
@{
*/


#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_mpi_kernels_decl.h"

#include "phist_gen_c.h"
#include "phist_mpi_kernels_decl.h"
#endif

#include "phist_gen_d.h"
#include "phist_mpi_kernels_decl.h"

#include "phist_gen_z.h"
#include "phist_mpi_kernels_decl.h"

#include "phist_gen_clean.h"

/*@}*/

#ifdef __cplusplus
} /*extern "C"*/
#endif

#endif
