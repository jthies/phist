#ifndef PHIST_KERNELS_H
#define PHIST_KERNELS_H

#include "phist_config.h"

#ifndef NO_INCLUDES_IN_HEADERS
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "phist_macros.h"

// note: the phist_typedefs.h file is provided in the subdirectory
// where the interface is implemented (e.g. ghost/, tpetra/).
#include "phist_typedefs.h"
#endif

//! @file phist_kernels.h
//! @brief Definition of an abstract interface for basic operations needed by iterative solvers.
//!
//! The kernels are based on a number of abstract concepts, which are
//! simply void* into a "kernel library":
//! - comm_t: encapsulates MPI, currently the interface does not really
//!   allow you to do much with this object, except use it to build other
//!   objects.
//! - map_t: encapsulates the distribution of points (e.g. rows of a matrix)
//!   over processors
//! - multi-vector (mvec_t): a dense matrix with N rows and m<<N columns, with a row-wise
//!   distribution over computational nodes. The storage may be either column-major or
//!   or column-major. In the latter case, the macro PHIST_MVECS_ROW_MAJOR should be 
//!   defined.
//! - sdMat_t - a small serial dense matrix, replicated on all MPI processes.
//! - crsMat_t: sparse matrix in *any* storage format (TODO: rename it)
//!
//! The type-specific objects for mvecs, sdMats and crsMats are called
//! Smvec_t, Dmvec_t, Cmvec_t etc. and defined in phist_kernels_def.h

#ifdef __cplusplus
extern "C" {
#endif

//! \name data-type independent kernel functions
//!@{

//! initialize kernel library. Should at least call MPI_Init if it has not been called
//! but is required.
void phist_kernels_init(int *argc, char*** argv, int* ierr);

//! finalize kernel library. Should at least call MPI_Finalize if it has not been called
//! but is required.
void phist_kernels_finalize(int* ierr);

//! creates a global comm object
void phist_comm_create(comm_ptr_t* comm, int* ierr);
//! delete a comm object. Only do this for comms obtained by phist_comm_create.
void phist_comm_delete(comm_ptr_t comm, int* ierr);
//! get the rank of the calling node
void phist_comm_get_rank(const_comm_ptr_t comm, int* rank, int* ierr);
//! get the number of MPI asks
void phist_comm_get_size(const_comm_ptr_t comm, int* size, int* ierr);
//! creates a map with default distribution of points
void phist_map_create(map_ptr_t* map, const_comm_ptr_t comm, gidx_t nglob, int *ierr);
//! delete a map object. Note that you should not do this if you got the map from
//! anything else than phist_map_create.
void phist_map_delete(map_ptr_t map, int *ierr);
//! returns the comm object used by a map
void phist_map_get_comm(const_map_ptr_t map, const_comm_ptr_t* comm, int* ierr);
//! returns the local number of elements in the map
void phist_map_get_local_length(const_map_ptr_t map, lidx_t* nloc, int* ierr);
//! returns the smallest global index in the map appearing on my partition. ierr is set to 1
//! in case the map is not contiguous, because in that case it may be that the
//! caller falsely assumes global elements [ilower ... iupper] are actually on this partition.
void phist_map_get_ilower(const_map_ptr_t map, gidx_t* ilower, int* ierr);
//! returns the largest global index in the map appearing on my partition. ierr is set to 1
//! in case the map is not contiguous, because in that case it may be that the
//! caller falsely assumes global elements [ilower ... iupper] are actually on this partition.
void phist_map_get_iupper(const_map_ptr_t map, gidx_t* iupper, int* ierr);

#ifdef __cplusplus
} // extern "C"
#endif

//!@}

//! include file for the basic operations (kernels)
//! in single/double, real/complex.
#ifdef PHIST_HAVE_SP
/*!\name single precision real kernel functions */
/*!@{*/
#include "phist_gen_s.h"
#include "phist_kernels_decl.h"
#include "phist_carp_decl.h"
/*!@}*/
/*!\name single precision complex kernel functions */
/*!@{*/
#include "phist_gen_c.h"
#include "phist_kernels_decl.h"
#include "phist_carp_decl.h"
/*!@}*/
#endif

/*!\name double precision real kernel functions */
/*!@{*/
#include "phist_gen_d.h"
#include "phist_kernels_decl.h"
#include "phist_carp_decl.h"
/*!@}*/
/*!\name double precision complex kernel functions */
/*!@{*/
#include "phist_gen_z.h"
#include "phist_kernels_decl.h"
#include "phist_carp_decl.h"
/*!@}*/
#endif
