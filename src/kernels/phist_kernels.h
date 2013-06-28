#ifndef PHIST_KERNELS_H
#define PHIST_KERNELS_H

//! note: the phist_typedefs.h file is provided in the subdirectory
//! where the interface is implemented (e.g. ghost/, tpetra/).
#include "phist_typedefs.h"

//! the kernels are based on a number of abstract concepts, which are
//! simply void* into a "kernel library":
//! - comm_t: encapsulates MPI, currently the interface does not really
//!   allow you to do much with this object, except use it to build other
//!   objects.
//! - map_t: encapsulates the distribution of points (e.g. rows of a matrix)
//!   over processors
//! - multi-vector (mvec_t): a dense matrix with N rows and m<<N columns, with a row-wise
//!   distribution over computational nodes.
//! - crsMat_t: sparse matrix in CRS format.
//!
//! The type-specific objects for mvecs and crsMats are called
//! Smvec_t, Dmvec_t, Cmvec_t etc. and defined in phist_kernels_def.h

#ifdef __cplusplus
extern "C" {
#endif

//! initialize kernel library. Should at least call MPI_Init if it has not been called
//! but is required.
void phist_kernels_init(int *argc, char*** argv, int* ierr);

//! finalize kernel library. Should at least call MPI_Finalize if it has not been called
//! but is required.
void phist_kernels_finalize(int* ierr);

//! creates a global comm object
void phist_comm_create(comm_ptr_t* comm, int* ierr);
//! get the rank of the calling node
void phist_comm_get_rank(const_comm_ptr_t comm, int* rank, int* ierr);
//! get the number of MPI asks
void phist_comm_get_size(const_comm_ptr_t comm, int* size, int* ierr);
//! creates a map with default distribution of points
void phist_map_create(map_ptr_t* map, const_comm_ptr_t comm, int nglob, int *ierr);
//! returns the comm object used by a map
void phist_map_get_comm(const_map_ptr_t map, const_comm_ptr_t* comm, int* ierr);

#ifdef __cplusplus
} // extern "C"
#endif

//! include file for the basic operations (kernels)
//! in single/double, real/complex.
#include "phist_gen_s.h"
#include "phist_kernels_decl.h"
#include "phist_gen_d.h"
#include "phist_kernels_decl.h"
#include "phist_gen_c.h"
#include "phist_kernels_decl.h"
#include "phist_gen_z.h"
#include "phist_kernels_decl.h"

#endif
