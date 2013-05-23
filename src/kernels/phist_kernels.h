#ifndef KERNELS_H
#define KERNELS_H

typedef void* comm_ptr_t;
typedef const void* const_comm_ptr_t;

typedef void* map_ptr_t;
typedef const void* const_map_ptr_t;


#ifdef __cplusplus
extern "C" {
#endif

//! initialize kernel library. Should at least call MPI_Init if it has not been called
//! but is required.
void phist_kernels_init(int *argc, char*** argv, int* ierr);

//! finalize kernel library. Should at least call MPI_Finalize if it has not been called
//! but is required.
void phist_kernels_finalize(int* ierr);

//!
void phist_comm_create(comm_ptr_t* comm, int* ierr);
//!
void phist_comm_get_rank(const_comm_ptr_t comm, int* rank, int* ierr);
//!
void phist_comm_get_size(const_comm_ptr_t comm, int* size, int* ierr);
//!
void phist_map_create(map_ptr_t* map, const_comm_ptr_t comm, int nglob, int *ierr);
//!
void phist_map_get_comm(const_map_ptr_t map, const_comm_ptr_t* comm, int* ierr);

#ifdef __cplusplus
} // extern "C"
#endif

//! include file for the basic operations (kernels)
//! in single/double, real/complex.
#include "gen_s.h"
#include "phist_kernels_decl.h"
#include "gen_d.h"
#include "phist_kernels_decl.h"
#include "gen_c.h"
#include "phist_kernels_decl.h"
#include "gen_z.h"
#include "phist_kernels_decl.h"

#endif
