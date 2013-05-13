#ifndef KERNELS_H
#define KERNELS_H

typedef void* comm_ptr_t;
typedef const void* const_comm_ptr_t;

typedef void* map_ptr_t;
typedef const void* const_map_ptr_t;

//!
void comm_create(comm_ptr_t* comm, int* ierr);
//!
void comm_get_rank(const_comm_ptr_t* comm, int* rank, int* ierr);
//!
void comm_get_size(const_comm_ptr_t* comm, int* size, int* ierr);
//!
void map_create(map_ptr_t* map, const_comm_ptr_t comm, int nglob, int *ierr);
//!
void map_get_comm(const_map_ptr_t map, const_comm_ptr_t* comm, int* ierr);


//! include file for the basic operations (kernels)
//! in single/double, real/complex.
#include "gen_s.h"
#include "kernels_decl.h"
#include "gen_d.h"
#include "kernels_decl.h"
#include "gen_c.h"
#include "kernels_decl.h"
#include "gen_z.h"
#include "kernels_decl.h"

#endif
