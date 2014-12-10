#ifndef PHIST_VOID_ALIASES_H
#define PHIST_VOID_ALIASES_H

//! wrapper for the MPI communicator
typedef void* comm_ptr_t;
//! pointer to const MPI communicator
typedef const void* const_comm_ptr_t;

//! a map is an object defining the distribution of points over compute nodes
typedef void* map_ptr_t;

//! pointer to const map
typedef const void* const_map_ptr_t;

#include "phist_gen_s.h"
#include "phist_void_aliases_decl.h"
#include "phist_gen_d.h"
#include "phist_void_aliases_decl.h"
#include "phist_gen_c.h"
#include "phist_void_aliases_decl.h"
#include "phist_gen_z.h"
#include "phist_void_aliases_decl.h"
#include "phist_gen_clean.h"

#endif
