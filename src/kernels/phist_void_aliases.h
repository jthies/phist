#ifndef PHIST_VOID_ALIASES_H
#define PHIST_VOID_ALIASES_H

//! wrapper for the MPI communicator
typedef void* phist_comm_ptr;
//! pointer to const MPI communicator
typedef const void* phist_const_comm_ptr;

//! a map is an object defining the distribution of points over compute nodes
typedef void* phist_map_ptr;

//! pointer to const map
typedef const void* phist_const_map_ptr;

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
