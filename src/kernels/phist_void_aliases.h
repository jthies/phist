/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
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

//! a context is an abstract object describing the shape of a matrix and
//! parallel distribution of a matrix. It contains all the information  
//! needed to create a new matrix operating on the same type of vector  
//! objects. In the Petra object model it contains pointers to the four 
//! maps (row, col, range and domain) of a sparse matrix.
typedef void* phist_context_ptr;

//! pointer to const context
typedef const void* phist_const_context_ptr;

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
