#include "phist_macros.h"
#include "phist_trilinos_macros.h"
#include "../phist_kernels.h"

#include "phist_typedefs.h"
#include "phist_ScalarTraits.hpp"
#include "BelosTsqrOrthoManager.hpp"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "ghost.h"
#include "ghost_vec.h"
#include "ghost_util.h"

extern "C" {

#ifndef PHIST_HAVE_MPI
typedef int MPI_Comm;
const int MPI_COMM_WORLD=0;
#endif

//                                                                                        
// in the petra object model that this abstract kernel interface follows, a map describes 
// the distribution of data among processors. It is used to e.g. define how vector entries
// or matrix rows are distributed among MPI processes and required to create a vector.    
//                                                                                        
// In ghost this is handled differently, the 'context' is the main object for describing  
// communication patterns, and it is owned by a sparse matrix. In order to create a vector
// we need the vtraits_t object, which also knows about number of vectors and             
// data type. So we define our own struct with a pointer to the context object and the    
// flags required to create the traits object in mvec_create(). 
typedef struct ghost_map_t
  {
  ghost_context_t* ctx;
  int vtraits_flags;
  } ghost_map_t;

// initialize ghost
void phist_kernels_init(int* argc, char*** argv, int* ierr)
  {
  *ierr=0;
  ghost_init(*argc, *argv);
  }

// finalize ghost
void phist_kernels_finalize(int* ierr)
  {
  ghost_finish();
  }


//! simply returns MPI_COMM_WORLD, the only MPI_Comm used in ghost.
void phist_comm_create(comm_ptr_t* vcomm, int* ierr)
  {
  *ierr=0;
  // concept doesn't exist in ghost, return MPI_Comm
  MPI_Comm* comm = new MPI_Comm;
  *comm=MPI_COMM_WORLD;
  *vcomm=(comm_ptr_t)comm;
  }

//!
void phist_comm_delete(comm_ptr_t vcomm, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(MPI_Comm,comm,vcomm,*ierr);
  delete comm;
  }

//!
void phist_comm_get_rank(const_comm_ptr_t vcomm, int* rank, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(MPI_Comm,comm,vcomm,*ierr);
  *rank=ghost_getRank(*comm);
  }
//!
void phist_comm_get_size(const_comm_ptr_t vcomm, int* size, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(MPI_Comm,comm,vcomm,*ierr);
  *size=ghost_getNumberOfRanks(*comm);
  }

//! this generates a default map with linear distribution of points among
//! processes. Vectors based on this map will have halo points so that   
//! they can be used as either X or Y in Y=A*X operations.
void phist_map_create(map_ptr_t* vmap, const_comm_ptr_t vcomm, gidx_t nglob, int *ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(MPI_Comm,comm,vcomm,*ierr);

  ghost_map_t* map = new ghost_map_t;
  
  map->ctx=ghost_createContext(nglob, nglob, GHOST_CONTEXT_DEFAULT, NULL,*comm,1.0);
  map->vtraits_flags = GHOST_VEC_LHS;
  *vmap=(map_ptr_t)(map);
  }

//!
void phist_map_delete(map_ptr_t vmap, int *ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(ghost_map_t,map,vmap,*ierr);
  // this is problematic because the object may be shared, so we don't delete
  // anything right now (it's just a tiny amount of data to leak for now)
  // delete map->ctx;
  delete map;
  vmap=NULL;
  }
  
//!
void phist_map_get_comm(const_map_ptr_t vmap, const_comm_ptr_t* vcomm, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const ghost_map_t,map,vmap,*ierr);
  *vcomm = (const_comm_ptr_t)(&map->ctx->mpicomm);
  }

//!
void phist_map_get_local_length(const_map_ptr_t vmap, int* nloc, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const ghost_map_t,map,vmap,*ierr);
  *nloc=map->ctx->communicator->lnrows[ghost_getRank(MPI_COMM_WORLD)];
  }

//! returns the smallest global index in the map appearing on my partition.
void phist_map_get_ilower(const_map_ptr_t vmap, int* ilower, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const ghost_map_t,map,vmap,*ierr);
  int me = ghost_getRank(map->ctx->mpicomm);
  *ilower = map->ctx->communicator->lfRow[me];
  }
//! returns the largest global index in the map appearing on my partition.
void phist_map_get_iupper(const_map_ptr_t vmap, int* iupper, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const ghost_map_t,map,vmap,*ierr);
  int me = ghost_getRank(map->ctx->mpicomm);
  *iupper = map->ctx->communicator->lfRow[me+1]-1;
  }


} //extern "C"

#include "phist_gen_s.h"
#include "kernels_def.hpp"

#include "phist_gen_d.h"
#include "kernels_def.hpp"

#include "phist_gen_c.h"
#include "kernels_def.hpp"

#include "phist_gen_z.h"
#include "kernels_def.hpp"

