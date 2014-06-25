#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include <iostream>
#include "phist_macros.h"
#include "../phist_kernels.h"

#include "phist_typedefs.h"
#include "typedefs.hpp"
#include "phist_ScalarTraits.hpp"

// these are from Trilinos, we need them to interface
// the TSQR library for orthogonalizing tall skinny matrices.
#ifdef PHIST_HAVE_BELOS
#include "phist_trilinos_macros.h"
#include "Ghost_TsqrAdaptor.hpp"
#include "Belos_GhostAdapter.hpp"
#include "BelosTsqrOrthoManager.hpp"
#endif
#ifdef PHIST_TIMEMONITOR
#include <Teuchos_TimeMonitor.hpp>
#endif

#include "phist_GhostMV.hpp"

#ifdef PHIST_HAVE_LIKWID
#include <likwid.h>
#endif

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <ghost.h>
#include <ghost/machine.h>
#include <limits>

namespace phist
  {
  int GhostMV::countObjects=0;
  }

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
// data type. So we define our own struct with a pointer to the context object and a temp-
// late for cloning the vtraits in mvec_create
typedef struct ghost_map_t
  {
  ghost_context_t* ctx;
  ghost_densemat_traits_t vtraits_template;
  } ghost_map_t;

// initialize ghost
void phist_kernels_init(int* argc, char*** argv, int* ierr)
  {
  *ierr=0;
  ghost_init(*argc, *argv);
  //ghost_pinThreads(GHOST_PIN_PHYS,NULL);
  char *str;
  ghost_string(&str);
  printf("%s\n",str);
  free(str); str = NULL;
  ghost_machine_string(&str);
  printf("%s\n",str);
  free(str); str = NULL;
#ifdef PHIST_HAVE_LIKWID
  LIKWID_MARKER_INIT;
  LIKWID_MARKER_START("phist<ghost>");
#endif
  }

// finalize ghost
void phist_kernels_finalize(int* ierr)
  {
#ifdef PHIST_HAVE_LIKWID
  LIKWID_MARKER_STOP("phist<ghost>");
  LIKWID_MARKER_CLOSE;
#endif
#ifdef PHIST_TIMEMONITOR
  Teuchos::TimeMonitor::summarize();
#endif
  ghost_finalize();
  *ierr=0;
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
  CAST_PTR_FROM_VOID(MPI_Comm,comm,vcomm,*ierr);
  delete comm;
  }

//!
void phist_comm_get_rank(const_comm_ptr_t vcomm, int* rank, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(MPI_Comm,comm,vcomm,*ierr);
  ghost_rank(rank,*comm);
  }
//!
void phist_comm_get_size(const_comm_ptr_t vcomm, int* size, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(MPI_Comm,comm,vcomm,*ierr);
  ghost_nrank(size,*comm);
  }

//! private helper function to create a vtraits object
ghost_densemat_traits_t phist_default_vtraits()
  {
  ghost_densemat_traits_t vtraits = GHOST_DENSEMAT_TRAITS_INITIALIZER;
  vtraits.nrows=0; // get from context
  vtraits.nrowshalo=0; // get from context
  vtraits.nrowspadded=0; // get from context
  vtraits.flags = (ghost_densemat_flags_t)(GHOST_DENSEMAT_DEFAULT|GHOST_DENSEMAT_HOST|GHOST_DENSEMAT_NO_HALO);
  vtraits.ncols=1;
#ifdef PHIST_MVECS_ROW_MAJOR
  vtraits.storage=GHOST_DENSEMAT_ROWMAJOR;
#else  
  vtraits.storage=GHOST_DENSEMAT_COLMAJOR;
#endif
  return vtraits;
  }

//! this generates a default map with linear distribution of points among
//! processes. Vectors based on this map will have halo points so that   
//! they can be used as either X or Y in Y=A*X operations.
void phist_map_create(map_ptr_t* vmap, const_comm_ptr_t vcomm, gidx_t nglob, int *ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(MPI_Comm,comm,vcomm,*ierr);

  ghost_map_t* map = new ghost_map_t;
  
  map->ctx=NULL;
//TODO: check ghost_err_t return codes everywhere like this:
  PHIST_CHK_GERR(ghost_context_create(&map->ctx,nglob, nglob, GHOST_CONTEXT_DEFAULT, NULL,GHOST_SPARSEMAT_SRC_NONE, *comm,1.0),*ierr);
  map->vtraits_template=phist_default_vtraits();
  // in ghost terminology, we look at LHS=A*RHS, the LHS is based on the
  // row distribution of A, the RHS has halo elements to allow importing from
  // neighbors. It is not possible to construct an RHS without a matrix, so if
  // we want to create one we need to get the 'map' from the matrix (e.g. the
  // domain map which is the same as the col map in ghost). A RHS vector can
  // be used as LHS, however.
  map->vtraits_template.flags = GHOST_DENSEMAT_NO_HALO;

  *vmap=(map_ptr_t)(map);
  }

//!
void phist_map_delete(map_ptr_t vmap, int *ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(ghost_map_t,map,vmap,*ierr);
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
  CAST_PTR_FROM_VOID(const ghost_map_t,map,vmap,*ierr);
  *vcomm = (const_comm_ptr_t)(&map->ctx->mpicomm);
  }

//!
void phist_map_get_local_length(const_map_ptr_t vmap, int* nloc, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(const ghost_map_t,map,vmap,*ierr);
  int me;
  ghost_rank(&me,map->ctx->mpicomm);
  *nloc=map->ctx->lnrows[me];
  }

//! returns the smallest global index in the map appearing on my partition.
void phist_map_get_ilower(const_map_ptr_t vmap, gidx_t* ilower, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(const ghost_map_t,map,vmap,*ierr);
  int me;
  ghost_rank(&me,map->ctx->mpicomm);
  *ilower = map->ctx->lfRow[me];
  }
//! returns the largest global index in the map appearing on my partition.
void phist_map_get_iupper(const_map_ptr_t vmap, gidx_t* iupper, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(const ghost_map_t,map,vmap,*ierr);
  int me;
  ghost_rank(&me,map->ctx->mpicomm);
  *iupper = map->ctx->lfRow[me+1]-1;
  }


} //extern "C"

// small helper function to preclude integer overflows (ghost allows 64 bit local indices, 
// but we don't right now)
  template<typename idx_t>
  int check_local_size(idx_t& i)
  {
    if (i>std::numeric_limits<lidx_t>::max()) return -1;
    return 0;
  }


#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "kernels_def.hpp"
#include "carp_def.hpp"

#include "phist_gen_c.h"
#include "kernels_def.hpp"
#include "carp_def.hpp"
#endif

#include "phist_gen_d.h"
#include "kernels_def.hpp"
#include "carp_def.hpp"

#include "phist_gen_z.h"
#include "kernels_def.hpp"
#include "carp_def.hpp"

