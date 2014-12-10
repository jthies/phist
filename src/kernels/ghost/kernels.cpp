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

#ifdef PHIST_HAVE_ANASAZI
#include "phist_trilinos_macros.h"
#include "Anasazi_GhostAdapter.hpp"
#include "phist_AnasaziOperatorTraits.hpp"
#include "AnasaziSVQBOrthoManager.hpp"
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
#include <ghost/thpool.h>
#include <ghost/pumap.h>
#include <ghost/locality.h>
#include <limits>
#include "phist_ghost_macros.hpp"

#if defined(PHIST_HAVE_BELOS)||defined(PHIST_HAVE_KOKKOS)
# if (!defined(GHOST_HAVE_LONGIDX_GLOBAL)) || defined(GHOST_HAVE_LONGIDX_LOCAL)
# warning "The interfaces between GHOST and Belos/TSQR cause problems unless you compile GHOST with LONGIDX_GLOBAL but *without* LONGIDX_LOCAL"
# endif
#endif

namespace phist
  {
  int GhostMV::countObjects=0;
  }

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
  ghost_permutation_t *permutation;
  } ghost_map_t;

// initialize ghost
extern "C" void phist_kernels_init(int* argc, char*** argv, int* ierr)
{
  *ierr=0;
  ghost_init(*argc, *argv);

  char *str = NULL;
  ghost_string(&str);
  PHIST_SOUT(PHIST_INFO,"%s\n",str);
  free(str); str = NULL;

  ghost_machine_string(&str);
  PHIST_SOUT(PHIST_INFO,"%s\n",str);
  free(str); str = NULL;

  ghost_thpool_t *thpool;
  int nnuma = 0;
  int npu = 0;

  ghost_thpool_get(&thpool);
  ghost_machine_nnuma(&nnuma);
  ghost_machine_npu(&npu,GHOST_NUMANODE_ANY);

  ghost_pumap_string(&str);
  PHIST_SOUT(PHIST_VERBOSE,"%s\n",str);
  free(str); str = NULL;

  PHIST_SOUT(PHIST_VERBOSE,"The thread pool consists of %d threads\n",thpool->nThreads);

#ifdef PHIST_HAVE_LIKWID
  LIKWID_MARKER_INIT;
#pragma omp parallel
  {
    LIKWID_MARKER_THREADINIT;
    LIKWID_MARKER_START("phist<ghost>");
  }
#endif
}

// finalize ghost
extern "C" void phist_kernels_finalize(int* ierr)
{
#ifdef PHIST_HAVE_LIKWID
#pragma omp parallel
  {
    LIKWID_MARKER_STOP("phist<ghost>");
  }
  LIKWID_MARKER_CLOSE;
#endif
PHIST_CXX_TIMER_SUMMARIZE;
  ghost_finalize();
  *ierr=0;
}


//! simply returns MPI_COMM_WORLD, the only MPI_Comm used in ghost.
extern "C" void phist_comm_create(comm_ptr_t* vcomm, int* ierr)
  {
  *ierr=0;
  // concept doesn't exist in ghost, return MPI_Comm
  MPI_Comm* comm = new MPI_Comm;
  *comm=MPI_COMM_WORLD;
  *vcomm=(comm_ptr_t)comm;
  }

//!
extern "C" void phist_comm_delete(comm_ptr_t vcomm, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(MPI_Comm,comm,vcomm,*ierr);
  delete comm;
  }

//!
extern "C" void phist_comm_get_rank(const_comm_ptr_t vcomm, int* rank, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(MPI_Comm,comm,vcomm,*ierr);
  ghost_rank(rank,*comm);
  }
//!
extern "C" void phist_comm_get_size(const_comm_ptr_t vcomm, int* size, int* ierr)
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
  vtraits.nrowsorig=0; // get from context
  vtraits.nrowshalo=0; // get from context
  vtraits.nrowspadded=0; // leave padding to ghost
  vtraits.ncols=0; // set in mvec_create
  vtraits.ncolsorig=0; // set in mvec_create
  vtraits.ncolspadded=0; // leave padding to ghost
  // ghost should set these correctly depending on GHOST_TYPE if we set HOST and DEVICE to 0
  int new_flags =   (int)vtraits.flags;
 //     new_flags|=    (int)GHOST_DENSEMAT_NO_HALO;
      new_flags&=   ~(int)GHOST_DENSEMAT_HOST;
      new_flags&=   ~(int)GHOST_DENSEMAT_DEVICE;
  vtraits.flags = (ghost_densemat_flags_t)new_flags;

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
extern "C" void phist_map_create(map_ptr_t* vmap, const_comm_ptr_t vcomm, gidx_t nglob, int *ierr)
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
  int new_flags=(int)map->vtraits_template.flags;
  new_flags |= (int)GHOST_DENSEMAT_NO_HALO;
  map->vtraits_template.flags=(ghost_densemat_flags_t)new_flags;
  // ghost should set these correctly depending on GHOST_TYPE if we set HOST and DEVICE to 0

  map->permutation=NULL; //permutation is defined by a matrix object.

  *vmap=(map_ptr_t)(map);
  }

//!
extern "C" void phist_map_delete(map_ptr_t vmap, int *ierr)
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
extern "C" void phist_map_get_comm(const_map_ptr_t vmap, const_comm_ptr_t* vcomm, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(const ghost_map_t,map,vmap,*ierr);
  *vcomm = (const_comm_ptr_t)(&map->ctx->mpicomm);
  }

//!
extern "C" void phist_map_get_local_length(const_map_ptr_t vmap, lidx_t* nloc, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(const ghost_map_t,map,vmap,*ierr);
  int me;
  ghost_rank(&me,map->ctx->mpicomm);
  *nloc=map->ctx->lnrows[me];
  }

//! returns the smallest global index in the map appearing on my partition.
extern "C" void phist_map_get_ilower(const_map_ptr_t vmap, gidx_t* ilower, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(const ghost_map_t,map,vmap,*ierr);
  int me;
  ghost_rank(&me,map->ctx->mpicomm);
  *ilower = map->ctx->lfRow[me];
  }
//! returns the largest global index in the map appearing on my partition.
extern "C" void phist_map_get_iupper(const_map_ptr_t vmap, gidx_t* iupper, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(const ghost_map_t,map,vmap,*ierr);
  int me;
  ghost_rank(&me,map->ctx->mpicomm);
  *iupper = map->ctx->lfRow[me]+map->ctx->lnrows[me]-1;
  }

// small helper function to preclude integer overflows (ghost allows 64 bit local indices, 
// but we don't right now)
  template<typename idx_t>
  int check_local_size(idx_t& i)
  {
    if (i>std::numeric_limits<lidx_t>::max()) return -1;
    return 0;
  }

#ifdef PHIST_TIMEMONITOR
extern "C" void phist_totalMatVecCount()
{
  ENTER_FCN(__FUNCTION__);
}
#endif

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

