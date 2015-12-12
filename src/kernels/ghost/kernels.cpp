#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <cstdio>
#include <cstdlib>

#include <iostream>
#include "phist_macros.h"
#include "../phist_kernels.h"
#include "phist_kernel_perfmodels.hpp"

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

#include "phist_ghost_internal.h"
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
#include <ghost/timing.h>
#include <limits>
#include <map>

#if defined(PHIST_HAVE_BELOS)||defined(PHIST_HAVE_KOKKOS)
# if defined(GHOST_HAVE_LONGIDX_LOCAL)
# warning "The interfaces between GHOST and Belos/TSQR cause problems unless you compile GHOST with LONGIDX_GLOBAL but *without* LONGIDX_LOCAL"
# endif
#endif

namespace phist 
{

  int GhostMV::countObjects=0;

  namespace ghost_internal
  {

void get_C_sigma(int* C, int* sigma, int flags, MPI_Comm comm)
{
  *C = PHIST_SELL_C;
  *sigma = PHIST_SELL_SIGMA;
  if( flags & PHIST_SPARSEMAT_OPT_SINGLESPMVM )
  {
    *C = 32;
    *sigma = 256;
  }
  if( flags & PHIST_SPARSEMAT_OPT_BLOCKSPMVM )
  {
    *C = 8;
    *sigma = 32;
  }

  // override with max(C,32) if anything runs on a CUDA device
  ghost_type_t gtype;
  ghost_type_get(&gtype);
  if (gtype==GHOST_TYPE_CUDA)
  {
    *C=std::max(*C,32);
//    *sigma=std::max(256,*sigma);
  }
  MPI_Allreduce(MPI_IN_PLACE,C,1,MPI_INT,MPI_MAX,comm);
  MPI_Allreduce(MPI_IN_PLACE,sigma,1,MPI_INT,MPI_MAX,comm);

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
 //     new_flags&=   ~(int)GHOST_DENSEMAT_HOST;
 //     new_flags&=   ~(int)GHOST_DENSEMAT_DEVICE;
      vtraits.flags = (ghost_densemat_flags_t)new_flags;

      vtraits.ncols=1;
#ifdef PHIST_MVECS_ROW_MAJOR
      vtraits.storage=GHOST_DENSEMAT_ROWMAJOR;
#else
      vtraits.storage=GHOST_DENSEMAT_COLMAJOR;
#endif
      return vtraits;
    }  
    
    ghost_map_t* MapGarbageCollector::new_map(const void* p)
    {
        ghost_map_t* m = new ghost_map_t;
        maps_[p].push_back(m);
        return m;
    }

    void MapGarbageCollector::delete_maps(void* p)
    {
      MapCollection::iterator it = maps_.find(p);
      if( it != maps_.end() )
      {
        for(int i = 0; i < it->second.size(); i++)
        {
          delete it->second[i];
        }
        maps_.erase(it);
      }
    }
    
  }//namespace ghost_internal
}//namespace phist

using namespace phist::ghost_internal;

//phist::ghost_internal::MapGarbageCollector phist::ghost_internal::mapGarbageCollector;

////////////////////////////////////////////////////////////
// public phist kernel interface implementation for ghost //
////////////////////////////////////////////////////////////

// initialize ghost
extern "C" void phist_kernels_init(int* argc, char*** argv, int* iflag)
{
  *iflag=0;
  // disable using Hyperthreads
  ghost_hwconfig_t hwconfig = GHOST_HWCONFIG_INITIALIZER; 
  const char* PHIST_NUM_THREADS=getenv("PHIST_NUM_THREADS");
  int num_threads= (PHIST_NUM_THREADS==NULL)? -1:atoi(PHIST_NUM_THREADS);

  // if the user sets both PHIST_NUM_THREADS and GHOST_TYPE=cuda,
  // avoid hanging ghost tasks because this process will get only one core.
  // So let ghost handle this situation even if PHIST_NUM_THREADS is set.
  // This will not handle the situation correctly where multiple MPI processes
  // are started on a node and rank 1 gets the GPU automatically, I think (cf. #128)
  const char* GHOST_TYPE=getenv("GHOST_TYPE");
  if (GHOST_TYPE && !strcasecmp(GHOST_TYPE,"cuda")) num_threads=-1;
  if (num_threads!=-1)
  {
    hwconfig.ncore = num_threads; 
  }
#ifdef PHIST_BUILD_MIC
  hwconfig.nsmt = 3; 
#else
  hwconfig.nsmt = 1; 
#endif
  ghost_hwconfig_set(hwconfig);

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

  PHIST_SOUT(PHIST_INFO,"The thread pool consists of %d threads\n",thpool->nThreads);

#if defined(PHIST_HAVE_KOKKOS)&&defined(PHIST_HAVE_BELOS)
   PHIST_SOUT(PHIST_INFO,"TSQR using node-type %s\n",node_t::name().c_str());
#else
   PHIST_SOUT(PHIST_INFO,"TSQR not available\n");
#endif

  PHIST_CHK_IERR(phist_kernels_common_init(argc,argv,iflag),*iflag);
}

// finalize ghost
extern "C" void phist_kernels_finalize(int* iflag)
{
  PHIST_CHK_IERR(phist_kernels_common_finalize(iflag),*iflag);
  char* ghostTimings = NULL;
  PHIST_CHK_GERR(ghost_timing_summarystring(&ghostTimings), *iflag);
  PHIST_SOUT(PHIST_INFO,"%s\n",ghostTimings);
  free(ghostTimings);
  ghost_finalize();
  *iflag=0;
}


//! simply returns MPI_COMM_WORLD, the only MPI_Comm used in ghost.
extern "C" void phist_comm_create(comm_ptr_t* vcomm, int* iflag)
  {
  *iflag=0;
  // concept doesn't exist in ghost, return MPI_Comm
  MPI_Comm* comm = new MPI_Comm;
  *comm=MPI_COMM_WORLD;
  *vcomm=(comm_ptr_t)comm;
  }

//!
extern "C" void phist_comm_delete(comm_ptr_t vcomm, int* iflag)
  {
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(MPI_Comm,comm,vcomm,*iflag);
  delete comm;
  }

//!
extern "C" void phist_comm_get_rank(const_comm_ptr_t vcomm, int* rank, int* iflag)
  {
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(MPI_Comm,comm,vcomm,*iflag);
  ghost_rank(rank,*comm);
  }
//!
extern "C" void phist_comm_get_size(const_comm_ptr_t vcomm, int* size, int* iflag)
  {
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(MPI_Comm,comm,vcomm,*iflag);
  ghost_nrank(size,*comm);
  }
#ifdef PHIST_HAVE_MPI
extern "C" void phist_comm_get_mpi_comm(const_comm_ptr_t vcomm, MPI_Comm* mpiComm, int* iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(MPI_Comm,comm,vcomm,*iflag);
  *mpiComm = *comm;
}
#endif
//! this generates a default map with linear distribution of points among
//! processes. Vectors based on this map will have halo points so that   
//! they can be used as either X or Y in Y=A*X operations.
extern "C" void phist_map_create(map_ptr_t* vmap, const_comm_ptr_t vcomm, gidx_t nglob, int *iflag)
  {
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(MPI_Comm,comm,vcomm,*iflag);

  ghost_map_t* map = new ghost_map_t;
  
  map->ctx=NULL;
//TODO: check ghost_err_t return codes everywhere like this:
  PHIST_CHK_GERR(ghost_context_create(&map->ctx,nglob, nglob, GHOST_CONTEXT_DEFAULT, NULL,GHOST_SPARSEMAT_SRC_NONE, *comm,1.0),*iflag);
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
extern "C" void phist_map_delete(map_ptr_t vmap, int *iflag)
  {
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(ghost_map_t,map,vmap,*iflag);
  // should be safe if calling order is respected? (e.g. create map, create stuff with map, destroy stuff, destroy map)
  ghost_context_destroy(map->ctx);
  delete map;
  vmap=NULL;
  }
  
//!
extern "C" void phist_map_get_comm(const_map_ptr_t vmap, const_comm_ptr_t* vcomm, int* iflag)
  {
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_map_t,map,vmap,*iflag);
  *vcomm = (const_comm_ptr_t)(&map->ctx->mpicomm);
  }

//!
extern "C" void phist_map_get_local_length(const_map_ptr_t vmap, lidx_t* nloc, int* iflag)
  {
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_map_t,map,vmap,*iflag);
  int me;
  ghost_rank(&me,map->ctx->mpicomm);
  *nloc=map->ctx->lnrows[me];
  }

//!
extern "C" void phist_map_get_global_length(const_map_ptr_t vmap, gidx_t* nglob, int* iflag)
  {
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_map_t,map,vmap,*iflag);
  int me;
  *nglob=map->ctx->gnrows;
  }

//! returns the smallest global index in the map appearing on my partition.
extern "C" void phist_map_get_ilower(const_map_ptr_t vmap, gidx_t* ilower, int* iflag)
  {
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_map_t,map,vmap,*iflag);
  int me;
  ghost_rank(&me,map->ctx->mpicomm);
  *ilower = map->ctx->lfRow[me];
  }
//! returns the largest global index in the map appearing on my partition.
extern "C" void phist_map_get_iupper(const_map_ptr_t vmap, gidx_t* iupper, int* iflag)
  {
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_map_t,map,vmap,*iflag);
  int me;
  ghost_rank(&me,map->ctx->mpicomm);
  *iupper = map->ctx->lfRow[me]+map->ctx->lnrows[me]-1;
  }


