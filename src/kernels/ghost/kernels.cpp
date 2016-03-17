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
# if defined(GHOST_IDX64_LOCAL)
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
  // if the user sets both to postive values, respect this choice
  // and do not override it by either flags or the presence of GPU
  // processes
  if (*C>0 && *sigma>0) return;

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
  ghost_type gtype;
  ghost_type_get(&gtype);
  if (gtype==GHOST_TYPE_CUDA)
  {
    *C=std::max(*C,32);
    *sigma=std::max(256,*sigma);
  }
  // if the user doesn´t set it in CMake or give a flag, it is -1, override with +1 (CRS)
  *C=std::max(*C,+1);
  *sigma=std::max(*sigma,+1);
  // everyone should have the max value found among MPI processes
  MPI_Allreduce(MPI_IN_PLACE,C,1,MPI_INT,MPI_MAX,comm);
  MPI_Allreduce(MPI_IN_PLACE,sigma,1,MPI_INT,MPI_MAX,comm);

}

int get_perm_flag(int outlev)
{
#ifdef USE_ZOLTAN
          PHIST_SOUT(outlev, "Trying to repartition the matrix with Zoltan\n");
          return GHOST_SPARSEMAT_ZOLTAN;
#elif defined(USE_SCOTCH)
          PHIST_SOUT(outlev, "Trying to repartition the matrix with SCOTCH\n");
          return GHOST_SPARSEMAT_SCOTCHIFY;
#else
          PHIST_SOUT(outlev, "SCOTCH not available, no matrix repartitioning\n");
          return 0;
#endif
}

 
    //! private helper function to create a vtraits object
    ghost_densemat_traits phist_default_vtraits()
    {
      ghost_densemat_traits vtraits = GHOST_DENSEMAT_TRAITS_INITIALIZER;
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
      vtraits.flags = (ghost_densemat_flags)new_flags;

      vtraits.ncols=1;
#ifdef PHIST_MVECS_ROW_MAJOR
      vtraits.storage=GHOST_DENSEMAT_ROWMAJOR;
#else
      vtraits.storage=GHOST_DENSEMAT_COLMAJOR;
#endif
      return vtraits;
    }  
    
    ghost_map* MapGarbageCollector::new_map(const void* p)
    {
      ghost_map* m = new ghost_map;
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
  ghost_hwconfig hwconfig = GHOST_HWCONFIG_INITIALIZER; 
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

  ghost_thpool *thpool;
  int nnuma = 0;
  int npu = 0;

  ghost_thpool_get(&thpool);
  ghost_machine_nnuma(&nnuma);
  ghost_machine_npu(&npu,GHOST_NUMANODE_ANY);

  ghost_pumap_string(&str);
  PHIST_ORDERED_OUT(PHIST_VERBOSE,"%s\n",str);
  free(str); str = NULL;

#if defined(PHIST_HAVE_KOKKOS)&&defined(PHIST_HAVE_BELOS)
   PHIST_SOUT(PHIST_INFO,"TSQR using node-type %s\n",node_type::name().c_str());
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
extern "C" void phist_comm_create(phist_comm_ptr* vcomm, int* iflag)
  {
  *iflag=0;
  // concept doesn't exist in ghost, return MPI_Comm
  MPI_Comm* comm = new MPI_Comm;
  *comm=MPI_COMM_WORLD;
  *vcomm=(phist_comm_ptr)comm;
  }

//!
extern "C" void phist_comm_delete(phist_comm_ptr vcomm, int* iflag)
  {
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(MPI_Comm,comm,vcomm,*iflag);
  delete comm;
  }

//!
extern "C" void phist_comm_get_rank(phist_const_comm_ptr vcomm, int* rank, int* iflag)
  {
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(MPI_Comm,comm,vcomm,*iflag);
  ghost_rank(rank,*comm);
  }
//!
extern "C" void phist_comm_get_size(phist_const_comm_ptr vcomm, int* size, int* iflag)
  {
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(MPI_Comm,comm,vcomm,*iflag);
  ghost_nrank(size,*comm);
  }
#ifdef PHIST_HAVE_MPI
extern "C" void phist_comm_get_mpi_comm(phist_const_comm_ptr vcomm, MPI_Comm* mpiComm, int* iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(MPI_Comm,comm,vcomm,*iflag);
  *mpiComm = *comm;
}
#endif
//! this generates a default map with linear distribution of points among
//! processes. Vectors based on this map will have halo points so that   
//! they can be used as either X or Y in Y=A*X operations.
extern "C" void phist_map_create(phist_map_ptr* vmap, phist_const_comm_ptr vcomm, phist_gidx nglob, int *iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(MPI_Comm,comm,vcomm,*iflag);

  ghost_map* map = new ghost_map;
  
  map->ctx=NULL;
  // passing in 0.0 here will lead to automatic load-balancing depending on the STREAM benchmark performance on each MPI 
  // rank. We certainly want to expose this cool feature to the user, but if we do it here a costly benchmark is run 
  // whenever a map is created. We probably want to create a static context object and reuse it instead.
  double proc_weight=1.0;
  PHIST_CHK_GERR(ghost_context_create(&map->ctx,nglob, nglob, GHOST_CONTEXT_DEFAULT, NULL,
        GHOST_SPARSEMAT_SRC_NONE, *comm,proc_weight),*iflag);
  map->vtraits_template=phist_default_vtraits();
  // in ghost terminology, we look at LHS=A*RHS, the LHS is based on the
  // row distribution of A, the RHS has halo elements to allow importing from
  // neighbors. It is not possible to construct an RHS without a matrix, so if
  // we want to create one we need to get the 'map' from the matrix (e.g. the
  // domain map which is the same as the col map in ghost). A RHS vector can
  // be used as LHS, however.
  int new_flags=(int)map->vtraits_template.flags;
  new_flags |= (int)GHOST_DENSEMAT_NO_HALO;
  map->vtraits_template.flags=(ghost_densemat_flags)new_flags;
  // ghost should set these correctly depending on GHOST_TYPE if we set HOST and DEVICE to 0

  map->permutation=NULL; //permutation is defined by a matrix object.

  *vmap=(phist_map_ptr)(map);
}

//!
extern "C" void phist_map_delete(phist_map_ptr vmap, int *iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(ghost_map,map,vmap,*iflag);
  // should be safe if calling order is respected? (e.g. create map, create stuff with map, destroy stuff, destroy map)
  ghost_context_destroy(map->ctx);
  delete map;
  vmap=NULL;
}
  
//!
extern "C" void phist_map_get_comm(phist_const_map_ptr vmap, phist_const_comm_ptr* vcomm, int* iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_map,map,vmap,*iflag);
  *vcomm = (phist_const_comm_ptr)(&map->ctx->mpicomm);
}

//!
extern "C" void phist_map_get_local_length(phist_const_map_ptr vmap, phist_lidx* nloc, int* iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_map,map,vmap,*iflag);
  int me;
  ghost_rank(&me,map->ctx->mpicomm);
  *nloc=map->ctx->lnrows[me];
}

//!
extern "C" void phist_map_get_global_length(phist_const_map_ptr vmap, phist_gidx* nglob, int* iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_map,map,vmap,*iflag);
  int me;
  *nglob=map->ctx->gnrows;
}

//! returns the smallest global index in the map appearing on my partition.
extern "C" void phist_map_get_ilower(phist_const_map_ptr vmap, phist_gidx* ilower, int* iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_map,map,vmap,*iflag);
  int me;
  ghost_rank(&me,map->ctx->mpicomm);
  *ilower = map->ctx->lfRow[me];
}
//! returns the largest global index in the map appearing on my partition.
extern "C" void phist_map_get_iupper(phist_const_map_ptr vmap, phist_gidx* iupper, int* iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_map,map,vmap,*iflag);
  int me;
  ghost_rank(&me,map->ctx->mpicomm);
  *iupper = map->ctx->lfRow[me]+map->ctx->lnrows[me]-1;
}

/* helper function that tells the matrix generation routines which 
   repartitioning flags they should pass to GHOST (depending on how
   GHOST was compiled/installed different options may be available)
   
   The order we try is Zoltan/Scotch/nothing
   
 */
   
#ifdef GHOST_HAVE_ZOLTAN
#define USE_ZOLTAN
#elif defined(GHOST_HAVE_SCOTCH)
#define USE_SCOTCH
#endif


/* use GHOST stream benchmarks instead of those in common/ */
#include "./bench_kernels.cpp"
