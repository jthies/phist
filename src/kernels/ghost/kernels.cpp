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

#include "phist_ghost_internal.h"
#include "phist_GhostMV.hpp"

#ifdef PHIST_HAVE_LIKWID
#include <likwid.h>
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

//! returns the local partition size based on a benchmark, the benchmark is run once when
//! this function is called first, subsequently the same weight will be returned unless
//! you set force_value>0. If force_vale<=0 is given, the benchmark is run anyway and the
//! newly measured value is returned in this and subsequent calls.
//! The resulting double can be passed as 'weight' parameter
//! to ghost_context_create (used in phist_map_create and sparseMat construction routines)
double get_proc_weight(double force_value)
{
  static double proc_weight=-1.0;
  int iflag;
  if (force_value>0.0) proc_weight=force_value;
  if (proc_weight<=0)
  {
    double max_bw, mean_bw;
    phist_bench_stream_triad(&mean_bw,&max_bw,&iflag);
    proc_weight=mean_bw*1.0e-9;
    if (iflag) proc_weight=1.0;
  }
  // do not return exactly 1.0 even if it was measured or an error occured on a single MPI process.
  // context_create below may hang if one process assumes it is in the default case of 1.0
  if (proc_weight==1.0) proc_weight+=1.0e-12;
  return proc_weight;
}


void get_C_sigma(int* C, int* sigma, int flags, MPI_Comm comm)
{
  // if the user sets both to postive values in the config file (via CMake), respect this choice
  // and do not override it by either flags or the presence of GPU processes. An exception is that
  // we strictly disallow reordering unless pHIST_SPARSEMAT_PERM_LOCAL is set
  static int C_stored=PHIST_SELL_C;
//  static int sigma_stored=PHIST_SELL_SIGMA;
  static int sigma_stored=1; // temporarily disable all sorting here until we merge the branch tests_with_sigma

  // only determine C and sigma once, then use these values subsequently for all maps/matrices. The
  // code below requires global reductions and we don't know if a user constructs many different matrices
  // in the course of a simulation.
  *C=C_stored;
  *sigma=sigma_stored;
  
  if (!(flags&PHIST_SPARSEMAT_PERM_LOCAL) ) *sigma=1;

  if (*C>0 && *sigma>0) return;

  if (*C<0)
  {
    if( flags & PHIST_SPARSEMAT_OPT_SINGLESPMVM )
    {
      *C = 32;
    }
    if( flags & PHIST_SPARSEMAT_OPT_BLOCKSPMVM )
    {
      *C = 8;
    }
  }

  // override with max(C,32) if anything runs on a CUDA device
  ghost_type gtype;
  ghost_type_get(&gtype);
  if (gtype==GHOST_TYPE_CUDA)
  {
    *C=std::max(*C,32);
  }
  // if the user doesn´t set it in CMake or give a flag, it is -1, override with +1 (CRS)
  *C=std::max(*C,+1);

  // everyone should have the max value found among MPI processes
  MPI_Allreduce(MPI_IN_PLACE,C,1,MPI_INT,MPI_MAX,comm);

  if (*sigma<0) *sigma=4*(*C);
  
  C_stored=*C;
  sigma_stored=*sigma;

}

int get_perm_flag(int iflag, int outlev)
{
  int oflag=GHOST_SPARSEMAT_DEFAULT;
  if (iflag&PHIST_SPARSEMAT_PERM_GLOBAL)
  {
#ifdef GHOST_HAVE_ZOLTAN
          PHIST_SOUT(outlev, "Trying to repartition the matrix with Zoltan\n");
          oflag|=GHOST_SPARSEMAT_ZOLTAN;
#elif defined(GHOST_HAVE_SCOTCH)
          PHIST_SOUT(outlev, "Trying to repartition the matrix with SCOTCH\n");
          oflag|=GHOST_SPARSEMAT_SCOTCHIFY;
#else
          PHIST_SOUT(outlev, "Zoltan or PT-Scotch not enabled in GHOST installation, no matrix repartitioning\n");
#endif
#ifdef GHOST_HAVE_SPMP
          PHIST_SOUT(outlev, "Enable local RCM reordering via SPMP\n");
          oflag|=GHOST_SPARSEMAT_RCM;
#endif
  }
  else if (iflag&PHIST_SPARSEMAT_PERM_LOCAL)
  {
#ifdef GHOST_HAVE_SPMP
          PHIST_SOUT(outlev, "Enable local RCM reordering via SPMP\n");
          oflag|=GHOST_SPARSEMAT_RCM;
#endif
  }
  if (iflag&PHIST_SPARSEMAT_OPT_CARP)
  {
        oflag|=GHOST_SOLVER_KACZ;
  }
  if (iflag&PHIST_SPARSEMAT_DIST2_COLOR)
  {
    oflag|=GHOST_SPARSEMAT_COLOR;
  }
  if (oflag!=GHOST_SPARSEMAT_DEFAULT) oflag|=GHOST_SPARSEMAT_PERMUTE;
  if (iflag&PHIST_SPARSEMAT_PERM_LOCAL == 0 )
  {
    PHIST_SOUT(PHIST_WARNING,"WARNING: based on your input flags, PHIST suggests to set permutation flags for the matrix.\n"
                             "         However, since PHIST_SPARSEMAT_PERM_LOCAL is missing from the input flags, I  can't set\n"
                             "         them. For optimal performance you should consider allowing at least local\n" 
                             "         permutations.\n");
    return 0;
  }
  return oflag;
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
      // do not allow GHOST to automatically allocate memory on the host CPU if a GPU process calls download().
      // If an mvec should be allocated on both host and devie, the PHIST user should specify
      // PHIST_REPLICATE_DEVICE_MEM when calling mvec_create.
      new_flags|=     (int)GHOST_DENSEMAT_NOT_RELOCATE;
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

    // this calls ghost_context_create with the given arguments and retries with a proc_weight of 1 if there are empty partitions
    void context_create(ghost_context **ctx, ghost_gidx gnrows, ghost_gidx gncols, 
        ghost_context_flags_t flags, void *matrixSource, ghost_sparsemat_src srcType, ghost_mpi_comm comm, double proc_weight, int* iflag)
    {
      *iflag=0;
      phist_gidx nglob_count;
      bool any_empty;
      *ctx=NULL;
      while (*ctx==NULL)
      {
        int rank,nproc;
        ghost_rank(&rank,comm);
        ghost_nrank(&nproc,comm);
        PHIST_ORDERED_OUT(PHIST_VERBOSE,comm,"PE%6d partition weight %4.2g\n",rank,proc_weight);
        PHIST_CHK_GERR(ghost_context_create(ctx,gnrows, gncols, flags, matrixSource,
            srcType, comm,proc_weight),*iflag);
        nglob_count=0;
        any_empty=false;
        for (int i=0;i<nproc;i++)
        {
          any_empty|=((*ctx)->lnrows[i]<=0);
          // do not count/compare number of rows if none was given
          if (gnrows!=0) nglob_count+=(*ctx)->lnrows[i];
        }
        if (nglob_count!=gnrows||(any_empty&&proc_weight!=1.0))
        {
          if (any_empty)          PHIST_SOUT(PHIST_WARNING,"empty partition in map/context\n");
          if (nglob_count!=gnrows) PHIST_SOUT(PHIST_WARNING,"number of nodes %ld in context does not match given nglob=%ld\n", (int64_t)nglob_count,(int64_t)gnrows);
          PHIST_SOUT(PHIST_WARNING,"GHOST did not give a correct context, %s\n",proc_weight==1.0?"aborting!":"retrying...");
          ghost_context_destroy(*ctx);
          *ctx=NULL;
          if (proc_weight==1.0)
          {
            *iflag=-1;
            break;
          }
          proc_weight=1.0;
        }
      }
      if (*iflag==0)
      {
        // set proc weight for subsequent contexts/maps
        get_proc_weight(proc_weight);
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
  PHIST_ORDERED_OUT(PHIST_VERBOSE,MPI_COMM_WORLD,"%s\n",str);
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

  // try to create the context object. We check wether
  // the weighting for the partitioning gives a feasible context. If the STREAM benchmark
  // returns very different values on the different MPI processes and the map is tiny (as
  // happens in our tests) the function may otherwise give empty partitions. In case the 
  // context looks strange we retry with equal process weights, accepting empty partitions.
  // If it still fails we return with an error code/message.
  double proc_weight=get_proc_weight();

  PHIST_CHK_IERR(phist::ghost_internal::context_create(&map->ctx,nglob, nglob, GHOST_CONTEXT_DEFAULT, NULL,
            GHOST_SPARSEMAT_SRC_NONE, *comm,proc_weight,iflag),*iflag);

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
  if (map->ctx) {
      ghost_rank(&me,map->ctx->mpicomm);
      *nloc=map->ctx->lnrows[me];
  } else {
      *nloc=map->vtraits_template.nrows;
  }
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


#ifdef GHOST_HAVE_CUDA
/* use GHOST stream benchmarks instead of those in common/ because the latter
   don't work on GPUs 
 */
#include "./bench_kernels.cpp"
#else
#include "../common/phist_bench_kernels.cpp"
#endif
