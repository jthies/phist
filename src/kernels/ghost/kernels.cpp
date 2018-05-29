/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_config.h"

#include <cstdio>
#include <cstdlib>

#include <iostream>
#include "phist_tools.h"
#include "../phist_kernels.h"
#include "phist_kernel_perfmodels.hpp"

#include "phist_typedefs.h"
#include "typedefs.hpp"

#include "phist_ghost_internal.h"

#include <ghost.h>
#include <ghost/machine.h>
#include <ghost/util.h>
#include <ghost/pumap.h>
#include <ghost/locality.h>
#include <ghost/timing.h>
#include <ghost/taskq.h>
#include <limits>
#include <map>

namespace phist 
{

  namespace ghost_internal
  {

//! returns a value to guide partitioning based on a benchmark of the memory bandwidth.  
//! The benchmark is run once when this function is called first, subsequently the same  
//! weight will be returned unless you set force_value>0. If force_vale<=0 is given, the 
//! benchmark is run anyway and the newly measured value is returned in this and         
//! subsequent calls. The resulting double can be passed as 'weight' parameter to ghost.
double get_proc_weight(double force_value)
{
  static double proc_weight=-1.0;
//  static double proc_weight=1.0; // variable partition sizes introduce lots of trouble with tests,
                                 // so we disable this feature for the moment (cf. issue #204)
  int iflag;
  if (force_value>0.0) proc_weight=force_value;
  if (proc_weight<=0)
  {
    double max_bw, mean_bw;
    phist_bench_stream_triad(&mean_bw,&max_bw,&iflag);
    proc_weight=mean_bw*1.0e-9;
    if (iflag) 
    {
      PHIST_OUT(PHIST_WARNING,"WARNING: stream benchmark for row distribution returned iflag=%d\n"
                              "(file %s, line %d). Setting proc_weight=1 on this rank!\n",
                              iflag,__FILE__,__LINE__);
      proc_weight=1.0;
    }
  }
  return proc_weight;
}


void get_C_sigma(int* C, int* sigma, int flags, MPI_Comm comm)
{
  // on the first call, figure out if there is a GPU process and set C=32 if so.
  // This requires global communication!!!
  static int any_GPUs=-1;
  if (any_GPUs<0)
  {
    ghost_type gtype;
    ghost_type_get(&gtype);
    if (gtype==GHOST_TYPE_CUDA)
    {
      any_GPUs=1;
    }
    else
    {
      any_GPUs=0;
    }

#ifdef PHIST_HAVE_MPI
    // everyone should have the max value found among MPI processes
    MPI_Allreduce(MPI_IN_PLACE,&any_GPUs,1,MPI_INT,MPI_SUM,comm);
#endif
  }

  // if the user sets both to postive values in the config file (via CMake), respect this choice
  // and do not override it by either flags or the presence of GPU processes. An exception is that
  // we strictly disallow reordering unless pHIST_SPARSEMAT_PERM_LOCAL is set, i.e. we set sigma=1
  // in that case.
  *C=PHIST_SELL_C;
  *sigma=PHIST_SELL_SIGMA;
  
  // workaround for #225: force sigma=1 unless the user provides a matrix (pair) himself, in which case
  // he or she is responsible for making them consistently ordered.
  if (*sigma!=1)
  {
    PHIST_SOUT(PHIST_WARNING,"To avoid possible inconsistent ordering between matrices, we force sigma=1\n"
                             " right now, see issue #225.\n");
    *sigma=+1;
  }
  
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
    if (any_GPUs)
    {
      *C=32;
    }
  }
  if (*C<0 || (flags&PHIST_SPARSEMAT_OPT_CARP)) 
  {
    *C=1;
    *sigma=1;
  }

  if (*sigma<0) *sigma=4*(*C);
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
//        PHIST_SOUT(outlev, "Enable reorderings for CARP kernel\n");
//        oflag|=GHOST_SOLVER_KACZ;
        PHIST_SOUT(PHIST_WARNING,"flag PHIST_SPARSEMAT_OPT_CARP ignored right now,\n"
                                 "since CARP is experimental in ghost and may mess up the matrix\n"
                                 "(file %s, line %d)\n",__FILE__,__LINE__);
  }
  if (iflag&PHIST_SPARSEMAT_DIST2_COLOR)
  {
    oflag|=GHOST_SPARSEMAT_COLOR;
  }

  if (oflag!=GHOST_SPARSEMAT_DEFAULT) oflag|=GHOST_SPARSEMAT_PERMUTE;
  if ( ((iflag&PHIST_SPARSEMAT_PERM_LOCAL)  == 0) && 
       ((iflag&PHIST_SPARSEMAT_PERM_GLOBAL) == 0) && 
        (oflag!=0) )
  {
    PHIST_SOUT(PHIST_WARNING,"WARNING: based on your input flags, PHIST suggests to set permutation flags for the matrix.\n"
                             "         However, since neither PHIST_SPARSEMAT_PERM_LOCAL nor PHIST_SPARSEMAT_PERM_GLOBAL\n"
                             "         are present in the input flags, I  can't set\n"
                             "         them. For optimal performance you should consider allowing at least local\n" 
                             "         permutations.\n");
    oflag=0;
  }


  if ( (iflag&PHIST_SPARSEMAT_OVERLAP_COMMUNICATION == 0) && 
       (PHIST_DEFAULT_SPMV_MODE==0) )
  {
    // save memory by not duplicating the matrix (otherwise it is stored as a whole and as local/non-local part)
    oflag|=GHOST_SPARSEMAT_NOT_STORE_SPLIT;
  }

  return oflag;
}

 
    //! private helper function to create a vtraits object
    ghost_densemat_traits default_vtraits()
    {
      ghost_densemat_traits vtraits = GHOST_DENSEMAT_TRAITS_INITIALIZER;
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
  int iflag_in=*iflag;
  bool quiet=(*iflag&PHIST_KERNELS_QUIET);
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

  if (std::getenv("OMP_SCHEDULE") == nullptr)
  {
    PHIST_OUT(PHIST_PERFWARNING,"GHOST uses OpenMP scheduling according to the environment variable 'OMP_SCHEDULE'.\n"
                                "We recommend setting it to e.g. 'static' because the default ('dynamic,1') usually\n"
                                "leads to rather poor performance.\n");
  }


  if (!quiet)
  {
    char *str = NULL;
    ghost_string(&str);
    PHIST_SOUT(PHIST_INFO,"%s\n",str);
    free(str); str = NULL;
    ghost_machine_string(&str);
    PHIST_SOUT(PHIST_INFO,"%s\n",str);
    free(str); str = NULL;

    int nnuma = 0;
    int npu = 0;
    ghost_machine_nnuma(&nnuma);
    ghost_machine_npu(&npu,GHOST_NUMANODE_ANY);

    ghost_pumap_string(&str);
    PHIST_ORDERED_OUT(PHIST_INFO,MPI_COMM_WORLD,"%s\n",str);
    free(str); str = NULL;
  }

  // initialize ghost's task-queue (required for tasks!)
  // (ghost does this also automatically when enqueuing the first task, but we might use
  //  som utility functions before!)
  PHIST_CHK_GERR(ghost_taskq_create(), *iflag);
  *iflag=iflag_in;
  PHIST_CHK_IERR(phist_kernels_common_init(argc,argv,iflag),*iflag);
}

// finalize ghost
extern "C" void phist_kernels_finalize(int* iflag)
{
  PHIST_CHK_IERR(phist_kernels_common_finalize(iflag),*iflag);
  bool quiet=(*iflag&PHIST_KERNELS_QUIET);
  if (!quiet)
  {
    char* ghostTimings = NULL;
    PHIST_CHK_GERR(ghost_timing_summarystring(&ghostTimings), *iflag);
    PHIST_SOUT(PHIST_INFO,"%s\n",ghostTimings);
    free(ghostTimings);
  }
  ghost_finalize();
  *iflag=0;
}


//! simply returns MPI_COMM_WORLD, the only MPI_Comm used in ghost.
extern "C" void phist_comm_create(phist_comm_ptr* vcomm, int* iflag)
{
  *iflag=0;
  // concept doesn't exist in ghost, return MPI_Comm
  MPI_Comm* comm = new MPI_Comm;
#ifdef PHIST_HAVE_MPI
  *comm=phist_get_default_comm();
#else
  *comm=(MPI_Comm)0;
#endif
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

  // create map object which is not part of a context.
  ghost_map* map=NULL;
  PHIST_CHK_GERR(ghost_map_create(&map,nglob,*comm,GHOST_MAP_ROW,GHOST_MAP_DEFAULT),*iflag);
  // allow distribution based on memory bandwidth benchmark:
  double weight=::phist::ghost_internal::get_proc_weight();
  int pad=std::max(32,PHIST_SELL_C);
  PHIST_CHK_GERR(ghost_map_create_distribution(map,NULL,weight,GHOST_MAP_DIST_NROWS,NULL),*iflag);
  
  // try to avoid empty partitions if the number of elements is small and the weights are unequal
  int me_empty=(map->dim==0)?1:0,any_empty;
#ifdef PHIST_HAVE_MPI
  PHIST_CHK_IERR(*iflag=MPI_Allreduce(&me_empty,&any_empty,1,MPI_INT,MPI_SUM,*comm),*iflag);
#else
  any_empty=me_empty;
#endif  
  if (any_empty)
  {
    weight=1.0;
    ghost_map_destroy(map);
    map=NULL;
    PHIST_CHK_GERR(ghost_map_create(&map,nglob,*comm,GHOST_MAP_ROW,GHOST_MAP_DEFAULT),*iflag);
    int pad=std::max(32,PHIST_SELL_C);
    PHIST_CHK_GERR(ghost_map_create_distribution(map,NULL,weight,GHOST_MAP_DIST_NROWS,NULL),*iflag);
  }

  map->dimpad=PAD(map->dimpad,pad);
  
  // in ghost terminology, we look at LHS=A*RHS, the LHS is based on the
  // row distribution of A, the RHS has halo elements to allow importing from
  // neighbors. It is not possible to construct an RHS without a matrix, so if
  // we want to create one we need to get the map from the matrix (e.g. the
  // domain map which is the same as the col map in ghost). An RHS vector can
  // be used as LHS, however.

  *vmap=(phist_map_ptr)(map);
}

//!
extern "C" void phist_map_delete(phist_map_ptr vmap, int *iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(ghost_map,map,vmap,*iflag);
  ghost_map_destroy(map);
  vmap=NULL;
}

//
extern "C" void phist_context_create(phist_context_ptr* vctx, 
                phist_const_map_ptr vrow_map,
                phist_const_map_ptr vcol_map,
                phist_const_map_ptr vrange_map, 
                phist_const_map_ptr vdomain_map, 
                int *iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_map,row_map,vrow_map,*iflag);
  ghost_map const* col_map=(ghost_map const*)vcol_map;
  ghost_map const* range_map=(ghost_map const*)vrange_map;
  ghost_map const* domain_map=(ghost_map const*)vdomain_map;
  
  if (domain_map==NULL) domain_map=row_map; // assume a square matrix
  if (range_map==NULL) 
  {
    range_map=row_map;  
  }
  // ghost doesn't allow range_map!=row_map!
  PHIST_CHK_IERR(*iflag=(range_map!=row_map)?PHIST_INVALID_INPUT:0,*iflag);
  
  ghost_context* ctx=NULL;
  
  ghost_gidx gnrows=row_map->gdim;
  ghost_gidx gncols=domain_map->gdim;
  ghost_context_create(&ctx,gnrows,gncols,GHOST_CONTEXT_DEFAULT,row_map->mpicomm,1.);
  
  PHIST_CHK_GERR(ghost_context_set_map(ctx,GHOST_MAP_ROW,(ghost_map*)row_map),*iflag);
  
  if (col_map!=NULL)
  {
    // this map will use the same communication as the matrix that originally created the maps
    // This means that the sparsity pattern must not include ghost nodes other than those in the
    // pattern of the matrix that originally created the context.
    PHIST_CHK_GERR(ghost_context_set_map(ctx,GHOST_MAP_COL,(ghost_map*)col_map),*iflag);
  }
  else
  {
    // the matrix has the same shape but may require different communication than the original one
    ghost_map_create_distribution(ctx->col_map,NULL,1.0,GHOST_MAP_DIST_NROWS,domain_map->ldim);
  }
  *vctx=(phist_context_ptr)ctx;
}

extern "C" void phist_context_delete(phist_context_ptr vctx, int* iflag)
{
  *iflag=0;
  PHIST_ENTER_FCN(__FUNCTION__);
  if (vctx==NULL) return;
  PHIST_CAST_PTR_FROM_VOID(ghost_context,ctx,vctx,*iflag);
  ghost_context_destroy(ctx);
  ctx=NULL;
}

//!
extern "C" void phist_map_get_comm(phist_const_map_ptr vmap, phist_const_comm_ptr* vcomm, int* iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_map,map,vmap,*iflag);
  *vcomm = (phist_const_comm_ptr)(&map->mpicomm);
}

//!
extern "C" void phist_map_get_local_length(phist_const_map_ptr vmap, phist_lidx* nloc, int* iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_map,map,vmap,*iflag);
   *nloc=map->dim;
}

//!
extern "C" void phist_map_get_global_length(phist_const_map_ptr vmap, phist_gidx* nglob, int* iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_map,map,vmap,*iflag);
  *nglob=map->gdim;
}

//! returns the smallest global index in the map appearing on my partition.
extern "C" void phist_map_get_ilower(phist_const_map_ptr vmap, phist_gidx* ilower, int* iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_map,map,vmap,*iflag);
  *ilower = map->offs;
}
//! returns the largest global index in the map appearing on my partition.
extern "C" void phist_map_get_iupper(phist_const_map_ptr vmap, phist_gidx* iupper, int* iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_map,map,vmap,*iflag);
  *iupper = map->offs + map->dim - 1;
}

extern "C" void phist_maps_compatible(phist_const_map_ptr vmap1, phist_const_map_ptr vmap2, int* iflag)
{
static bool first_time=true;
if (first_time)
{
  PHIST_SOUT(PHIST_WARNING,"phist_maps_compatible is *NOT IMPLEMENTED!* for GHOST 1.1, returning 0 anyway for development purposes.\n");
  first_time=false;
}
  *iflag=0;
  return;
  #if FIX_ME
  PHIST_CAST_PTR_FROM_VOID(const phist_ghost_map,map1,vmap1,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const phist_ghost_map,map2,vmap2,*iflag);
  *iflag=-1;
  // same object?
  if (map1==map2) {*iflag=0; return;}
  
  const ghost_densemat_traits& vtraits1 = map1->vtraits_template;
  const ghost_densemat_traits& vtraits2 = map2->vtraits_template;
  
  int permuted1 = vtraits1.flags&GHOST_DENSEMAT_PERMUTED;
  int permuted2 = vtraits2.flags&GHOST_DENSEMAT_PERMUTED;

  const ghost_context *ctx1 = map1->ctx;  
  const ghost_context *ctx2 = map2->ctx;

  const ghost_lidx *lperm1 = map1->perm_local;
  const ghost_lidx *liperm1 = map1->perm_local_inv;
  const ghost_lidx *lperm2 = map2->perm_local;
  const ghost_lidx *liperm2 = map2->perm_local_inv;
  const ghost_gidx *gperm1 = map1->perm_global;
  const ghost_gidx *giperm1 = map1->perm_global_inv;
  const ghost_gidx *gperm2 = map2->perm_global;
  const ghost_gidx *giperm2 = map2->perm_global_inv;
  
  // compare contexts and permutation objects as far as we need the info to be consistent
  if (ctx1!=ctx2)
  {
    if ( ctx1->row_map->gdim  != ctx2->row_map->gdim  ||
         ctx1->mpicomm != ctx2->mpicomm )
    {
      *iflag=-1;
      return;
    }
  }
  
  // check if the two maps contain the same vector permutation info.
  // If not, this is not a game breaker: the maps are still compatible
  // if at most one of them is permuted.
  bool same_lperm=true;
  bool same_gperm=true;
  
  if (lperm1!=lperm2)
  {
    if ( lperm1!=lperm2        ||
       liperm1!=liperm2 )
    {
      same_lperm=false;
    }
  }

  if (gperm1!=gperm2)
  {
    if ( gperm1!=gperm2        ||
       giperm1!=giperm2 )
    {
      same_gperm=false;
    }
  }

  if (!(permuted1||permuted2))
  {
    *iflag=0;
    return;
  }

  if ( (same_lperm==false || same_gperm==false) &&
       (permuted1&&permuted2) )
  {
    *iflag=-1;
    return;
  }
    
  // contexts, permutations and traits/flags are the same
  if (&vtraits1==&vtraits2 && same_lperm && same_gperm)
  {
    // this means if the vectors are permuted, they are permuted the same way and do
    // not need to be permuted via mvec_to_mvec in order to be compatible.
    *iflag=0;
    return;
  }

  // since ghost_vtraits doesn't implement operator==, we have to compare manually
  if (map1->map->dim == map2->map->dim &&
      permuted1      == permuted2      &&
      same_lperm     && same_gperm)
  {
    // these maps will create vectors with the same distribution and permutation
    *iflag=0;
    return;
  }
  
  // one of the vectors is permuted, the other isn't. Check if the permutation is global or local
  if (gperm1!=NULL || gperm2!=NULL) 
  {
    *iflag=2;
  }
  else if (lperm1!=NULL || lperm2!=NULL)
  {
    *iflag=1;
  }
  else
  {
    // I think we never should get to this point!
    *iflag=-2;
  }
  return;
#endif
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


/* use GHOST stream benchmarks instead of those in common/ because the latter
   don't work on GPUs and tended to fail with GHOST even without CUDA
 */
#include "./bench_kernels.cpp"
/*
#include "../common/phist_bench_kernels.cpp"
*/
