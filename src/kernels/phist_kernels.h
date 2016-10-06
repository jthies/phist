#ifndef PHIST_KERNELS_H
#define PHIST_KERNELS_H

#include "phist_config.h"
#include "phist_kernel_flags.h"

#ifndef DOXYGEN
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#else
typedef int MPI_Comm;
#endif
#include "phist_macros.h"

// note: the phist_typedefs.h file is provided in the subdirectory
// where the interface is implemented (e.g. ghost/, tpetra/).
#include "phist_typedefs.h"
#include "phist_tasks.h"

#ifdef PHIST_HAVE_GHOST
#include <ghost/types.h>
#include <ghost/sparsemat.h>
#else
typedef phist_lidx ghost_lidx;
typedef phist_gidx ghost_gidx;
#endif

#endif

//! \defgroup kernels Kernel function interface
//@{

//! @file phist_kernels.h
//! @brief Definition of an abstract interface for basic operations needed by iterative solvers.
//!
//! The kernels are based on a number of abstract concepts, which are
//! simply void* into a "kernel library":
//! - phist_comm: encapsulates MPI, currently the interface does not really
//!   allow you to do much with this object, except use it to build other
//!   objects.
//! - phist_map: encapsulates the distribution of points (e.g. rows of a matrix)
//!   over processors
//! - multi-vector (mvec): a dense matrix with N rows and m<<N columns, with a row-wise
//!   distribution over computational nodes. The storage may be either column-major or
//!   or row-major. In the latter case, the macro PHIST_MVECS_ROW_MAJOR must be 
//!   defined in the phist_config.h file.
//! - sdMat - a small dense matrix, replicated on all MPI processes. sdMats are always
//!   stored in column-major order but may have a leading dimension larger than the number of rows.
//! - sparseMat: sparse matrix in *any* storage format
//!
//! The type-specific objects for mvecs, sdMats and sparseMats are called
//! phist_Smvec, phist_Dmvec, phist_Cmvec etc. and their related functions defined in 
//! phist_kernels_decl.h

#ifdef __cplusplus
extern "C" {
#endif

//! \name data-type independent kernel functions
//!@{

//! initialize kernel library. Should at least call MPI_Init if it has not been called
//! but is required.
void phist_kernels_init(int *argc, char*** argv, int* iflag);

//! finalize kernel library. Should at least call MPI_Finalize if it has not been called
//! but is required.
void phist_kernels_finalize(int* iflag);

//! creates a global comm object
void phist_comm_create(phist_comm_ptr* comm, int* iflag);
//! delete a comm object. Only do this for comms obtained by phist_comm_create.
void phist_comm_delete(phist_comm_ptr comm, int* iflag);
//! get the rank of the calling node
void phist_comm_get_rank(phist_const_comm_ptr comm, int* rank, int* iflag);
//! get the number of MPI asks
void phist_comm_get_size(phist_const_comm_ptr comm, int* size, int* iflag);
//! get MPI comm. If the kernel lib does not use MPI, NULL is returned.
void phist_comm_get_mpi_comm(phist_const_comm_ptr comm, MPI_Comm* mpiComm,int* iflag);

//! creates a map with default distribution of points
void phist_map_create(phist_map_ptr* map, phist_const_comm_ptr comm, phist_gidx nglob, int *iflag);
//! delete a map object. Note that you should not do this if you got the map from
//! anything else than phist_map_create.
void phist_map_delete(phist_map_ptr map, int *iflag);
//! returns the comm object used by a map
void phist_map_get_comm(phist_const_map_ptr map, phist_const_comm_ptr* comm, int* iflag);
//! returns the local number of elements in the map
void phist_map_get_local_length(phist_const_map_ptr map, phist_lidx* nloc, int* iflag);
//! returns the global number of elements in the map.
void phist_map_get_global_length(phist_const_map_ptr map, phist_gidx* nglob, int* iflag);
//! returns the smallest global index in the map appearing on my partition. iflag is set to 1
//! in case the map is not contiguous, because in that case it may be that the
//! caller falsely assumes global elements [ilower ... iupper] are actually on this partition.
void phist_map_get_ilower(phist_const_map_ptr map, phist_gidx* ilower, int* iflag);
//! returns the largest global index in the map appearing on my partition. iflag is set to 1
//! in case the map is not contiguous, because in that case it may be that the
//! caller falsely assumes global elements [ilower ... iupper] are actually on this partition.
void phist_map_get_iupper(phist_const_map_ptr map, phist_gidx* iupper, int* iflag);
//! set *iflag= 0 if the two maps represent the same distribution and permutation, 
//!     *iflag=+1 if map2 is a local permutation of map1
//!     *iflag=+2 if map2 is a global permutation of map1
//! *iflag=-1 otherwise (the maps are incommpatible)
//! if *iflag>=0, mvec_to_mvec can be used to transform vectors with the one map to vectors with the other.
//! 
//! This function may require communication, so it has to be called by all processes in the map's comm object.
//! 
void phist_maps_compatible(phist_const_map_ptr map1, phist_const_map_ptr map2, int* iflag);
//!@}

//! \name STREAM benchmarks for measuring main memory bandwidth. These are implemented in the common/
//! subdirectory for those kernel libraries that use only MPI+threads
//!@{

//! simple load dominated streaming benchmark
void phist_bench_stream_load(double* mean_bw, double *max_bw, int* iflag);
//! simple store dominated streaming benchmark
void phist_bench_stream_store(double* mean_bw, double* max_bw, int* iflag);
//! simple stream triad benchmark
void phist_bench_stream_triad(double* mean_bw, double* max_bw, int* iflag);
//!@}

//! this function should not be called by the user but by each kernel lib in kernels_init()
void phist_kernels_common_init(int *argc, char*** argv, int* iflag);

//! this function should not be called by the user but by each kernel lib in kernels_finalize()
void phist_kernels_common_finalize(int* iflag);

#ifdef __cplusplus
} // extern "C"
#endif

//! \name functions to fill mvecs and sparseMats
//!@{
#ifdef PHIST_HAVE_GHOST
typedef ghost_sparsemat_rowfunc phist_sparseMat_rowFunc;
#else
typedef int (*phist_sparseMat_rowFunc)(ghost_gidx, ghost_lidx *,
        ghost_gidx *, void *, void *);
#endif        

typedef int (*phist_mvec_elemFunc)(ghost_gidx, ghost_lidx, void *, void *);

//@}

/* this allows to use only the type-independent part of the header for 
   interface and doc generation
 */
#ifndef DOXYGEN_NO_TG


//! include file for the basic operations (kernels)
//! in single/double, real/complex.
#ifdef PHIST_HAVE_SP
//! \name single precision real kernel functions
//!@{
#include "phist_gen_s.h"
#include "phist_kernels_decl.h"
#include "phist_kernels_common_decl.h"
#include "phist_kernels_fused_decl.h"
#include "phist_kernels_carp_decl.h"
//!@}
/*!\name single precision complex kernel functions */
//!@{
#include "phist_gen_c.h"
#include "phist_kernels_decl.h"
#include "phist_kernels_common_decl.h"
#include "phist_kernels_fused_decl.h"
#include "phist_kernels_carp_decl.h"
//!@}
#endif /* PHIST_HAVE_SP */

/*!\name double precision real kernel functions */
//!@{
#include "phist_gen_d.h"
#include "phist_kernels_decl.h"
#include "phist_kernels_common_decl.h"
#include "phist_kernels_fused_decl.h"
#include "phist_kernels_carp_decl.h"
//!@}
/*!\name double precision complex kernel functions */
//!@{
#include "phist_gen_z.h"
#include "phist_kernels_decl.h"
#include "phist_kernels_common_decl.h"
#include "phist_kernels_fused_decl.h"
#include "phist_kernels_carp_decl.h"
//!@}
#include "phist_gen_clean.h"
#endif /* DOXYGEN_NO_TG */

//@}

#endif 

