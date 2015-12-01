/*! \file magma/kernels.cpp
 * wraps implementation of magma kernel routines
 * \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
 * \author "Jonas Thies <Jonas.Thies@DLR.de>
 *
*/

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "phist_macros.h"
#include "phist_kernel_perfmodels.hpp"
#include "../phist_kernels.h"

#include "phist_typedefs.h"
#include "typedefs.hpp"
#include "phist_ScalarTraits.hpp"

#ifdef PHIST_HAVE_LIKWID
#include <likwid.h>
#endif

extern "C" {


// initialize
void phist_kernels_init(int* argc, char*** argv, int* iflag)
{
  *iflag=0;

  PHIST_CHK_IERR( *iflag = MPI_Init(argc, argv), *iflag);

  PHIST_CHK_IERR( *iflag = magma_init(), *iflag);

  magma_print_environment();

  PHIST_CHK_IERR(phist_kernels_common_init(argc,argv,iflag),*iflag);
}

// finalize magma kernels
void phist_kernels_finalize(int* iflag)
{
  *iflag=0;
  PHIST_CHK_IERR( *iflag = magma_finalize(), *iflag);
  PHIST_CHK_IERR( *iflag = MPI_Finalize(), *iflag);
  PHIST_CHK_IERR(phist_kernels_common_finalize(iflag),*iflag);
}

void phist_comm_create(comm_ptr_t* vcomm, int* iflag)
{
  *iflag = 0;
  *vcomm = NULL;
}
void phist_comm_delete(comm_ptr_t vcomm, int* iflag)
{
  *iflag = 0;
}
void phist_comm_get_rank(const_comm_ptr_t vcomm, int* rank, int* iflag)
{
#ifdef PHIST_HAVE_MPI
  PHIST_CHK_IERR( *iflag = MPI_Comm_rank(MPI_COMM_WORLD, rank), *iflag);
#endif
}
void phist_comm_get_size(const_comm_ptr_t vcomm, int* size, int* iflag)
{
#ifdef PHIST_HAVE_MPI
  PHIST_CHK_IERR( *iflag = MPI_Comm_size(MPI_COMM_WORLD, size), *iflag);
#endif
}
#ifdef PHIST_HAVE_MPI
void phist_comm_get_mpi_comm(const_comm_ptr_t comm, MPI_Comm* mpiComm, int* iflag)
{
  *iflag=0;
  *mpiComm = MPI_COMM_WORLD;
}
#endif
void phist_map_create(map_ptr_t* vmap, const_comm_ptr_t vcomm, gidx_t nglob, int *iflag)
{
  *iflag = 0;
  gidx_t *n = new gidx_t;
  *n = nglob;
  *vmap = n;
}
void phist_map_delete(map_ptr_t vmap, int *iflag)
{
  *iflag = 0;
  delete (gidx_t*)vmap;
}
void phist_map_get_comm(const_map_ptr_t vmap, const_comm_ptr_t* vcomm, int* iflag)
{
  *iflag = 0;
  *vcomm = NULL;
}
void phist_map_get_local_length(const_map_ptr_t vmap, lidx_t* nloc, int* iflag)
{
  *iflag = 0;
  *nloc = (lidx_t) *((gidx_t*)vmap);
}
void phist_map_get_global_length(const_map_ptr_t vmap, gidx_t* nglob, int* iflag)
{
  *iflag = 0;
  *nglob = (lidx_t) *((gidx_t*)vmap);
}
void phist_map_get_ilower(const_map_ptr_t vmap, gidx_t* ilower, int* iflag)
{
  *iflag = 0;
  *ilower = 0;
}
void phist_map_get_iupper(const_map_ptr_t vmap, gidx_t* iupper, int* iflag)
{
  *iflag = 0;
  *iupper = (lidx_t) *((gidx_t*)vmap) - 1;
}


} //extern "C"

#include "phist_gen_s.h"
#include "kernels_def.hpp"
#include "../common/carp_noimpl.c"
#include "../common/kernels_no_inplace_VC.cpp"
#include "../common/kernels_no_fused.cpp"

#include "phist_gen_d.h"
#include "kernels_def.hpp"
#include "../common/carp_noimpl.c"
#include "../common/kernels_no_inplace_VC.cpp"
#include "../common/kernels_no_fused.cpp"

#include "phist_gen_c.h"
#include "kernels_def.hpp"
#include "../common/carp_noimpl.c"
#include "../common/kernels_no_inplace_VC.cpp"
#include "../common/kernels_no_fused.cpp"

#include "phist_gen_z.h"
#include "kernels_def.hpp"
#include "../common/carp_noimpl.c"
#include "../common/kernels_no_inplace_VC.cpp"
#include "../common/kernels_no_fused.cpp"

