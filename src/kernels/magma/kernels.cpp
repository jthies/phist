/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
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

void phist_comm_create(phist_comm_ptr* vcomm, int* iflag)
{
  *iflag = 0;
  *vcomm = NULL;
}
void phist_comm_delete(phist_comm_ptr vcomm, int* iflag)
{
  *iflag = 0;
}
void phist_comm_get_rank(phist_const_comm_ptr vcomm, int* rank, int* iflag)
{
#ifdef PHIST_HAVE_MPI
  PHIST_CHK_IERR( *iflag = MPI_Comm_rank(MPI_COMM_WORLD, rank), *iflag);
#endif
}
void phist_comm_get_size(phist_const_comm_ptr vcomm, int* size, int* iflag)
{
#ifdef PHIST_HAVE_MPI
  PHIST_CHK_IERR( *iflag = MPI_Comm_size(MPI_COMM_WORLD, size), *iflag);
#endif
}
#ifdef PHIST_HAVE_MPI
void phist_comm_get_mpi_comm(phist_const_comm_ptr comm, MPI_Comm* mpiComm, int* iflag)
{
  *iflag=0;
  *mpiComm = MPI_COMM_WORLD;
}
#endif
void phist_map_create(map_ptr* vmap, phist_const_comm_ptr vcomm, phist_gidx nglob, int *iflag)
{
  *iflag = 0;
  phist_gidx *n = new phist_gidx;
  *n = nglob;
  *vmap = n;
}
void phist_map_delete(map_ptr vmap, int *iflag)
{
  *iflag = 0;
  delete (phist_gidx*)vmap;
}
void phist_map_get_comm(phist_const_map_ptr vmap, phist_const_comm_ptr* vcomm, int* iflag)
{
  *iflag = 0;
  *vcomm = NULL;
}
void phist_map_get_local_length(phist_const_map_ptr vmap, phist_lidx* nloc, int* iflag)
{
  *iflag = 0;
  *nloc = (phist_lidx) *((phist_gidx*)vmap);
}
void phist_map_get_global_length(phist_const_map_ptr vmap, phist_gidx* nglob, int* iflag)
{
  *iflag = 0;
  *nglob = (phist_lidx) *((phist_gidx*)vmap);
}
void phist_map_get_ilower(phist_const_map_ptr vmap, phist_gidx* ilower, int* iflag)
{
  *iflag = 0;
  *ilower = 0;
}
void phist_map_get_iupper(phist_const_map_ptr vmap, phist_gidx* iupper, int* iflag)
{
  *iflag = 0;
  *iupper = (phist_lidx) *((phist_gidx*)vmap) - 1;
}


} //extern "C"

#include "../common/phist_bench_kernels.cpp"

#include "phist_gen_s.h"
#include "kernels_def.hpp"
#include "../common/default_mvec_get_data_def.hpp"
#include "kernels_no_carp.cpp"
#include "kernels_no_inplace_VC.cpp"
#include "kernels_no_fused.cpp"

#include "phist_gen_d.h"
#include "kernels_def.hpp"
#include "../common/default_mvec_get_data_def.hpp"
#include "kernels_no_carp.cpp"
#include "kernels_no_inplace_VC.cpp"
#include "kernels_no_fused.cpp"

#include "phist_gen_c.h"
#include "kernels_def.hpp"
#include "../common/default_mvec_get_data_def.hpp"
#include "kernels_no_carp.cpp"
#include "kernels_no_inplace_VC.cpp"
#include "kernels_no_fused.cpp"

#include "phist_gen_z.h"
#include "kernels_def.hpp"
#include "../common/default_mvec_get_data_def.hpp"
#include "kernels_no_carp.cpp"
#include "kernels_no_inplace_VC.cpp"
#include "../common/kernels_no_VC_add_WD.cpp"
#include "kernels_no_fused.cpp"

