/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
/*! \file petsc/kernels.cpp
 * wraps implementation of petsc kernel routines
 * \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
 * \author "Jonas Thies <j.thies@tudelft.nl>
 *
*/

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <petscsys.h>

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
  PHIST_CHK_IERR( *iflag = PetscInitialize(argc,argv,NULL,NULL),*iflag);
  PHIST_CHK_IERR(phist_kernels_common_init(argc,argv,iflag),*iflag);
}

// finalize kernels
void phist_kernels_finalize(int* iflag)
{
  PHIST_CHK_IERR(phist_kernels_common_finalize(iflag),*iflag);
  PHIST_CHK_IERR( *iflag = PetscFinalize(), *iflag);
}

void phist_comm_create(phist_comm_ptr* vcomm, int* iflag)
{
  MPI_Comm *comm = NULL;
  PHIST_CHK_IERR( *iflag = PetscNew(&comm), *iflag);
  PHIST_CHK_IERR( *iflag = PetscCommDuplicate(phist_get_default_comm(), comm, NULL), *iflag);
  *vcomm = comm;
}
void phist_comm_delete(phist_comm_ptr vcomm, int* iflag)
{
  MPI_Comm *comm = (MPI_Comm*)vcomm;
  PHIST_CHK_IERR( *iflag = PetscCommDestroy(comm), *iflag);
  PHIST_CHK_IERR( *iflag = PetscFree(comm), *iflag);
}
void phist_comm_get_rank(phist_const_comm_ptr vcomm, int* rank, int* iflag)
{
#ifdef PHIST_HAVE_MPI
  PHIST_CHK_IERR( *iflag = MPI_Comm_rank(*((const MPI_Comm*)vcomm), rank), *iflag);
#endif
}
void phist_comm_get_size(phist_const_comm_ptr vcomm, int* size, int* iflag)
{
#ifdef PHIST_HAVE_MPI
  PHIST_CHK_IERR( *iflag = MPI_Comm_size(*((const MPI_Comm*)vcomm), size), *iflag);
#endif
}
#ifdef PHIST_HAVE_MPI
void phist_comm_get_mpi_comm(phist_const_comm_ptr comm, MPI_Comm* mpiComm, int* iflag)
{
  *iflag=0;
  *mpiComm = *(MPI_Comm*)comm;
}
#endif
namespace
{
  struct MapData
  {
    phist_gidx nglob, ilower;
    phist_lidx nlocal;
    phist_const_comm_ptr vcomm;
  };
  bool operator==(const MapData& m1, const MapData& m2)
  {
    if( m1.nglob != m2.nglob ) return false;
    if( m1.nlocal != m2.nlocal ) return false;
    if( m1.ilower != m2.ilower ) return false;
    // TODO communicator comparison ...
    return true;
  }
}
void phist_map_create(phist_map_ptr* vmap, phist_const_comm_ptr vcomm, phist_gidx nglob, int *iflag)
{
  int irank, nprocs;
  PHIST_CHK_IERR( phist_comm_get_rank(vcomm,&irank,iflag), *iflag);
  PHIST_CHK_IERR( phist_comm_get_size(vcomm,&nprocs,iflag), *iflag);

  MapData *map = NULL;
  PHIST_CHK_IERR( *iflag = PetscNew(&map), *iflag);
  map->nglob = nglob;
  map->vcomm = vcomm;
  map->nlocal = nglob/nprocs;
  map->ilower = irank*map->nlocal;
  if( irank < nglob % nprocs )
  {
    map->nlocal++;
    map->ilower += irank;
  }
  else
  {
    map->ilower += nglob % nprocs;
  }
  *vmap = map;
}
void phist_map_delete(phist_map_ptr vmap, int *iflag)
{
  PHIST_CHK_IERR( *iflag = PetscFree(vmap), *iflag);
}
void phist_map_get_comm(phist_const_map_ptr vmap, phist_const_comm_ptr* vcomm, int* iflag)
{
  *iflag = 0;
  const MapData *map = (const MapData*)vmap;
  *vcomm = map->vcomm;
}
void phist_map_get_local_length(phist_const_map_ptr vmap, phist_lidx* nloc, int* iflag)
{
  *iflag = 0;
  const MapData *map = (const MapData*)vmap;
  *nloc = map->nlocal;
}
void phist_map_get_global_length(phist_const_map_ptr vmap, phist_gidx* nglob, int* iflag)
{
  *iflag = 0;
  const MapData *map = (const MapData*)vmap;
  *nglob = map->nglob;
}
void phist_map_get_ilower(phist_const_map_ptr vmap, phist_gidx* ilower, int* iflag)
{
  *iflag = 0;
  const MapData *map = (const MapData*)vmap;
  *ilower = map->ilower;
}
void phist_map_get_iupper(phist_const_map_ptr vmap, phist_gidx* iupper, int* iflag)
{
  *iflag = 0;
  const MapData *map = (const MapData*)vmap;
  *iupper = map->ilower + map->nlocal - 1;
}
void phist_maps_compatible(phist_const_map_ptr vmap1, phist_const_map_ptr vmap2, int* iflag)
{
  const MapData *map1 = (const MapData*)vmap1;
  const MapData *map2 = (const MapData*)vmap2;
  if( *map1 == *map2 )
    *iflag = 0;
  else
    *iflag = -1;
}


} //extern "C"

#include "../common/phist_bench_kernels.cpp"
#include "../common/default_context.cpp"
