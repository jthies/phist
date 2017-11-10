/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_config.h"
#include "phist_tools.h"
#include "../phist_kernels.h"
#include "phist_trilinos_macros.hpp"
#include "Epetra_config.h"
#ifdef PHIST_HAVE_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_SerialComm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Import.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"

#include "Teuchos_StandardCatchMacros.hpp"
#include "EpetraExt_CrsMatrixIn.h"


#ifdef PHIST_HAVE_BELOS
#include "Belos_config.h"
# ifdef HAVE_BELOS_TSQR
# include "Epetra_TsqrAdaptor.hpp"
# include "./BelosEpetraAdapter.hpp"
# include "BelosTsqrOrthoManager.hpp"
# endif
#endif


#include "epetra_helpers.cpp"

#include "hg_repart.hpp"

namespace {
  int mpiInitializedBefore = 0;
}

// initialize kernel library. Should at least call MPI_Init if it has not been called
// but is required.
extern "C" void phist_kernels_init(int* argc, char*** argv, int* iflag)
{
  *iflag=0;
#ifdef PHIST_HAVE_MPI
  PHIST_CHK_IERR( *iflag = MPI_Initialized(&mpiInitializedBefore), *iflag);
  if (!mpiInitializedBefore)
  {
    *iflag=MPI_Init(argc,argv);
  }
#endif
  phist_kernels_common_init(argc,argv,iflag);
}

// finalize kernel library. Should at least call MPI_Finalize if it has not been called
// but is required.
extern "C" void phist_kernels_finalize(int* iflag)
{
    phist_kernels_common_finalize(iflag);
#ifdef PHIST_HAVE_MPI
  if (!mpiInitializedBefore)
{
  *iflag=MPI_Finalize();
}
#endif  
}


//!
extern "C" void phist_comm_create(phist_comm_ptr* vcomm, int* iflag)
{
  *iflag=0;
#ifdef PHIST_HAVE_MPI
  Epetra_Comm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_Comm *comm = new Epetra_SerialComm;
#endif
  if (comm==NULL) *iflag=-1;
  *vcomm=(phist_comm_ptr)(comm);
}

//!
extern "C" void phist_comm_delete(phist_comm_ptr vcomm, int* iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(Epetra_Comm,comm,vcomm,*iflag);
  delete comm;
  vcomm=NULL;
}

#ifdef PHIST_HAVE_MPI
extern void phist_comm_get_mpi_comm(phist_const_comm_ptr vcomm, MPI_Comm* mpiComm, int* iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_MpiComm,comm,vcomm,*iflag);
  *mpiComm=comm->Comm();
}
#endif

//!
extern "C" void phist_comm_get_rank(phist_const_comm_ptr vcomm, int* rank, int* iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_Comm,comm,vcomm,*iflag);
  *rank=comm->MyPID();
}
//!
extern "C" void phist_comm_get_size(phist_const_comm_ptr vcomm, int* size, int* iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_Comm,comm,vcomm,*iflag);
  *size=comm->NumProc();
}
//!
extern "C" void phist_map_create(phist_map_ptr* vmap, phist_const_comm_ptr vcomm, phist_gidx nglob, int *iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_Comm,comm,vcomm,*iflag);
  Epetra_BlockMap* map;
  PHIST_TRY_CATCH(map = new Epetra_Map(nglob,0,*comm),*iflag);
  *vmap=(phist_map_ptr)(map);
}

//!
extern "C" void phist_map_delete(phist_map_ptr vmap, int *iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(Epetra_BlockMap,map,vmap,*iflag);
  delete map;
  vmap=NULL;
}
  
//!
extern "C" void phist_map_get_comm(phist_const_map_ptr vmap, phist_const_comm_ptr* vcomm, int* iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_BlockMap,map,vmap,*iflag);
  *vcomm = (phist_const_comm_ptr)(&map->Comm());  
}

//!
extern "C" void phist_map_get_local_length(phist_const_map_ptr vmap, phist_lidx* nloc, int* iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_BlockMap,map,vmap,*iflag);
  *nloc = map->NumMyElements();
}

//!
extern "C" void phist_map_get_global_length(phist_const_map_ptr vmap, phist_gidx* nglob, int* iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_BlockMap,map,vmap,*iflag);
#ifdef EPETRA_NO_64BIT_GLOBAL_INDICES
  *nglob = map->NumGlobalElements();
#else
  *nglob = map->NumGlobalElements64();
#endif
}

//! returns the smallest global index in the map appearing on my partition. iflag is set to 1
//! in case the map is not contiguous, because in that case it may be that the
//! caller falsely assumes global elements [ilower ... iupper] are actually on this partition.
extern "C" void phist_map_get_ilower(phist_const_map_ptr vmap, phist_gidx* ilower, int* iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_BlockMap,map,vmap,*iflag);
  if (map->LinearMap()==false) *iflag=1;
#ifdef EPETRA_NO_64BIT_GLOBAL_INDICES
  *ilower = map->MinMyGID();
#else
  *ilower = map->MinMyGID64();
#endif
}
//! returns the largest global index in the map appearing on my partition. iflag is set to 1
//! in case the map is not contiguous, because in that case it may be that the
//! caller falsely assumes global elements [ilower ... iupper] are actually on this partition.
extern "C" void phist_map_get_iupper(phist_const_map_ptr vmap, phist_gidx* iupper, int* iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_BlockMap,map,vmap,*iflag);
  if (map->LinearMap()==false) *iflag=1;
#ifdef EPETRA_NO_64BIT_GLOBAL_INDICES
  *iupper = map->MaxMyGID();
#else
  *iupper = map->MaxMyGID64();
#endif
}

extern "C" void phist_maps_compatible(phist_const_map_ptr vmap1, phist_const_map_ptr vmap2, int* iflag)
{
  PHIST_CAST_PTR_FROM_VOID(const Epetra_BlockMap,map1,vmap1,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Epetra_BlockMap,map2,vmap2,*iflag);
  *iflag=-1;
  // same object?
  if (map1==map2) {*iflag=0; return;}
  // maps identical?
  if (map1->SameAs(*map2)) {*iflag=0; return;}

  // check the maps use the same comm
#ifdef PHIST_HAVE_MPI
  MPI_Comm mpi_comm1,mpi_comm2;
  *iflag=0;
  PHIST_CHK_IERR(phist_comm_get_mpi_comm(&map1->Comm(),&mpi_comm1,iflag),*iflag);
  PHIST_CHK_IERR(phist_comm_get_mpi_comm(&map1->Comm(),&mpi_comm2,iflag),*iflag);
  *iflag=-1;
  if (mpi_comm1!=mpi_comm2) {return;}
#endif  

  // same global size and local size? We can't tell if there is a global
  // permutation going on (actually we could, but it involves creating the
  // import object, which is too much effort for this query).
  if (map1->PointSameAs(*map2)) {*iflag=2; return;}
  // if the MPI communicator is the same, Epetra can always transform vectors from
  // one map to the other using an appropriate Import or Export object.

  // check some simple heuristics, index base, global length, global min and max index.
  // mvec_to_mvec will take care of most conversions with Epetra anyway.
  if ( map1->NumGlobalElements() == map2->NumGlobalElements() &&
       map1->IndexBase()    == map2->IndexBase()    &&
       map1->MinAllGID()    == map2->MinAllGID()    && 
       map1->MaxAllGID()    == map2->MaxAllGID() )
  {
    *iflag=2;
  }

  // if not, return the -1 (incompatible maps)
  *iflag=-1;
  return;
}

#include "../common/default_context.cpp"

#include "../common/phist_bench_kernels.cpp"

#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "../common/kernels_no_impl.c"

#include "phist_gen_c.h"
#include "../common/kernels_no_impl.cpp"
#endif
#include "phist_gen_z.h"
#include "../common/kernels_no_impl.cpp"

#include "phist_gen_d.h"
#include "kernels_def.hpp"
#include "../common/default_context_def.hpp"
#include "carp_def.hpp"
#include "../common/kernels_no_io.cpp"
#include "../common/kernels_no_gpu.cpp"
#include "../common/kernels_no_fused.cpp"
#include "../common/kernels_no_inplace_VC.cpp"
#include "../common/kernels_no_VC_add_WD.cpp"
