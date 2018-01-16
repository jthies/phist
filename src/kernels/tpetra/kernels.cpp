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
#include "phist_typedefs.h"
#include "phist_tpetra_typedefs.hpp"
#include "../phist_kernels.h"
#include "phist_trilinos_macros.hpp"

#include "Kokkos_Core.hpp"
#include "Tpetra_Core.hpp"

#include "Kokkos_hwloc.hpp"

#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_DefaultComm.hpp"
#ifdef PHIST_HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"
#include "Tpetra_MpiPlatform.hpp"
#else
#include "Tpetra_SerialPlatform.hpp"
#endif
#include "Tpetra_Map_decl.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Tpetra_MatrixIO.hpp"

#include <fstream>
#include <sstream>
#include <cstdlib>


using namespace phist::tpetra;

namespace 
{
static int myMpiSession = 0;
} // Anonymous namespace

extern "C" void phist_kernels_init(int* argc, char*** argv, int* iflag)
{
  // Check if PHIST_NUM_THREADS is set and use it
  // Else do the same for OMP_NUM_THREADS
  // if both are not available, default to 1
  int numThreads = std::getenv("PHIST_NUM_THREADS") != nullptr ? 
                        std::strtol(std::getenv("PHIST_NUM_THREADS"), nullptr, 10) 
                      :
                        std::getenv("OMP_NUM_THREADS") != nullptr ?
                            std::strtol(std::getenv("OMP_NUM_THREADS"), nullptr, 10) 
                          : 
                            -1;  

int numNuma = -1;

#ifdef PHIST_TRY_TO_PIN_THREADS
if (Kokkos::hwloc::available())
{
  // note: by default, Kokkos will use Hyperthreads as well, we don't do that if
  // the installation gives us the freedom to pin threads.
  numNuma = Kokkos::hwloc::get_available_numa_count();
  if (numThreads==-1) numThreads=Kokkos::hwloc::get_available_cores_per_numa()*numNuma;
}
#endif

#if TRILINOS_MAJOR_MINOR_VERSION>=121300
  Kokkos::InitArguments args{numThreads, numNuma};
#else
  Kokkos::InitArguments args; args.num_threads=numThreads; args.num_numa=numNuma;
#endif
  PHIST_TRY_CATCH(Kokkos::initialize(args), *iflag);
  MPI_Init(argc, argv);

#ifdef PHIST_TRY_TO_PIN_THREADS
if (Kokkos::hwloc::available() && Kokkos::hwloc::can_bind_threads())
{
  PHIST_OUT(PHIST_WARNING,"Tpetra/Kokkos/OpenMP: pinning of threads is not implemented!\n"
                          "                      Suggesting %d threads, %d NUMA domains on this MPI rank.\n",
                          numThreads,numNuma);
  //bool bind_this_thread( const std::pair<unsigned,unsigned> );
}
#endif

  PHIST_CHK_IERR(phist_kernels_common_init(argc, argv, iflag), *iflag);

  *iflag = PHIST_SUCCESS;
}
            
extern "C" void phist_kernels_finalize(int* iflag)
{
  PHIST_CHK_IERR(phist_kernels_common_finalize(iflag),*iflag);
  // Does not throw
  Tpetra::finalize();
}

//!
extern "C" void phist_comm_create(phist_comm_ptr* vcomm, int* iflag)
{
  *vcomm = (phist_comm_ptr)Tpetra::getDefaultComm().get();  
  *iflag = PHIST_SUCCESS;
}

//!
extern "C" void phist_comm_delete(phist_comm_ptr vcomm, int* iflag)
{
  *iflag = PHIST_SUCCESS;
}

#ifdef PHIST_HAVE_MPI
void phist_comm_get_mpi_comm(phist_const_comm_ptr vcomm, MPI_Comm* mpiComm, int* iflag)
{
  PHIST_CAST_PTR_FROM_VOID(const comm_type, comm, vcomm, *iflag);
  const auto teuchosMpiComm = dynamic_cast<const Teuchos::MpiComm<int>* >(comm);
  if (teuchosMpiComm == nullptr) // invalid input, return
  {
    *mpiComm = MPI_COMM_SELF;
    *iflag = PHIST_INVALID_INPUT;
    return;
  }

  auto wrappedMpiComm = teuchosMpiComm->getRawMpiComm();

  if (wrappedMpiComm != Teuchos::null )
  {
    *mpiComm = (*wrappedMpiComm)();
    return;
  }

  *mpiComm = MPI_COMM_SELF;
  *iflag = PHIST_INVALID_INPUT;
}
#endif

//!
extern "C" void phist_comm_get_rank(phist_const_comm_ptr vcomm, int* rank, int* iflag)
{
  PHIST_CAST_PTR_FROM_VOID(const comm_type, comm, vcomm, *iflag);
  *rank = comm->getRank();
  *iflag = PHIST_SUCCESS;
}
//!
extern "C" void phist_comm_get_size(phist_const_comm_ptr vcomm, int* size, int* iflag)
{
  PHIST_CAST_PTR_FROM_VOID(const comm_type, comm, vcomm, *iflag);
  *size = comm->getSize();
  *iflag = PHIST_SUCCESS;
}

//!
extern "C" void phist_map_create(phist_map_ptr* vmap, phist_const_comm_ptr vcomm, phist_gidx nglob, int *iflag)
{
  PHIST_CAST_PTR_FROM_VOID(const comm_type, comm, vcomm, *iflag);
  auto phistCommPtr = Teuchos::rcp(comm, false);
  auto mapPtr = Tpetra::createUniformContigMap<phist_lidx, phist_gidx>(nglob, phistCommPtr);
  *vmap = (phist_map_ptr)(mapPtr.release().get());
  *iflag = PHIST_SUCCESS;
}

//!
extern "C" void phist_map_delete(phist_map_ptr vmap, int *iflag)
{
  PHIST_CAST_PTR_FROM_VOID(const map_type, map, vmap, *iflag);
  delete map;
  vmap = nullptr;
  *iflag = PHIST_SUCCESS;
}

//!
extern "C" void phist_map_get_comm(phist_const_map_ptr vmap, phist_const_comm_ptr* vcomm, int* iflag)
{
  PHIST_CAST_PTR_FROM_VOID(const map_type, map, vmap, *iflag);
  auto comm = map->getComm();
  *vcomm = (phist_const_comm_ptr)(comm.get());
  *iflag = PHIST_SUCCESS;
}

//!
extern "C" void phist_map_get_local_length(phist_const_map_ptr vmap, phist_lidx* nloc, int* iflag)
{
  PHIST_CAST_PTR_FROM_VOID(const map_type, map, vmap, *iflag);
  *nloc = map->getNodeNumElements();
  *iflag = PHIST_SUCCESS;
}

//!
extern "C" void phist_map_get_global_length(phist_const_map_ptr vmap, phist_gidx* nglob, int* iflag)
{
  PHIST_CAST_PTR_FROM_VOID(const map_type, map, vmap, *iflag);
  *nglob = map->getGlobalNumElements();
  *iflag = PHIST_SUCCESS;
}

//! returns the smallest global index in the map appearing on my partition. iflag is set to 1
//! in case the map is not contiguous, because in that case it may be that the
//! caller falsely assumes global elements [ilower ... iupper] are actually on this partition.
extern "C" void phist_map_get_ilower(phist_const_map_ptr vmap, phist_gidx* ilower, int* iflag)
{
  PHIST_CAST_PTR_FROM_VOID(const map_type, map, vmap, *iflag);

  *iflag = map->isContiguous() ? PHIST_SUCCESS : PHIST_WARNING;
  *ilower = map->getMinGlobalIndex();
}
//! returns the largest global index in the map appearing on my partition. iflag is set to 1
//! in case the map is not contiguous, because in that case it may be that the
//! caller falsely assumes global elements [ilower ... iupper] are actually on this partition.
extern "C" void phist_map_get_iupper(phist_const_map_ptr vmap, phist_gidx* iupper, int* iflag)
{
  PHIST_CAST_PTR_FROM_VOID(const map_type, map, vmap, *iflag);

  *iflag = map->isContiguous() ? PHIST_SUCCESS : PHIST_WARNING;
  *iupper = map->getMaxGlobalIndex();
}

// allow verifying that maps are compatible (the same or permuted)
extern "C" void phist_maps_compatible(phist_const_map_ptr vmap1, phist_const_map_ptr vmap2, int* iflag)
{
  PHIST_CAST_PTR_FROM_VOID(const map_type, map1, vmap1, *iflag);
  PHIST_CAST_PTR_FROM_VOID(const map_type, map2, vmap2, *iflag);

  // Same object or same map
  if (map1 == map2 || *map1 == *map2)
  {
    *iflag = PHIST_SUCCESS;
    return;
  }

  auto comm1 = map1->getComm();
  auto comm2 = map2->getComm();

  if (comm1 != comm2)
  {
    #ifdef PHIST_HAVE_MPI
      MPI_Comm mpiComm1;
      PHIST_CHK_NEG_IERR(phist_comm_get_mpi_comm(comm1.get(), &mpiComm1, iflag), *iflag);

      MPI_Comm mpiComm2;
      PHIST_CHK_NEG_IERR(phist_comm_get_mpi_comm(comm2.get(), &mpiComm2, iflag), *iflag);

      if (mpiComm1 != mpiComm2)
      {
        *iflag = PHIST_FUNCTIONAL_ERROR;
        return;
      }
    #else
      *iflag = PHIST_NOT_IMPLEMENTED;
      return;
    #endif
  }

    bool isCompat = map1->isCompatible(*map2);
    bool isLocallySame = map1->locallySameAs(*map2);

    #ifdef PHIST_HAVE_MPI
      MPI_Comm mpiComm;
      PHIST_CHK_IERR(phist_comm_get_mpi_comm(comm1.get(), &mpiComm, iflag), *iflag);
      PHIST_CHK_MPIERR(*iflag = MPI_Allreduce(MPI_IN_PLACE, &isLocallySame, 
                                              1, MPI_LOGICAL, MPI_LOR, mpiComm),
                                              *iflag);
    #endif
    
    if (isCompat && !isLocallySame)
    {
      *iflag = PHIST_WARNING;
      return;
    }

    if (isCompat)
    {
      *iflag = PHIST_INFO;
      return;
    }

    // TODO: Look up how Frank does this nicely
    if (
              map1->getIndexBase() == map2->getIndexBase() 
          &&
              map1->getGlobalNumElements() == map2->getGlobalNumElements()
          &&
              map1->getMaxGlobalIndex() == map2->getMaxGlobalIndex()
          &&
              map1->getMinGlobalIndex() == map2->getMinGlobalIndex()
       )
    {
      *iflag = PHIST_INFO;
      return;
     }

  *iflag = PHIST_FUNCTIONAL_ERROR;
}


#include "../common/default_context.cpp"
#include "../common/phist_bench_kernels.cpp"
