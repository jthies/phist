#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "phist_macros.h"
#include "phist_typedefs.h"
#include "phist_tpetra_typedefs.hpp"
#include "../phist_kernels.h"
#include "phist_trilinos_macros.h"
#include "phist_ScalarTraits.hpp"

#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_DefaultComm.hpp"
#ifdef PHIST_HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"
#endif
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Tpetra_MatrixIO.hpp"
#include "Tpetra_DefaultPlatform.hpp"

#ifdef PHIST_HAVE_LIKWID
#include <likwid.h>
#endif

#include <fstream>
#include <sstream>


using namespace phist::tpetra;

namespace {
static int myMpiSession=0;
}

// initialize kernel library. Should at least call MPI_Init if it has not been called
// but is required.
extern "C" void phist_kernels_init(int* argc, char*** argv, int* iflag)
{
  *iflag=0;
#ifdef PHIST_HAVE_MPI
  int mpi_initialized;
  MPI_Initialized(&mpi_initialized);
  myMpiSession= mpi_initialized? 0: 1;
  if (myMpiSession==1)
  {
    *iflag=MPI_Init(argc,argv);
  }
#endif
  // describe yourself!
  std::ostringstream oss;
  Tpetra::DefaultPlatform::getDefaultPlatform().describe(oss,Teuchos::VERB_MEDIUM);
  PHIST_SOUT(PHIST_INFO,"Tpetra platform:\n%s\n", oss.str().c_str());
  PHIST_CHK_IERR(phist_kernels_common_init(argc,argv,iflag),*iflag);
}
      
  // finalize kernel library. Should at least call MPI_Finalize if it has not been called
  // but is required.
  extern "C" void phist_kernels_finalize(int* iflag)
  {
    *iflag=0;
    PHIST_CHK_IERR(phist_kernels_common_finalize(iflag),*iflag);
#ifdef PHIST_HAVE_MPI
  if (myMpiSession==1)
  {
    *iflag=MPI_Finalize();
  }
#endif
}
            

//!
extern "C" void phist_comm_create(phist_comm_ptr* vcomm, int* iflag)
{
  *iflag=0;
  *vcomm = (phist_comm_ptr)(Teuchos::DefaultComm<int>::getComm().get());
}

//!
extern "C" void phist_comm_delete(phist_comm_ptr vcomm, int* iflag)
{
  *iflag=0;
  PHIST_TOUCH(vcomm);
  // note - as comm_create returns a raw pointer to the default comm, we should not delete it.
}

#ifdef PHIST_HAVE_MPI
void phist_comm_get_mpi_comm(phist_const_comm_ptr vcomm, MPI_Comm* mpiComm, int* iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const comm_type,comm,vcomm,*iflag);
  const Teuchos::MpiComm<int>* Teuchos_mpiComm=dynamic_cast<const Teuchos::MpiComm<int>*>(comm);
  if (Teuchos_mpiComm==NULL)
  {
    *mpiComm = MPI_COMM_SELF;
    *iflag=PHIST_INVALID_INPUT;
    return;
  }
  Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > wrapped_mpiComm
        = Teuchos_mpiComm->getRawMpiComm();
  if (wrapped_mpiComm==Teuchos::null)
  {
    *mpiComm = MPI_COMM_SELF;
    *iflag=PHIST_INVALID_INPUT;
  }
  else
  {
    *mpiComm = (*wrapped_mpiComm)();
  }
}
#endif

//!
extern "C" void phist_comm_get_rank(phist_const_comm_ptr vcomm, int* rank, int* iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const comm_type,comm,vcomm,*iflag);
  *rank=comm->getRank();
}
//!
extern "C" void phist_comm_get_size(phist_const_comm_ptr vcomm, int* size, int* iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const comm_type,comm,vcomm,*iflag);
  *size=comm->getSize();
}

// This function creates the default compute node to be used for vectors and matrices.
// It is not exposed to the user because the Kokkos node concept is specific for Tpetra
// and doesn't have a direct equivalent in epetra or ghost. The function checks if a file
// node.xml exists in which the node parameters like "Num Threads" can be set, otherwise
// it just uses default parameters.
extern "C" void phist_tpetra_node_create(node_type** node, phist_const_comm_ptr vcomm, int* iflag)
{
  // print messages only once
  int outputLevel;
  {
    static int globalMessagesPrinted = 0;
    int messagesPrinted;
#pragma omp atomic capture
    messagesPrinted = globalMessagesPrinted++;
    outputLevel = messagesPrinted > 0 ? PHIST_DEBUG : PHIST_INFO;
  }

  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const comm_type,comm,vcomm,*iflag);
  Teuchos::RCP<Teuchos::ParameterList> nodeParams=Teuchos::rcp(new Teuchos::ParameterList);
  // check if the file exists
  bool haveNodeFile=false;
  if (comm->getRank()==0)
    {
    std::ifstream test("phist_node.xml");
    if (test) haveNodeFile=true;
    }
  
  comm->broadcast(0,sizeof(bool),(char*)&haveNodeFile);

  if (haveNodeFile)
  {
    bool status=true;
    try 
    {
      Teuchos::updateParametersFromXmlFileAndBroadcast("phist_node.xml",nodeParams.ptr(),*comm);
    } TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,status);
    if (status==false)
      {
      *iflag=PHIST_CAUGHT_EXCEPTION; 
      return;
      }
  }
  else
  {
    PHIST_SOUT(outputLevel,"File phist_node.xml not found, using default node settings in tpetra\n");
  }
  const char* PHIST_NUM_THREADS=getenv("PHIST_NUM_THREADS");
  if (PHIST_NUM_THREADS!=NULL)
  {
    int nThreads=atoi(PHIST_NUM_THREADS);
    if (nThreads>0)
    {
      PHIST_SOUT(outputLevel,"taking #threads from env variable PHIST_NUM_THREADS\n");
      nodeParams->set("Num Threads",nThreads);
    }
  }
  PHIST_SOUT(outputLevel,"# threads requested: %d\n",nodeParams->get("Num Threads",0));
  *node = new node_type(*nodeParams);
}

//!
extern "C" void phist_map_create(phist_map_ptr* vmap, phist_const_comm_ptr vcomm, phist_gidx nglob, int *iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const comm_type,comm,vcomm,*iflag);
  
  Teuchos::RCP<const comm_type> phist_comm_ptr = Teuchos::rcp(comm,false);
node_type* node;
  PHIST_CHK_IERR(phist_tpetra_node_create(&node,vcomm,iflag),*iflag);
  Teuchos::RCP<node_type> node_ptr = Teuchos::rcp(node,true);
  Teuchos::RCP<const map_type> map_ptr =
        Tpetra::createUniformContigMapWithNode<phist_lidx,phist_gidx,node_type>
                (nglob, phist_comm_ptr,node_ptr);
  *vmap = (phist_map_ptr)(map_ptr.release().get());
}

//!
extern "C" void phist_map_delete(phist_map_ptr vmap, int *iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(map_type,map,vmap,*iflag);
  delete map;
  vmap=NULL;
}

//!
extern "C" void phist_map_get_comm(phist_const_map_ptr vmap, phist_const_comm_ptr* vcomm, int* iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const map_type,map,vmap,*iflag);
  Teuchos::RCP<const comm_type> comm = map->getComm();
  *vcomm = (phist_const_comm_ptr)(comm.get());
}

//!
extern "C" void phist_map_get_local_length(phist_const_map_ptr vmap, phist_lidx* nloc, int* iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const map_type,map,vmap,*iflag);
  *nloc = map->getNodeNumElements();
}

//!
extern "C" void phist_map_get_global_length(phist_const_map_ptr vmap, phist_gidx* nglob, int* iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const map_type,map,vmap,*iflag);
  *nglob = map->getGlobalNumElements();
}

//! returns the smallest global index in the map appearing on my partition. iflag is set to 1
//! in case the map is not contiguous, because in that case it may be that the
//! caller falsely assumes global elements [ilower ... iupper] are actually on this partition.
extern "C" void phist_map_get_ilower(phist_const_map_ptr vmap, phist_gidx* ilower, int* iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const map_type,map,vmap,*iflag);
  if (map->isContiguous()==false) *iflag=1;
  *ilower = map->getMinGlobalIndex();
}
//! returns the largest global index in the map appearing on my partition. iflag is set to 1
//! in case the map is not contiguous, because in that case it may be that the
//! caller falsely assumes global elements [ilower ... iupper] are actually on this partition.
extern "C" void phist_map_get_iupper(phist_const_map_ptr vmap, phist_gidx* iupper, int* iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const map_type,map,vmap,*iflag);
  if (map->isContiguous()==false) *iflag=1;
  *iupper = map->getMaxGlobalIndex();
}

// allow verifying that maps are compatible (the same or permuted)
extern "C" void phist_maps_compatible(phist_const_map_ptr vmap1, phist_const_map_ptr vmap2, int* iflag)
{
  PHIST_CAST_PTR_FROM_VOID(const map_type,map1,vmap1,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const map_type,map2,vmap2,*iflag);
  *iflag=-1;
  // same object?
  if (map1==map2) {*iflag=0; return;}
  // maps identical?
  if (*map1==*map2) {*iflag=0; return;}
  // check wether the maps have the same communicator, some of the functions below
  // need the same number of processes etc. We simply demand the same communicator
  // and return -1 otherwise.
  Teuchos::RCP<const comm_type> comm1=map1->getComm();
  Teuchos::RCP<const comm_type> comm2=map2->getComm();
  if (comm1!=comm2) 
  {
#ifdef PHIST_HAVE_MPI
  MPI_Comm mpi_comm1, mpi_comm2;
  PHIST_CHK_NEG_IERR(phist_comm_get_mpi_comm(comm1.get(),&mpi_comm1,iflag),*iflag);
  PHIST_CHK_NEG_IERR(phist_comm_get_mpi_comm(comm2.get(),&mpi_comm2,iflag),*iflag);
  if (mpi_comm1!=mpi_comm2)
  {
    return;
  }
#else
    // can't check communicator compatibility for now
    PHIST_CHK_IERR(*iflag=-99,*iflag);
    return;
#endif
  }

  // query how compatible the maps are
  
  // same global size and distribution?
  
  bool compat=map1->isCompatible(*map2);
  bool locallySame=map1->locallySameAs(*map2);

#ifdef PHIST_HAVE_MPI  
  MPI_Comm mpi_comm;
  PHIST_CHK_IERR(phist_comm_get_mpi_comm(comm1.get(),&mpi_comm,iflag),*iflag);
  PHIST_CHK_MPIERR(*iflag=MPI_Allreduce(MPI_IN_PLACE,&locallySame,1,MPI_LOGICAL,MPI_LOR,mpi_comm),*iflag);
#endif

  if (compat&&!locallySame)
  {
    *iflag=1; // just a local permutation
    return;
  }
  
  // "compatible" in Tpetra is a strong word, it means that the partition sizes can't change,
  // so e.g. load balancing is not included. If the maps are compatible, return 2 here (global perm), 
  // otherwise do some more checks.
  if (compat) {*iflag=2; return;}
  if (map1->getIndexBase()         == map2->getIndexBase()         &&
      map1->getGlobalNumElements() == map2->getGlobalNumElements() &&
      map1->getMaxGlobalIndex() == map2->getMaxGlobalIndex() &&
      map1->getMinGlobalIndex() == map2->getMinGlobalIndex() )
  {
    *iflag=2; // assume that the maps are a global permutation (repartitioning+reordering) of each other.
              // the check is not 100% fool proof, you could have map1=[0 2 | 4 6] and map2=[0 1 | 3 6] and
              // falsely return +2 here, but then subsequent Tpetra calls will complain. Anyway, mvec_to_mvec
              // will probably still work, and we don't give the user any interface to construct such funny maps,
              // so if he does it's in a way his own fault.
    return;
  }
  return;
}


#include "../common/phist_bench_kernels.cpp"
