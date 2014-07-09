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
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Tpetra_MatrixIO.hpp"

#include "BelosTpetraAdapter.hpp"
#include "BelosTsqrOrthoManager.hpp"

#ifdef PHIST_TIMEMONITOR
#include "phist_timemonitor.hpp"
namespace phist_TimeMonitor
{
  Timer::TimeDataMap Timer::_timingResults;
}
#endif

#ifdef PHIST_HAVE_LIKWID
#include <likwid.h>
#endif

#include <fstream>


using namespace phist::tpetra;


static int myMpiSession=0;

// initialize kernel library. Should at least call MPI_Init if it has not been called
// but is required.
extern "C" void phist_kernels_init(int* argc, char*** argv, int* ierr)
{
#ifdef PHIST_HAVE_MPI
  int mpi_initialized;
  MPI_Initialized(&mpi_initialized);
  myMpiSession= mpi_initialized? 0: 1;
  if (myMpiSession==1)
  {
    *ierr=MPI_Init(argc,argv);
  }
#endif
#ifdef PHIST_HAVE_LIKWID
  LIKWID_MARKER_INIT;
  LIKWID_MARKER_START("phist<tpetra>");
#endif
}
      
  // finalize kernel library. Should at least call MPI_Finalize if it has not been called
  // but is required.
  extern "C" void phist_kernels_finalize(int* ierr)
  {
#ifdef PHIST_HAVE_LIKWID
    LIKWID_MARKER_STOP("phist<tpetra>");
    LIKWID_MARKER_CLOSE;
#endif
#ifdef PHIST_TIMEMONITOR
  phist_TimeMonitor::Timer::summarize();
#endif
#ifdef PHIST_HAVE_MPI
  if (myMpiSession==1)
  {
    *ierr=MPI_Finalize();
  }
#endif
}
            

//!
extern "C" void phist_comm_create(comm_ptr_t* vcomm, int* ierr)
{
  *ierr=0;
  *vcomm = (comm_ptr_t)(Teuchos::DefaultComm<int>::getComm().get());
}

//!
extern "C" void phist_comm_delete(comm_ptr_t vcomm, int* ierr)
{
  *ierr=0;
  TOUCH(vcomm);
  // note - as comm_create returns a raw pointer to the default comm, we should not delete it.
}


//!
extern "C" void phist_comm_get_rank(const_comm_ptr_t vcomm, int* rank, int* ierr)
{
  *ierr=0;
  CAST_PTR_FROM_VOID(const comm_t,comm,vcomm,*ierr);
  *rank=comm->getRank();
}
//!
extern "C" void phist_comm_get_size(const_comm_ptr_t vcomm, int* size, int* ierr)
{
  *ierr=0;
  CAST_PTR_FROM_VOID(const comm_t,comm,vcomm,*ierr);
  *size=comm->getSize();
}

// This function creates the default compute node to be used for vectors and matrices.
// It is not exposed to the user because the Kokkos node concept is specific for Tpetra
// and doesn't have a direct equivalent in epetra or ghost. The function checks if a file
// node.xml exists in which the node parameters like "Num Threads" can be set, otherwise
// it just uses default parameters.
extern "C" void phist_tpetra_node_create(node_t** node, const_comm_ptr_t vcomm, int* ierr)
{
  *ierr=0;
  CAST_PTR_FROM_VOID(const comm_t,comm,vcomm,*ierr);
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
      *ierr=PHIST_CAUGHT_EXCEPTION; 
      return;
      }
  }
  else
  {
    PHIST_SOUT(PHIST_WARNING,"File phist_node.xml not found, using default node settings in tpetra\n");
  }

  PHIST_SOUT(PHIST_VERBOSE,"# threads requested: %d\n",nodeParams->get("Num Threads",0));
  *node = new node_t(*nodeParams);
}

//!
extern "C" void phist_map_create(map_ptr_t* vmap, const_comm_ptr_t vcomm, gidx_t nglob, int *ierr)
{
  *ierr=0;
  CAST_PTR_FROM_VOID(const comm_t,comm,vcomm,*ierr);
  
  Teuchos::RCP<const comm_t> comm_ptr = Teuchos::rcp(comm,false);
node_t* node;
  PHIST_CHK_IERR(phist_tpetra_node_create(&node,vcomm,ierr),*ierr);
  Teuchos::RCP<node_t> node_ptr = Teuchos::rcp(node,true);
  Teuchos::RCP<const map_t> map_ptr =
        Tpetra::createUniformContigMapWithNode<lidx_t,gidx_t,node_t>
                (nglob, comm_ptr,node_ptr);
  *vmap = (map_ptr_t)(map_ptr.release().get());
}

//!
extern "C" void phist_map_delete(map_ptr_t vmap, int *ierr)
{
  *ierr=0;
  CAST_PTR_FROM_VOID(map_t,map,vmap,*ierr);
  delete map;
  vmap=NULL;
}

//!
extern "C" void phist_map_get_comm(const_map_ptr_t vmap, const_comm_ptr_t* vcomm, int* ierr)
{
  *ierr=0;
  CAST_PTR_FROM_VOID(const map_t,map,vmap,*ierr);
  Teuchos::RCP<const comm_t> comm = map->getComm();
  *vcomm = (const_comm_ptr_t)(comm.get());
}

//!
extern "C" void phist_map_get_local_length(const_map_ptr_t vmap, lidx_t* nloc, int* ierr)
{
  *ierr=0;
  CAST_PTR_FROM_VOID(const map_t,map,vmap,*ierr);
  *nloc = map->getNodeNumElements();
}

//! returns the smallest global index in the map appearing on my partition. ierr is set to 1
//! in case the map is not contiguous, because in that case it may be that the
//! caller falsely assumes global elements [ilower ... iupper] are actually on this partition.
extern "C" void phist_map_get_ilower(const_map_ptr_t vmap, gidx_t* ilower, int* ierr)
{
  *ierr=0;
  CAST_PTR_FROM_VOID(const map_t,map,vmap,*ierr);
  if (map->isContiguous()==false) *ierr=1;
  *ilower = map->getMinGlobalIndex();
}
//! returns the largest global index in the map appearing on my partition. ierr is set to 1
//! in case the map is not contiguous, because in that case it may be that the
//! caller falsely assumes global elements [ilower ... iupper] are actually on this partition.
extern "C" void phist_map_get_iupper(const_map_ptr_t vmap, gidx_t* iupper, int* ierr)
{
  *ierr=0;
  CAST_PTR_FROM_VOID(const map_t,map,vmap,*ierr);
  if (map->isContiguous()==false) *ierr=1;
  *iupper = map->getMaxGlobalIndex();
}

#ifdef PHIST_TIMEMONITOR
extern "C" void phist_totalMatVecCount()
{
  ENTER_FCN(__FUNCTION__);
}
#endif

#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "kernels_def.hpp"
#include "carp_def.hpp"

#include "phist_gen_c.h"
#include "kernels_def.hpp"
#include "carp_def.hpp"
#endif

#include "phist_gen_d.h"
#include "kernels_def.hpp"
#include "carp_def.hpp"

#include "phist_gen_z.h"
#include "kernels_def.hpp"
#include "carp_def.hpp"
