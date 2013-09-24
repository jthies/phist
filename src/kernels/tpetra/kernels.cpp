#include "phist_macros.h"
#include "typedefs.hpp"
#include "../phist_kernels.h"
#include "phist_trilinos_macros.h"
#include "phist_ScalarTraits.hpp"

#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_RCP.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Tpetra_MatrixIO.hpp"

#include "BelosTpetraAdapter.hpp"
#include "BelosTsqrOrthoManager.hpp"



extern "C" {

// initialize kernel library. Should at least call MPI_Init if it has not been called
// but is required.
void phist_kernels_init(int* argc, char*** argv, int* ierr)
  {
#ifdef PHIST_HAVE_MPI
  *ierr=MPI_Init(argc,argv);
#endif
  }
      
  // finalize kernel library. Should at least call MPI_Finalize if it has not been called
  // but is required.
  void phist_kernels_finalize(int* ierr)
    {
#ifdef PHIST_HAVE_MPI
  *ierr=MPI_Finalize();
#endif
    }
            

//!
void phist_comm_create(comm_ptr_t* vcomm, int* ierr)
  {
  *ierr=0;
  *vcomm = (comm_ptr_t)(Teuchos::DefaultComm<int>::getComm().get());
  }

//!
void phist_comm_delete(comm_ptr_t vcomm, int* ierr)
  {
  *ierr=0;
  // note - as comm_create returns a raw pointer to the default comm, we should not delete it.
  }


//!
void phist_comm_get_rank(const_comm_ptr_t vcomm, int* rank, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const comm_t,comm,vcomm,*ierr);
  *rank=comm->getRank();
  }
//!
void phist_comm_get_size(const_comm_ptr_t vcomm, int* size, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const comm_t,comm,vcomm,*ierr);
  *size=comm->getSize();
  }
//!
void phist_map_create(map_ptr_t* vmap, const_comm_ptr_t vcomm, gidx_t nglob, int *ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const comm_t,comm,vcomm,*ierr);
  map_t* map = new map_t(nglob,0,Teuchos::rcp(comm,false));
  *vmap = (map_ptr_t)(map);
  }

//!
void phist_map_delete(map_ptr_t vmap, int *ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(map_t,map,vmap,*ierr);
  delete map;
  vmap=NULL;
  }

//!
void phist_map_get_comm(const_map_ptr_t vmap, const_comm_ptr_t* vcomm, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const map_t,map,vmap,*ierr);
  Teuchos::RCP<const comm_t> comm = map->getComm();
  *vcomm = (const_comm_ptr_t)(comm.get());
  }

//!
void phist_map_get_local_length(const_map_ptr_t vmap, int* nloc, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const map_t,map,vmap,*ierr);
  *nloc = map->getNodeNumElements();
  }

//! returns the smallest global index in the map appearing on my partition. ierr is set to 1
//! in case the map is not contiguous, because in that case it may be that the
//! caller falsely assumes global elements [ilower ... iupper] are actually on this partition.
void phist_map_get_ilower(const_map_ptr_t vmap, int* ilower, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const map_t,map,vmap,*ierr);
  if (map->isContiguous()==false) *ierr=1;
  *ilower = map->getMinGlobalIndex();
  }
//! returns the largest global index in the map appearing on my partition. ierr is set to 1
//! in case the map is not contiguous, because in that case it may be that the
//! caller falsely assumes global elements [ilower ... iupper] are actually on this partition.
void phist_map_get_iupper(const_map_ptr_t vmap, int* iupper, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const map_t,map,vmap,*ierr);
  if (map->isContiguous()==false) *ierr=1;
  *iupper = map->getMaxGlobalIndex();
  }

} // extern "C"

#include "phist_gen_s.h"
#include "kernels_def.hpp"

#include "phist_gen_d.h"
#include "kernels_def.hpp"

#include "phist_gen_c.h"
#include "kernels_def.hpp"

#include "phist_gen_z.h"
#include "kernels_def.hpp"
