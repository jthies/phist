#include "typedefs.hpp"
#include "../phist_kernels.h"
#include "../cpp_macros.h"

#include "Teuchos_DefaultComm.hpp"

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
void phist_map_create(map_ptr_t* vmap, const_comm_ptr_t vcomm, int nglob, int *ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const comm_t,comm,vcomm,*ierr);
  map_t* map = new map_t(nglob,0,Teuchos::rcp(comm,false));
  *vmap = (map_ptr_t)(map);
  }
//!
void phist_map_get_comm(const_map_ptr_t vmap, const_comm_ptr_t* vcomm, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const map_t,map,vmap,*ierr);
  Teuchos::RCP<const comm_t> comm = map->getComm();
  *vcomm = (const_comm_ptr_t)(comm.get());
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
