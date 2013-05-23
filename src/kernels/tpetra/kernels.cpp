#include "typedefs.hpp"
#include "../phist_kernels.h"

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
void phist_comm_get_rank(const_comm_ptr_t comm, int* rank, int* ierr)
  {
  *ierr=-99;
  }
//!
void phist_comm_get_size(const_comm_ptr_t comm, int* size, int* ierr)
  {
  *ierr=-99;
  }
//!
void phist_map_create(map_ptr_t* map, const_comm_ptr_t comm, int nglob, int *ierr)
  {
  *ierr=-99;
  }
//!
void phist_map_get_comm(const_map_ptr_t map, const_comm_ptr_t* comm, int* ierr)
  {
  *ierr=-99;
  }

} // extern "C"

#include "../gen_s.h"
#include "kernels_def.hpp"

#include "../gen_d.h"
#include "kernels_def.hpp"

#include "../gen_c.h"
#include "kernels_def.hpp"

#include "../gen_z.h"
#include "kernels_def.hpp"
