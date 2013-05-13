#include "typedefs.hpp"
#include "../kernels.h"

//!
void comm_create(comm_ptr_t* vcomm, int* ierr)
  {
  *ierr=0;
  *vcomm = (comm_ptr_t)(Teuchos::DefaultComm<int>::getComm().get());
  }
//!
void comm_get_rank(const_comm_ptr_t* comm, int* rank, int* ierr)
  {
  *ierr=-99;
  }
//!
void comm_get_size(const_comm_ptr_t* comm, int* size, int* ierr);
  {
  *ierr=-99;
  }
//!
void map_create(map_ptr_t* map, const_comm_ptr_t comm, int nglob, int *ierr);
  {
  *ierr=-99;
  }
//!
void map_get_comm(const_map_ptr_t map, const_comm_ptr_t* comm, int* ierr)
  {
  *ierr=-99;
  }


#include "../gen_s.h"
#include "kernels_def.hpp"

#include "../gen_d.h"
#include "kernels_def.hpp"

#include "../gen_c.h"
#include "kernels_def.hpp"

#include "../gen_z.h"
#include "kernels_def.hpp"
