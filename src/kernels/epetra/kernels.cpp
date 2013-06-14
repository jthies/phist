#include "../phist_kernels.h"
#include "../cpp_macros.h"
#include "Epetra_config.h"
#ifdef PHIST_HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "epetra_helpers.cpp"

extern "C" {

// initialize kernel library. Should at least call MPI_Init if it has not been called
// but is required.
void phist_kernels_init(int* argc, char*** argv, int* ierr)
  {
  *ierr=0;
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
#ifdef PHIST_HAVE_MPI
  Epetra_Comm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_Comm *comm = new Epetra_SerialComm;
#endif
  if (comm==NULL) *ierr=-1;
  *vcomm=(comm_ptr_t)(comm);
  }
//!
void phist_comm_get_rank(const_comm_ptr_t vcomm, int* rank, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Epetra_Comm,comm,vcomm,*ierr);
  *rank=comm->MyPID();
  }
//!
void phist_comm_get_size(const_comm_ptr_t vcomm, int* size, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Epetra_Comm,comm,vcomm,*ierr);
  *size=comm->NumProc();
  }
//!
void phist_map_create(map_ptr_t* vmap, const_comm_ptr_t vcomm, int nglob, int *ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Epetra_Comm,comm,vcomm,*ierr);
  Epetra_BlockMap* map;
  _TRY_CATCH_(map = new Epetra_Map(nglob,0,*comm),*ierr);
  *vmap=(map_ptr_t)(map);
  }
//!
void phist_map_get_comm(const_map_ptr_t vmap, const_comm_ptr_t* vcomm, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Epetra_BlockMap,map,vmap,*ierr);
  *vcomm = (const_comm_ptr_t)(&map->Comm());  
  }

#include "phist_gen_s.h"
#include "../kernels_noimpl.c"

#include "phist_gen_c.h"
#include "../kernels_noimpl.c"

#include "phist_gen_z.h"
#include "../kernels_noimpl.c"

} //extern "C"

#include "phist_gen_d.h"
#include "kernels_def.hpp"

