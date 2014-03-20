#include "../phist_kernels.h"
#include "phist_trilinos_macros.h"
#include "phist_typedefs.h"
#include "Epetra_config.h"
#ifdef PHIST_HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#ifdef PHIST_TIMEMONITOR
#include <Teuchos_TimeMonitor.hpp>
#endif
#include "epetra_helpers.cpp"
#include "phist_ScalarTraits.hpp"

#ifdef PHIST_HAVE_LIKWID
#include <likwid.h>
#endif

extern "C" {

// initialize kernel library. Should at least call MPI_Init if it has not been called
// but is required.
void phist_kernels_init(int* argc, char*** argv, int* ierr)
  {
  *ierr=0;
#ifdef PHIST_HAVE_MPI
  *ierr=MPI_Init(argc,argv); 
#endif
#ifdef PHIST_HAVE_LIKWID
  LIKWID_MARKER_INIT;
  LIKWID_MARKER_START("phist<epetra>");
#endif
  }

// finalize kernel library. Should at least call MPI_Finalize if it has not been called
// but is required.
void phist_kernels_finalize(int* ierr)
  {
#ifdef PHIST_HAVE_LIKWID
  LIKWID_MARKER_STOP("phist<epetra>");
  LIKWID_MARKER_CLOSE;
#endif
#ifdef PHIST_TIMEMONITOR
  Teuchos::TimeMonitor::summarize();
#endif
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
void phist_comm_delete(comm_ptr_t vcomm, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(Epetra_Comm,comm,vcomm,*ierr);
  delete comm;
  vcomm=NULL;
  }

//!
void phist_comm_get_rank(const_comm_ptr_t vcomm, int* rank, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(const Epetra_Comm,comm,vcomm,*ierr);
  *rank=comm->MyPID();
  }
//!
void phist_comm_get_size(const_comm_ptr_t vcomm, int* size, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(const Epetra_Comm,comm,vcomm,*ierr);
  *size=comm->NumProc();
  }
//!
void phist_map_create(map_ptr_t* vmap, const_comm_ptr_t vcomm, gidx_t nglob, int *ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(const Epetra_Comm,comm,vcomm,*ierr);
  Epetra_BlockMap* map;
  TRY_CATCH(map = new Epetra_Map(nglob,0,*comm),*ierr);
  *vmap=(map_ptr_t)(map);
  }

//!
void phist_map_delete(map_ptr_t vmap, int *ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(Epetra_BlockMap,map,vmap,*ierr);
  delete map;
  vmap=NULL;
  }
  
//!
void phist_map_get_comm(const_map_ptr_t vmap, const_comm_ptr_t* vcomm, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(const Epetra_BlockMap,map,vmap,*ierr);
  *vcomm = (const_comm_ptr_t)(&map->Comm());  
  }

//!
void phist_map_get_local_length(const_map_ptr_t vmap, int* nloc, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(const Epetra_BlockMap,map,vmap,*ierr);
  *nloc = map->NumMyElements();
  }

//! returns the smallest global index in the map appearing on my partition. ierr is set to 1
//! in case the map is not contiguous, because in that case it may be that the
//! caller falsely assumes global elements [ilower ... iupper] are actually on this partition.
void phist_map_get_ilower(const_map_ptr_t vmap, gidx_t* ilower, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(const Epetra_BlockMap,map,vmap,*ierr);
  if (map->LinearMap()==false) *ierr=1;
  *ilower = map->MinMyGID();
  }
//! returns the largest global index in the map appearing on my partition. ierr is set to 1
//! in case the map is not contiguous, because in that case it may be that the
//! caller falsely assumes global elements [ilower ... iupper] are actually on this partition.
void phist_map_get_iupper(const_map_ptr_t vmap, gidx_t* iupper, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(const Epetra_BlockMap,map,vmap,*ierr);
  if (map->LinearMap()==false) *ierr=1;
  *iupper = map->MaxMyGID();
  }

#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "../kernels_noimpl.c"
#include "../carp_noimpl.c"

#include "phist_gen_c.h"
#include "../kernels_noimpl.c"
#include "../carp_noimpl.c"
#endif
#include "phist_gen_z.h"
#include "../kernels_noimpl.c"
#include "../carp_noimpl.c"

} //extern "C"

#include "phist_gen_d.h"
#include "kernels_def.hpp"
#include "carp_def.hpp"

