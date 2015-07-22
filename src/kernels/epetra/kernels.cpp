#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "../phist_kernels.h"
#include "phist_trilinos_macros.h"
#include "phist_typedefs.h"
#include "phist_typedefs.h"
#include "Epetra_config.h"
#ifdef PHIST_HAVE_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_SerialComm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"

#include "Teuchos_StandardCatchMacros.hpp"
#include "EpetraExt_CrsMatrixIn.h"

#include "BelosEpetraAdapter.hpp"
#include "BelosTsqrOrthoManager.hpp"

#include "epetra_helpers.cpp"
#include "phist_ScalarTraits.hpp"

#ifdef PHIST_HAVE_LIKWID
#include <likwid.h>
#endif

extern "C" {

// initialize kernel library. Should at least call MPI_Init if it has not been called
// but is required.
extern "C" void phist_kernels_init(int* argc, char*** argv, int* iflag)
  {
  *iflag=0;
#ifdef PHIST_HAVE_MPI
  *iflag=MPI_Init(argc,argv); 
#endif
#ifdef PHIST_HAVE_LIKWID
  LIKWID_MARKER_INIT;
  LIKWID_MARKER_START("phist<epetra>");
#endif
  }

// finalize kernel library. Should at least call MPI_Finalize if it has not been called
// but is required.
extern "C" void phist_kernels_finalize(int* iflag)
  {
#ifdef PHIST_HAVE_LIKWID
  LIKWID_MARKER_STOP("phist<epetra>");
  LIKWID_MARKER_CLOSE;
#endif
#if defined(PHIST_TIMEMONITOR) || defined(PHIST_TIMEMONITOR_PERLINE)
PHIST_CXX_TIMER_SUMMARIZE;
#endif
#ifdef PHIST_HAVE_MPI
  *iflag=MPI_Finalize();
#endif  
  }


//!
extern "C" void phist_comm_create(comm_ptr_t* vcomm, int* iflag)
  {
  *iflag=0;
#ifdef PHIST_HAVE_MPI
  Epetra_Comm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_Comm *comm = new Epetra_SerialComm;
#endif
  if (comm==NULL) *iflag=-1;
  *vcomm=(comm_ptr_t)(comm);
  }

//!
extern "C" void phist_comm_delete(comm_ptr_t vcomm, int* iflag)
  {
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(Epetra_Comm,comm,vcomm,*iflag);
  delete comm;
  vcomm=NULL;
  }

#ifdef PHIST_HAVE_MPI
extern void phist_comm_get_mpi_comm(const_comm_ptr_t vcomm, MPI_Comm* mpiComm, int* iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_MpiComm,comm,vcomm,*iflag);
  *mpiComm=comm->Comm();
}
#endif

//!
extern "C" void phist_comm_get_rank(const_comm_ptr_t vcomm, int* rank, int* iflag)
  {
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_Comm,comm,vcomm,*iflag);
  *rank=comm->MyPID();
  }
//!
extern "C" void phist_comm_get_size(const_comm_ptr_t vcomm, int* size, int* iflag)
  {
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_Comm,comm,vcomm,*iflag);
  *size=comm->NumProc();
  }
//!
extern "C" void phist_map_create(map_ptr_t* vmap, const_comm_ptr_t vcomm, gidx_t nglob, int *iflag)
  {
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_Comm,comm,vcomm,*iflag);
  Epetra_BlockMap* map;
  PHIST_TRY_CATCH(map = new Epetra_Map(nglob,0,*comm),*iflag);
  *vmap=(map_ptr_t)(map);
  }

//!
extern "C" void phist_map_delete(map_ptr_t vmap, int *iflag)
  {
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(Epetra_BlockMap,map,vmap,*iflag);
  delete map;
  vmap=NULL;
  }
  
//!
extern "C" void phist_map_get_comm(const_map_ptr_t vmap, const_comm_ptr_t* vcomm, int* iflag)
  {
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_BlockMap,map,vmap,*iflag);
  *vcomm = (const_comm_ptr_t)(&map->Comm());  
  }

//!
extern "C" void phist_map_get_local_length(const_map_ptr_t vmap, lidx_t* nloc, int* iflag)
  {
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_BlockMap,map,vmap,*iflag);
  *nloc = map->NumMyElements();
  }

//! returns the smallest global index in the map appearing on my partition. iflag is set to 1
//! in case the map is not contiguous, because in that case it may be that the
//! caller falsely assumes global elements [ilower ... iupper] are actually on this partition.
extern "C" void phist_map_get_ilower(const_map_ptr_t vmap, gidx_t* ilower, int* iflag)
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
extern "C" void phist_map_get_iupper(const_map_ptr_t vmap, gidx_t* iupper, int* iflag)
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

extern "C" void phist_bench_stream_load(double* bw, int* iflag)
{
  *iflag = PHIST_NOT_IMPLEMENTED;
}
extern "C" void phist_bench_stream_store(double* bw, int* iflag)
{
  *iflag = PHIST_NOT_IMPLEMENTED;
}
extern "C" void phist_bench_stream_triad(double* bw, int* iflag)
{
  *iflag = PHIST_NOT_IMPLEMENTED;
}

#ifdef PHIST_TIMEMONITOR
extern "C" void phist_totalMatVecCount()
{
  PHIST_ENTER_FCN(__FUNCTION__);
}
#endif


#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "../kernels_noimpl.c"
#include "../kernels_nogpu.c"
#include "../carp_noimpl.c"

#include "phist_gen_c.h"
#include "../kernels_noimpl.c"
#include "../kernels_nogpu.c"
#include "../carp_noimpl.c"
#endif
#include "phist_gen_z.h"
#include "../kernels_noimpl.c"
#include "../kernels_nogpu.c"
#include "../carp_noimpl.c"

} //extern "C"

#include "phist_gen_d.h"
#include "kernels_def.hpp"
#include "carp_def.hpp"
#include "kernels_no_fused.cpp"

