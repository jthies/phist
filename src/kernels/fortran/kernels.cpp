#include "phist_macros.h"
#include "phist_trilinos_macros.h"
#include "../phist_kernels.h"

#include "phist_typedefs.h"
#include "typedefs.hpp"
#include "phist_ScalarTraits.hpp"

#ifdef PHIST_TIMEMONITOR
#include <Teuchos_TimeMonitor.hpp>
#endif

#ifdef PHIST_HAVE_LIKWID
#include <likwid.h>
#endif

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

extern "C" {

#ifndef PHIST_HAVE_MPI
#error "fortran kernels only work with MPI"
#endif

// initialize
void phist_kernels_init(int* argc, char*** argv, int* ierr)
{
  *ierr=0;
  PHIST_CHK_IERR( *ierr = MPI_Init(argc, argv), *ierr);
#ifdef PHIST_HAVE_LIKWID
  LIKWID_MARKER_INIT;
  LIKWID_MARKER_START("phist<ghost>");
#endif
}

// finalize ghost
void phist_kernels_finalize(int* ierr)
{
#ifdef PHIST_HAVE_LIKWID
  LIKWID_MARKER_STOP("phist<ghost>");
  LIKWID_MARKER_CLOSE;
#endif
#ifdef PHIST_TIMEMONITOR
  Teuchos::TimeMonitor::summarize();
#endif
  PHIST_CHK_IERR( *ierr = MPI_Finalize(), *ierr);
  *ierr=0;
}

void phist_comm_create(comm_ptr_t* vcomm, int* ierr);
void phist_comm_delete(comm_ptr_t vcomm, int* ierr);
void phist_comm_get_rank(const_comm_ptr_t vcomm, int* rank, int* ierr);
void phist_comm_get_size(const_comm_ptr_t vcomm, int* size, int* ierr);
void phist_map_create(map_ptr_t* vmap, const_comm_ptr_t vcomm, gidx_t nglob, int *ierr);
void phist_map_delete(map_ptr_t vmap, int *ierr);
void phist_map_get_comm(const_map_ptr_t vmap, const_comm_ptr_t* vcomm, int* ierr);
void phist_map_get_local_length(const_map_ptr_t vmap, int* nloc, int* ierr);
void phist_map_get_ilower(const_map_ptr_t vmap, int* ilower, int* ierr);
void phist_map_get_iupper(const_map_ptr_t vmap, int* iupper, int* ierr);


#include "phist_gen_s.h"
#include "../kernels_noimpl.c"

#include "phist_gen_c.h"
#include "../kernels_noimpl.c"

#include "phist_gen_z.h"
#include "../kernels_noimpl.c"

} //extern "C"

#include "phist_gen_d.h"
#include "kernels_def.hpp"

