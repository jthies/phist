#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

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

#include <malloc.h>

extern "C" {

#ifndef PHIST_HAVE_MPI
#error "fortran kernels only work with MPI"
#endif


// comment in for glibc/gcc and memory alignment problems...
/*
// force all memory allocations to be 16 byte aligned!
static void*(*old_malloc_hook)(size_t,const void*);
static void * my_malloc_hook (size_t size, const void *caller)
{
   void *result;
   // Restore all old hooks
   __malloc_hook = old_malloc_hook;
   // Call recursively
   result = memalign (16, size);
   // Save underlying hooks
   old_malloc_hook = __malloc_hook;
   // printf might call malloc, so protect it too.
   printf ("malloc (%u) returns %p\n", (unsigned int) size, result);
   // Restore our own hooks
   __malloc_hook = my_malloc_hook;
   return result;
}
static void my_init_hook (void)
{
  old_malloc_hook = __malloc_hook;
  __malloc_hook = my_malloc_hook;
}
// Override initializing hook from the C library.
void (*__malloc_initialize_hook) (void) = my_init_hook;
*/



// initialize
void phist_kernels_init(int* argc, char*** argv, int* ierr)
{
  *ierr=0;
  int initialized = 0;

  PHIST_CHK_IERR( *ierr = MPI_Initialized(&initialized), *ierr);
  if (!initialized) {
      PHIST_CHK_IERR( *ierr = MPI_Init(argc, argv), *ierr);
  }
#ifdef PHIST_HAVE_LIKWID
  LIKWID_MARKER_INIT;
#pragma omp parallel
  {
    LIKWID_MARKER_THREADINIT;
    LIKWID_MARKER_START("phist<fortran>");
  }
#endif
}

// finalize fortran
void phist_kernels_finalize(int* ierr)
{
#ifdef PHIST_HAVE_LIKWID
#pragma omp parallel
  {
    LIKWID_MARKER_STOP("phist<fortran>");
  }
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
void phist_map_get_ilower(const_map_ptr_t vmap, gidx_t* ilower, int* ierr);
void phist_map_get_iupper(const_map_ptr_t vmap, gidx_t* iupper, int* ierr);

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
#include "../carp_noimpl.c"

