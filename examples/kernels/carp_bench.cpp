#include "phist_config.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#ifdef PHIST_HAVE_LIKWID
#include "likwid.h"
#endif

#include "phist_macros.h"
#include "phist_enums.h"
#include "phist_kernels.h"
#include "phist_operator.h"
#include "phist_gen_d.h"
#include "phist_driver_utils.h"
#include "phist_ScalarTraits.hpp"

typedef phist::ScalarTraits<_ST_> st;

#define NUM_CARP 360

// benchmark application for the CARP kernel in DP with one or more vectors and
// with or without a real or complex shift.
int main(int argc, char** argv)
{
  int rank, num_proc;
  int iflag;
  int verbose;

  comm_ptr_t comm = NULL;
  sparseMat_ptr_t A = NULL;
  const_map_ptr_t map = NULL;
  mvec_ptr_t xr = NULL,xi=NULL,rhs=NULL;
  
  PHIST_ICHK_IERR(phist_kernels_init(&argc,&argv,&iflag),iflag);
PHIST_MAIN_TASK_BEGIN
  PHIST_ICHK_IERR(phist_comm_create(&comm,&iflag),iflag);

  PHIST_ICHK_IERR(phist_comm_get_rank(comm, &rank,&iflag),iflag);
  PHIST_ICHK_IERR(phist_comm_get_size(comm, &num_proc,&iflag),iflag);

  verbose= (rank==0);

  if (argc<2)
  {
    if (verbose) fprintf(stdout,"Usage: %s <matrix filename|problem string> [<nvecs> [<shift-type> [mat_flags]]]\n"
                                "       apply a couple of times X<-DKSWP(A-sI,X)\n"
                                "       nvecs defaults to 1\n"
                                "       shift_type 0: no shift, 1/2: real/complex shift (default 0)\n"
                                "       mat_flags: flags passed to the matrix creation routines. Default: PHIST_SPARSEMAT_REPART|PHIST_SPARSEMAT_DIST2_COLOR.\n"
                                "       \n"
                                "       The flag can be constructed by bitwise OR of the following:\n"
                                PHIST_SPARSEMAT_FLAGS_DESCRIPTION
                                "       \n",
        argv[0]);
    // print usage message for creating/reading a matrix
    SUBR(create_matrix)(NULL,NULL,"usage",&iflag);        
    return 1;
  }

  const char* matname=argv[1];
  
  int nvecs, shift_type;
  nvecs=1;
  shift_type=0;
  if (argc>2) nvecs=atoi(argv[2]);
  if (argc>3) shift_type=atoi(argv[3]);
  
  int sparseMat_flags=PHIST_SPARSEMAT_OPT_CARP|PHIST_SPARSEMAT_REPARTITION;
  if (argc>4) sparseMat_flags=atoi(argv[4]);

  iflag=sparseMat_flags;
  PHIST_ICHK_IERR(SUBR(create_matrix)(&A,comm,matname,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(sparseMat_get_domain_map)(A, &map,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(mvec_create)(&xr,map,nvecs,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(mvec_random)(xr,&iflag),iflag);
  
  int numShifts=1;
  double sigma_r=shift_type? 42.9:0.0;
  double sigma_i=shift_type>1?-3.14159:0.0;
  double omega=1.5;
  double* nrms_ai2i=NULL;
  void* work=NULL;
  
  PHIST_ICHK_IERR(SUBR(carp_setup)(A,numShifts,&sigma_r,&sigma_i,&nrms_ai2i,&work,&iflag),iflag);
  
  if (shift_type==2)
  {
    PHIST_ICHK_IERR(SUBR(mvec_create)(&xi,map,nvecs,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(mvec_random)(xi,&iflag),iflag);
  }
  else
  {
    xi=NULL;
  }
  

  // we only study the case rhs=NULL here because it is the dominant situation in CARP-CG
  // (the rhs is only passed in in the first sweep or for correction steps).
  for(int i = 0; i < NUM_CARP/nvecs; i++)
  {
    PHIST_ICHK_IERR(SUBR(carp_sweep)(A,numShifts,&sigma_r,&sigma_i,rhs,&xr,&xi,
    nrms_ai2i, work, &omega, &iflag),iflag);
  }
  
  PHIST_ICHK_IERR(SUBR(carp_destroy)(A,numShifts,nrms_ai2i,work,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(mvec_delete)(xr,&iflag),iflag);
  if (xi!=NULL)
  {
    PHIST_ICHK_IERR(SUBR(mvec_delete)(xi,&iflag),iflag);
  }
  PHIST_ICHK_IERR(SUBR(sparseMat_delete)(A,&iflag),iflag);

  double max_bw;
  PHIST_ICHK_IERR(phist_bench_stream_triad(&max_bw, &iflag),iflag);
  PHIST_SOUT(PHIST_INFO,"Maximum memory bandwidth (STREAM_TRIAD): %g GB/s\n",max_bw*1.0e-9);
PHIST_MAIN_TASK_END
  PHIST_ICHK_IERR(phist_kernels_finalize(&iflag),iflag);

  return iflag;
}
