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

  phist_comm_ptr comm = NULL;
  sparseMat_ptr A = NULL;
  phist_const_map_ptr map = NULL;
  mvec_ptr xr = NULL,xi=NULL,rhs=NULL;
  
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
  
  int sparseMat_flags=PHIST_SPARSEMAT_OPT_CARP|PHIST_SPARSEMAT_PERM_GLOBAL;
  if (argc>4) sparseMat_flags=atoi(argv[4]);

  iflag=sparseMat_flags;
  PHIST_ICHK_IERR(SUBR(create_matrix)(&A,comm,matname,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(sparseMat_get_domain_map)(A, &map,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(mvec_create)(&xr,map,nvecs,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(mvec_random)(xr,&iflag),iflag);
  
  double sigma_r=shift_type? 42.9:0.0;
  double sigma_i=shift_type>1?-3.14159:0.0;
  double omega=1.5;
  void* work=NULL;
  
  if (shift_type==2)
  {
    PHIST_ICHK_IERR(SUBR(mvec_create)(&xi,map,nvecs,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(mvec_random)(xi,&iflag),iflag);
  }
  else
  {
    xi=NULL;
  }
  
  ST shifts[nvecs];
  MT shifts_r[nvecs],shifts_i[nvecs],omegas[nvecs];
  for (int i=0;i<nvecs;i++)
  {
    shifts_r[i]=sigma_r;
    shifts_i[i]=sigma_i;
    shifts[i]=shifts_r[i]+shifts_i[i]*st::cmplx_I();
    omegas[i]=omega;
  }

#ifdef IS_COMPLEX
  PHIST_ICHK_IERR(SUBR(carp_setup)(A,nvecs,shifts_r,shifts_i,&work,&iflag),iflag);
#else
  if (shift_type==2)
  {
    PHIST_ICHK_IERR(SUBR(carp_setup_rc)(A,nvecs,shifts_r,shifts_i,&work,&iflag),iflag);
  }
  else
  {
    PHIST_ICHK_IERR(SUBR(carp_setup)(A,nvecs,shifts,&work,&iflag),iflag);
  }
#endif
  
  // we only study the case rhs=NULL here because it is the dominant situation in CARP-CG
  // (the rhs is only passed in in the first sweep or for correction steps).
  for(int i = 0; i < NUM_CARP/nvecs; i++)
  {
#ifdef IS_COMPLEX
    PHIST_ICHK_IERR(SUBR(carp_sweep)(A,shifts_r,shifts_i,rhs,xr,xi,
        work, omegas, &iflag),iflag);
#else
    if (shift_type==2)
    {
      PHIST_ICHK_IERR(SUBR(carp_sweep_rc)(A,shifts_r,shifts_i,rhs,xr,xi,
        work, omegas, &iflag),iflag);
    }
    else
    {
      PHIST_ICHK_IERR(SUBR(carp_sweep)(A,shifts,rhs,xr,
        work, omegas, &iflag),iflag);
    }

#endif
  }
  
  PHIST_ICHK_IERR(SUBR(carp_destroy)(A,work,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(mvec_delete)(xr,&iflag),iflag);
  if (xi!=NULL)
  {
    PHIST_ICHK_IERR(SUBR(mvec_delete)(xi,&iflag),iflag);
  }
  PHIST_ICHK_IERR(SUBR(sparseMat_delete)(A,&iflag),iflag);

#ifndef PHIST_KERNEL_LIB_MAGMA
  double max_bw,mean_bw;
  PHIST_ICHK_IERR(phist_bench_stream_triad(&mean_bw,&max_bw, &iflag),iflag);
  PHIST_SOUT(PHIST_INFO,"Maximum memory bandwidth (STREAM_TRIAD): %g GB/s\n",max_bw*1.0e-9);
#endif
PHIST_MAIN_TASK_END
  PHIST_ICHK_IERR(phist_kernels_finalize(&iflag),iflag);

  return iflag;
}
