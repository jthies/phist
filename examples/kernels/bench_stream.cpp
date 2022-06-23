/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "phist_kernels.h"
#include "phist_macros.h"

// macro to do several barriers to make sure output is finished on a rank
#define IO_BARRIER \
for (int _io_barrier_i=0; _io_barrier_i<20; _io_barrier_i++) \
{fflush(stdout); MPI_Barrier(MPI_COMM_WORLD);}

/* This driver just runs three STREAM benchmarks to measure the memory bandwidth
   achieved by each MPI Process. The benchmarks are LOAD, STORE and TRIAD.
   */
int main(int argc, char** argv)
{
#ifdef PHIST_HAVE_MPI
  int iflag = 0;
  PHIST_ICHK_IERR(phist_kernels_init(&argc,&argv,&iflag),iflag);
PHIST_MAIN_TASK_BEGIN

#ifndef PHIST_KERNEL_LIB_MAGMA
  double max_bw_load,  max_bw_store,  max_bw_copy,  max_bw_triad;
  double mean_bw_load, mean_bw_store, mean_bw_copy, mean_bw_triad;
  PHIST_ICHK_IERR(phist_bench_stream_load(&mean_bw_load,&max_bw_load,&iflag),iflag);
  PHIST_ICHK_IERR(phist_bench_stream_store(&mean_bw_store,&max_bw_store,&iflag),iflag);
  PHIST_ICHK_IERR(phist_bench_stream_copy(&mean_bw_copy,&max_bw_copy,&iflag),iflag);
  PHIST_ICHK_IERR(phist_bench_stream_triad(&mean_bw_triad,&max_bw_triad,&iflag),iflag);
#else
  mean_bw_load=0.0;  max_bw_load=0.0;
  mean_bw_store=0.0; max_bw_store=0.0;
  mean_bw_copy=0.0; max_bw_copy=0.0;
  mean_bw_triad=0.0; max_bw_triad=0.0;
  PHIST_SOUT(PHIST_WARNING,"benchmark not implemented correctly with MAGMA\n");
#endif

// output in GB/s
max_bw_load/=1.0e9;
max_bw_store/=1.0e9;
max_bw_copy/=1.0e9;
max_bw_triad/=1.0e9;

mean_bw_load/=1.0e9;
mean_bw_store/=1.0e9;
mean_bw_copy/=1.0e9;
mean_bw_triad/=1.0e9;


// print measured bandwidths in a table. Note that each MPI process 
// may run a different CPU component or GPU, so the full picture may
// be interesting.
int rank, size, comm_size;
MPI_Comm_rank(MPI_COMM_WORLD,&rank);
MPI_Comm_size(MPI_COMM_WORLD,&comm_size);
size=comm_size;

// collect the results for printing them
double all_mean_store[comm_size], all_mean_load[comm_size], all_mean_copy[comm_size], all_mean_triad[comm_size];

MPI_Gather(&mean_bw_load,1,MPI_DOUBLE,
           all_mean_load,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

MPI_Gather(&mean_bw_store,1,MPI_DOUBLE,
           all_mean_store,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

MPI_Gather(&mean_bw_copy,1,MPI_DOUBLE,
           all_mean_copy,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

MPI_Gather(&mean_bw_triad,1,MPI_DOUBLE,
           all_mean_triad,1,MPI_DOUBLE,0,MPI_COMM_WORLD);


if (size>10)
{
  PHIST_SOUT(PHIST_WARNING,"showing results of first 10 MPI ranks only\n");
  size=10;
}

if (rank==0)
{
  PHIST_SOUT(PHIST_INFO,"\nSTREAM BENCHMARK RESULTS (GB/s) per MPI rank\n");
  PHIST_SOUT(PHIST_INFO,"=============");
  for (int i=0; i<size; i++) PHIST_SOUT(PHIST_INFO,"========");
  PHIST_SOUT(PHIST_INFO,"\n| benchmark |");
  for (int i=0; i<size; i++) PHIST_SOUT(PHIST_INFO,"  P%d   |",i);
  PHIST_SOUT(PHIST_INFO,"\n-------------");
  for (int i=0; i<size; i++) PHIST_SOUT(PHIST_INFO,"--------");

  PHIST_SOUT(PHIST_INFO,"\n|  load     |");
  for (int i=0; i<size; i++)
  {
    fprintf(stdout," %5.3g |",all_mean_load[i]);
  }
  PHIST_SOUT(PHIST_INFO,"\n|  store    |");
  for (int i=0; i<size; i++)
  {
    fprintf(stdout," %5.3g |",all_mean_store[i]);
  }

  PHIST_SOUT(PHIST_INFO,"\n|  copy    |");
  for (int i=0; i<size; i++)
  {
    fprintf(stdout," %5.3g |",all_mean_copy[i]);
  }

  PHIST_SOUT(PHIST_INFO,"\n|  triad    |");
  for (int i=0; i<size; i++)
  {
    fprintf(stdout," %5.3g |",all_mean_triad[i]);
  }

  PHIST_SOUT(PHIST_INFO,"\n=============");
  for (int i=0; i<size; i++) PHIST_SOUT(PHIST_INFO,"========");

  // deduce process weights from the results

  double bw_triad_total=0.0;
  for (int i=0; i<comm_size; i++) bw_triad_total+=all_mean_triad[i];
  PHIST_SOUT(PHIST_INFO,"\n\nsuggested weights based on stream triad benchmark:\n");
  for (int i=0; i<comm_size; i++)
  {
    PHIST_SOUT(PHIST_INFO,"P%d: %4.2f\n",i,all_mean_triad[i]/bw_triad_total);
  }
}

PHIST_MAIN_TASK_END
  PHIST_ICHK_IERR(phist_kernels_finalize(&iflag),iflag);
#endif
}
