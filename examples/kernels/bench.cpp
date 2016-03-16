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

int main(int argc, char** argv)
{
  int iflag = 0;
  PHIST_ICHK_IERR(phist_kernels_init(&argc,&argv,&iflag),iflag);
PHIST_MAIN_TASK_BEGIN

#ifndef PHIST_KERNEL_LIB_MAGMA
  double bw_load,bw_store,bw_triad;
  PHIST_ICHK_IERR(phist_bench_stream_load(&bw_load,&iflag),iflag);
  PHIST_ICHK_IERR(phist_bench_stream_store(&bw_store,&iflag),iflag);
  PHIST_ICHK_IERR(phist_bench_stream_triad(&bw_triad,&iflag),iflag);
#else
  bw_load=0.0;
  bw_store=0.0;
  bw_triad=0.0;
  PHIST_SOUT(PHIST_WARNING,"benchmark not implemented correctly with MAGMA\n");
#endif

// output in GB/s
bw_load/=1.0e9;
bw_store/=1.0e9;
bw_triad/=1.0e9;


// print measured bandwidths in a table. Note that each MPI process 
// may run a different CPU component or GPU, so the full picture may
// be interesting.
int rank, size, comm_size;
MPI_Comm_rank(MPI_COMM_WORLD,&rank);
MPI_Comm_size(MPI_COMM_WORLD,&comm_size);
size=comm_size;

if (size>8)
{
  PHIST_SOUT(PHIST_WARNING,"showing results of first 10 MPI ranks only\n");
  size=10;
}

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
  if (rank==i) fprintf(stdout," %5.3g |",bw_load);
  IO_BARRIER;
}
PHIST_SOUT(PHIST_INFO,"\n|  store    |");
for (int i=0; i<size; i++)
{ 
  if (rank==i) fprintf(stdout," %5.3g |",bw_store);
  IO_BARRIER;
}

PHIST_SOUT(PHIST_INFO,"\n|  triad    |");
for (int i=0; i<size; i++)
{ 
  if (rank==i) fprintf(stdout," %5.3g |",bw_triad);
  IO_BARRIER;
}

PHIST_SOUT(PHIST_INFO,"\n=============");
for (int i=0; i<size; i++) PHIST_SOUT(PHIST_INFO,"========");

// deduce process weights from the results
double all_store[comm_size], all_load[comm_size], all_triad[comm_size];

MPI_Gather(&bw_load,1,MPI_DOUBLE,
           all_load,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

MPI_Gather(&bw_store,1,MPI_DOUBLE,
           all_store,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

MPI_Gather(&bw_triad,1,MPI_DOUBLE,
           all_triad,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

if (rank==0)
{
  double bw_triad_total=0.0;
  for (int i=0; i<comm_size; i++) bw_triad_total+=all_triad[i];
  PHIST_SOUT(PHIST_INFO,"\n\nsuggested weights based on stream triad benchmark:\n");
  for (int i=0; i<comm_size; i++)
  {
    PHIST_SOUT(PHIST_INFO,"P%d: %4.2f\n",i,all_triad[i]/bw_triad_total);
  }
}

PHIST_MAIN_TASK_END
  PHIST_ICHK_IERR(phist_kernels_finalize(&iflag),iflag);
  }
