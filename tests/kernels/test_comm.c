#include "essex_kernels.h"
#include "essex_macros.h"
#include <stdio.h>
#include <assert.h>

#include <mpi.h>

#include "essex_test_helpers.h"

int main(int argc, char** argv)
  {
  int rank, act_rank;
  int size, act_size;
  int ierr;

  comm_ptr_t comm;
  essex_bad_object *bad_obj;

  _ESSEX_TEST_HANDLER_(essex_kernels_init(&argc,&argv,&ierr),ierr,0);
  MPI_Comm_rank(MPI_COMM_WORLD,&act_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&act_size);

  _ESSEX_TEST_HANDLER_(essex_comm_get_rank(bad_obj,&rank,&ierr),ierr,_ESSEX_BAD_CAST_);
  _ESSEX_TEST_HANDLER_(essex_comm_get_size(bad_obj,&size,&ierr),ierr,_ESSEX_BAD_CAST_);

  _ESSEX_TEST_HANDLER_(essex_comm_create(&comm,&ierr),ierr,0);
  _ESSEX_TEST_HANDLER_(essex_comm_get_rank(comm,&rank,&ierr),ierr,0);
  _ESSEX_TEST_HANDLER_(essex_comm_get_size(comm,&size,&ierr),ierr,0);

assert(rank==act_rank);
assert(size==act_size);

fprintf(stderr,"%d (%d)\n",rank,size);

//  essex_comm_delete(&comm);
  _ESSEX_TEST_HANDLER_(essex_kernels_finalize(&ierr),ierr,0);
  }
