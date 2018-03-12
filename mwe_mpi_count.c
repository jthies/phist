#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>
#include <mpi.h>

/* minimal working example */

#include "scamac.h"

// split integer range [a...b-1] in n nearly equally sized pieces [ia...ib-1], for i=0,...,n-1
void split_range(ScamacIdx a, ScamacIdx b, ScamacIdx n, ScamacIdx i, ScamacIdx *ia, ScamacIdx *ib) {
  ScamacIdx m = (b-a-1)/n + 1;
  ScamacIdx d = n-(n*m -(b-a));
  if (i < d) {
    *ia = m*i + a;
    *ib = m*(i+1) + a;
  } else {
    *ia = m*d + (i-d)*(m-1) + a;
    *ib = m*d + (i-d+1)*(m-1) + a;
  }
}

int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv);

  int mpi_world_size, mpi_rank;

  MPI_Comm_size(MPI_COMM_WORLD, &mpi_world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  
  ScamacErrorCode err;
  ScamacGenerator * my_gen;
    
  if (mpi_rank == 0) {
    // read matrix name & parameters from command line
    if (argc>1) {// assume, that MPI_Init has cleared mpirun parameters
      char *errstr = NULL;
      err=scamac_parse_argstr(argv[1], &my_gen, &errstr);
      if (err) {
        printf("Problem with example string:\n%s\n",errstr);
        SCAMAC_CHKERR_MPI(err);
      }

      err = scamac_generator_check(my_gen, &errstr);
      if (err) {
        printf("*** Problem with example parameters:\n%s\n",errstr);
        if (SCAMAC_DISWARN(err)) {//real error
          SCAMAC_CHKERR_MPI(err);
        } else {
          printf("Continue anyways -- you have been warned\n");
        }
      }
      SCAMAC_TRY_MPI(scamac_generator_finalize(my_gen));
      
      // print some information about the example
      printf("Example\n-------\n%s\n",scamac_generator_query_name(my_gen));

      char * data_s=NULL;
      SCAMAC_TRY_MPI(scamac_generator_parameter_desc(my_gen, "desc", &data_s));
      printf("\nParameters\n----------\n%s\n\n",data_s);
      free(data_s);  
      
    } else {
      printf("Error: Expected non-empty argument string\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
     
    // Everything's fine with the example. Now tell the other processes what to do.
    char * argstr;
    // Obtain argument string ...
    SCAMAC_TRY_MPI(scamac_generator_parameter_desc(my_gen, "argstr", &argstr));
    // ... and pass to other processes.
    int strsize = strlen(argstr)+1;
    MPI_Bcast(&strsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(argstr, strsize * sizeof *argstr, MPI_CHAR, 0, MPI_COMM_WORLD);
    free(argstr);
    
  } else { // receive argstr
    
    // Receive argument string from root process.    
    int strsize;
    MPI_Bcast(&strsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    char * argstr = malloc( strsize * sizeof *argstr);
    MPI_Bcast(argstr, strsize * sizeof *argstr, MPI_CHAR, 0, MPI_COMM_WORLD);
    
    // and obtain generator, just as for rank=0. No errors should occur now.
    SCAMAC_TRY_MPI(scamac_parse_argstr(argstr, &my_gen, NULL));
    SCAMAC_TRY_MPI(scamac_generator_finalize(my_gen));
    
    free(argstr);
  }
  
  
  // now each process has a finalized generator
  
  ScamacIdx nrow = scamac_generator_query_nrow(my_gen);
  ScamacIdx maxnzrow = scamac_generator_query_maxnzrow(my_gen);

  double t1 = MPI_Wtime();
    
  ScamacWorkspace * my_ws;
  SCAMAC_TRY_MPI(scamac_workspace_alloc(my_gen, &my_ws));
  
  ScamacIdx *cind = malloc(maxnzrow * sizeof *cind);
  double *val;
  if (scamac_generator_query_valtype(my_gen) == SCAMAC_VAL_REAL) {
    val = malloc(maxnzrow * sizeof *val);
  } else {// SCAMAC_VAL_COMPLEX
    val = malloc(2*maxnzrow * sizeof *val);
  }

  ScamacIdx idx;
  ScamacIdx ia,ib;
  // current process will generate rows ia ... ib-1
  split_range(0,nrow, mpi_world_size, mpi_rank, &ia, &ib);

  ScamacIdx my_nrow = ib-ia;
  ScamacIdx my_nz = 0;

  for (idx=ia; idx<ib; idx++) {
    ScamacIdx k;
    SCAMAC_TRY_MPI(scamac_generate_row(my_gen, my_ws, idx, SCAMAC_DEFAULT, &k, cind, val));
    my_nz += k;
  }

  free(cind);
  free(val);
  SCAMAC_TRY_MPI(scamac_workspace_free(my_ws));
  SCAMAC_TRY_MPI(scamac_generator_destroy(my_gen));

  double t2 = MPI_Wtime();
  
  double elapsedt = t2-t1;
    
  ScamacIdx * all_nrow = NULL;
  ScamacIdx * all_nz = NULL;
  double * all_dt = NULL;
  

  if (mpi_rank==0) {
    printf("Gather results ...\n");
    all_nrow = malloc(mpi_world_size * sizeof * all_nrow);        
    all_nz = malloc(mpi_world_size * sizeof * all_nz);
    all_dt = malloc(mpi_world_size * sizeof * all_dt);
  }
#ifdef SCAMAC_INDEX_TYPE_int64
  MPI_Gather( &my_nrow,  1, MPI_INT64_T, all_nrow, 1, MPI_INT64_T, 0, MPI_COMM_WORLD);
  MPI_Gather( &my_nz,    1, MPI_INT64_T, all_nz,   1, MPI_INT64_T, 0, MPI_COMM_WORLD);
#elif defined SCAMAC_INDEX_TYPE_int32
  MPI_Gather( &my_nrow,  1, MPI_INT32_T, all_nrow, 1, MPI_INT32_T, 0, MPI_COMM_WORLD);
  MPI_Gather( &my_nz,    1, MPI_INT32_T, all_nz,   1, MPI_INT32_T, 0, MPI_COMM_WORLD);
#endif
  MPI_Gather( &elapsedt, 1, MPI_DOUBLE,   all_dt,   1, MPI_DOUBLE,   0, MPI_COMM_WORLD);
  
  if (mpi_rank==0) {
    ScamacIdx tot_nz = 0;
    ScamacIdx tot_nrow = 0;
    double max_dt = 0.0;
        
    int irank;
    for (irank=0;irank<mpi_world_size;irank++) {
      tot_nrow += all_nrow[irank];
      tot_nz += all_nz[irank];
      max_dt = fmax(max_dt,all_dt[irank]);
      printf("========================================\nProcess number: %d [of %d processes]\n generated %"SCAMACPRIDX" rows with %"SCAMACPRIDX" non-zeros\n in %12.3f seconds.\n",
            irank, mpi_world_size, all_nrow[irank], all_nz[irank], all_dt[irank]);
    }
             
    printf("\n========================================\nGenerated a total of %"SCAMACPRIDX" rows with %"SCAMACPRIDX" non-zeros\n in %12.3f seconds [%d rows/second]\n",
         tot_nrow, tot_nz, max_dt, (int) ((double) nrow/max_dt));
         
    
      
  }

  MPI_Finalize();

  return 0;
}
