#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>

#include "scamac.h"

#ifdef _OPENMP
#include <omp.h>
#else
#include "scamac_omp.h"
#endif


/* minimal working example with OpenMP*/

int main(int argc, char *argv[]) {

  ScamacErrorCode err;
  ScamacGenerator * my_gen;
  char *errstr = NULL;

  if (argc>1) {
    // read matrix name & parameters from command line
    err = scamac_parse_argstr(argv[1], &my_gen, &errstr);
    if (err) {
      printf("Problem with example string:\n%s\n",errstr);
      SCAMAC_CHKERR(err);
    }
  } else {
    printf("Error: Expected non-empty argument string\n");
    exit(EXIT_FAILURE);
  }

  err = scamac_generator_check(my_gen, &errstr);
  if (err) {
    printf("*** Problem with example parameters:\n%s\n",errstr);
    if (SCAMAC_DISWARN(err)) {
      SCAMAC_CHKERR(err);
    } else {
      printf("Continue anyways -- you have been warned\n");
    }
  }
  SCAMAC_TRY(scamac_generator_finalize(my_gen));

  printf("Example\n-------\n%s\n",scamac_generator_query_name(my_gen));

  char * data_s;
  SCAMAC_TRY(scamac_generator_parameter_desc(my_gen, "desc", &data_s));
  printf("\nParameters\n----------\n%s\n\n",data_s);
  free(data_s);

  ScamacIdx nrow = scamac_generator_query_nrow(my_gen);
  ScamacIdx maxnzrow = scamac_generator_query_maxnzrow(my_gen);

  ScamacIdx tot_nrow=0;
  ScamacIdx n_nz=0;
  double elapsedt=0.0;


  #pragma omp parallel
  {
    double t1 = omp_get_wtime();

    ScamacWorkspace * my_ws;
    SCAMAC_TRY(scamac_workspace_alloc(my_gen, &my_ws));


    ScamacIdx *cind = malloc(maxnzrow * sizeof *cind);
    double *val;
    if (scamac_generator_query_valtype(my_gen) == SCAMAC_VAL_REAL) {
      val = malloc(maxnzrow * sizeof *val);
    } else {// SCAMAC_VAL_COMPLEX
      val = malloc(2*maxnzrow * sizeof *val);
    }

    ScamacIdx idx;
    ScamacIdx my_nrow = 0;
    ScamacIdx my_nz = 0;

    #pragma omp for schedule(dynamic)
    for (idx=0; idx<nrow; idx++) {
      ScamacIdx k;
      SCAMAC_TRY(scamac_generate_row(my_gen, my_ws, idx, SCAMAC_DEFAULT, &k, cind, val));
      my_nz += k;
      my_nrow++;
    }

    free(cind);
    free(val);
    SCAMAC_TRY(scamac_workspace_free(my_ws));

    double t2 = omp_get_wtime();

    #pragma omp critical
    {
      printf("========================================\nThread number: %d [of %d threads]\n generated %"SCAMACPRIDX" rows with %"SCAMACPRIDX" non-zeros\n in %12.3f seconds.\n",
             omp_get_thread_num(),omp_get_num_threads(), my_nrow, my_nz, t2-t1);

      n_nz += my_nz;
      tot_nrow += my_nrow;
      elapsedt = fmax(elapsedt,t2-t1);
    }

  }

  printf("\n========================================\nGenerated a total of %"SCAMACPRIDX" rows (=%"SCAMACPRIDX") with %"SCAMACPRIDX" non-zeros\n in %12.3f seconds [%d rows/second]\n",
         tot_nrow, nrow, n_nz, elapsedt, (int) ((double) nrow/elapsedt));


  SCAMAC_TRY(scamac_generator_destroy(my_gen));

  return 0;
}
