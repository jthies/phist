#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#else
#include "scamac_omp.h"
#endif

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

  char *examplestr = argv[1];

  ScamacErrorCode err;
  ScamacGenerator * my_gen;
  char *errstr = NULL;

  err = scamac_parse_argstr(examplestr, &my_gen, &errstr);
  if (err) {
    printf("Problem with example string:\n%s\n",errstr);
    exit(EXIT_FAILURE);
  }

  err = scamac_generator_check(my_gen, &errstr);
  if (err) {
    printf("Problem with example parameters:\n%s\n",errstr);
    if (err != SCAMAC_EWARNING) {
      exit(EXIT_FAILURE);
    } else {
      printf("Continue anyways -- you have been warned\n");
    }
  }
  SCAMAC_TRY(scamac_generator_finalize(my_gen));

  printf("Example\n-------\n%s\n",scamac_generator_query_name(my_gen));

  char * data_s;
  SCAMAC_TRY(scamac_generator_parameter_desc(my_gen, &data_s));

  printf("\nParameters\n----------\n%s\n\n",data_s);

  ScamacIdx nrow = scamac_generator_query_nrow(my_gen);
  ScamacIdx maxnzrow = scamac_generator_query_maxnzrow(my_gen);

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
    ScamacIdx ia,ib;
    // current thread will generate rows ia ... ib-1
    split_range(0,nrow, omp_get_num_threads(), omp_get_thread_num(), &ia, &ib);

    ScamacIdx my_nz = 0;

    for (idx=ia; idx<ib; idx++) {
      ScamacIdx k;
      SCAMAC_TRY(scamac_generate_row(my_gen, my_ws, idx, SCAMAC_DEFAULT, &k, cind, val));
      my_nz += k;
    }

    free(cind);
    free(val);
    SCAMAC_TRY(scamac_workspace_free(my_ws));

    double t2 = omp_get_wtime();

    #pragma omp critical
    {
      printf("========================================\nThread number: %d [of %d threads]\n generated %"SCAMACPRIDX" rows with %"SCAMACPRIDX" non-zeros\n in %12.3f seconds.\n",
             omp_get_thread_num(),omp_get_num_threads(), ib-ia, my_nz, t2-t1);

      n_nz += my_nz;
      elapsedt = fmax(elapsedt,t2-t1);
    }

  }

  printf("\n========================================\nGenerated a total of %"SCAMACPRIDX" rows with %"SCAMACPRIDX" non-zeros\n in %12.3f seconds [%d rows/second]\n",
         nrow, n_nz, elapsedt, (int) ((double) nrow/elapsedt));


  SCAMAC_TRY(scamac_generator_destroy(my_gen));

  return 0;
}
