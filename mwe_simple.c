#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>

/* minimal working example */

#include "scamac.h"

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

  double t1 = (double) clock() / (double) CLOCKS_PER_SEC;

  ScamacWorkspace * my_ws;
  SCAMAC_TRY(scamac_workspace_alloc(my_gen, &my_ws));

  ScamacIdx *cind = malloc(maxnzrow * sizeof *cind);
  double *val;
  if (scamac_generator_query_valtype(my_gen) == SCAMAC_VAL_REAL) {
    val = malloc(maxnzrow * sizeof *val);
  } else {// SCAMAC_VAL_COMPLEX
    val = malloc(2*maxnzrow * sizeof *val);
  }

  ScamacIdx n_nz=0;
  ScamacIdx idx;
  for (idx=0; idx<nrow; idx++) {
    ScamacIdx k;
    SCAMAC_TRY(scamac_generate_row(my_gen, my_ws, idx, SCAMAC_DEFAULT, &k, cind, val));
    n_nz += k;
  }

  double t2 = (double) clock() / (double) CLOCKS_PER_SEC;
  double elapsedt = t2-t1;

  printf("\n========================================\nGenerated a total of %"SCAMACPRIDX" rows with %"SCAMACPRIDX" non-zeros\n in %12.3f seconds [%d rows/second]\n",
         nrow, n_nz, elapsedt, (int) ((double) nrow/elapsedt));

  free(cind);
  free(val);
  SCAMAC_TRY(scamac_workspace_free(my_ws));
  SCAMAC_TRY(scamac_generator_destroy(my_gen));


  return 0;
}
