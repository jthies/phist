#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>
#include <unistd.h>

/* minimal working example */

#include "scamac.h"
#include "scamac_tools.h"


void print_progress_bar(double progress) {
  if (isatty(fileno(stdout))) {
    if (progress<=1.0) {
      int barWidth = 70;

      fprintf(stdout,"[");
      int pos = barWidth * progress;
      int i;
      for (i = 0; i < barWidth; ++i) {
        if (i < pos) fprintf(stdout,"=");
        else if (i == pos) fprintf(stdout,">");
        else fprintf(stdout," ");
      }
      fprintf(stdout,"] %3d%%\r", (int) (100.0 * progress) );
      fflush(stdout);

    }
    if (progress==1.0) {
      fprintf(stdout,"\n");
    }
  } else {
    fprintf(stdout,"[ progress: %3d%% ]\n", (int) (100.0 * progress) );
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
  free(data_s);

  ScamacIdx nrow, n_nz=0;
  nrow = scamac_generator_query_nrow(my_gen);

  scamac_matrix_pattern_st * pat;
  SCAMAC_TRY(scamac_pattern_alloc(1024,1024, &pat));
  SCAMAC_TRY(scamac_pattern_empty(pat, scamac_generator_query_nrow(my_gen), scamac_generator_query_ncol(my_gen), scamac_generator_query_valtype(my_gen)));


  ScamacWorkspace * my_ws;
  SCAMAC_TRY(scamac_workspace_alloc(my_gen, &my_ws));

  ScamacIdx *cind;
  double *val;
  SCAMAC_TRY(scamac_alloc_cind_val(my_gen, SCAMAC_DEFAULT, &cind, &val));

  ScamacIdx idx;

  double t1 = (double) clock() / (double) CLOCKS_PER_SEC;


  for (idx=0; idx<nrow; idx++) {
    ScamacIdx k;
    SCAMAC_TRY(scamac_generate_row(my_gen, my_ws, idx, SCAMAC_DEFAULT, &k, cind, val));
    SCAMAC_TRY(scamac_pattern_update(pat, idx, k, cind));
    n_nz += k;
    if ( idx % 10000 == 0) {
      print_progress_bar((double) idx/nrow);
    }
  }

  print_progress_bar(1.0);

  free(cind);
  free(val);
  SCAMAC_TRY(scamac_workspace_free(my_ws));

  double t2 = (double) clock() / (double) CLOCKS_PER_SEC;

  double elapsedt = t2-t1;

  printf("\n========================================\nGenerated a total of %"SCAMACPRIDX" rows with %"SCAMACPRIDX" non-zeros\n in %12.3f seconds [%d rows/second]\n",
         nrow, n_nz, elapsedt, (int) ((double) nrow/elapsedt));


  SCAMAC_TRY(scamac_generator_destroy(my_gen));

  SCAMAC_TRY(scamac_plot_pattern(pat, 5, "mwe.pattern.png"));
  SCAMAC_TRY(scamac_pattern_free(pat));

  return 0;
}
