#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#else
#include "scamac_omp.h"
#endif

/* minimal working example */

#include "scamac.h"
#include "scamac_tools.h"

// split integer range [a...b-1] in n nearly equal pieces [ia...ib-1], for i=0,...,n-1
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

  ScamacIdx my_nrow = scamac_generator_query_nrow(my_gen);

  scamac_matrix_statistics_st st;
  scamac_statistics_empty(&st, scamac_generator_query_nrow(my_gen), scamac_generator_query_ncol(my_gen), scamac_generator_query_valtype(my_gen));

  ScamacIdx nrow=0, n_nz=0;
  double elapsedt=0.0;

  #pragma omp parallel
  {
    double t1 = omp_get_wtime();

    ScamacWorkspace * my_ws;
    SCAMAC_TRY(scamac_workspace_alloc(my_gen, &my_ws));

    ScamacIdx *cind;
    double *val;
    SCAMAC_TRY(scamac_alloc_cind_val(my_gen, SCAMAC_DEFAULT, &cind, &val));

    //val=NULL;

    ScamacIdx idx;
    ScamacIdx ia,ib;
    // current thread will generate rows ia ... ib-1
    split_range(0,my_nrow, omp_get_num_threads(), omp_get_thread_num(), &ia, &ib);

    ScamacIdx my_nz = 0;

    scamac_matrix_statistics_st st_thread;

    SCAMAC_TRY(scamac_statistics_empty(&st_thread, scamac_generator_query_nrow(my_gen), scamac_generator_query_ncol(my_gen), scamac_generator_query_valtype(my_gen)));

    for (idx=ia; idx<ib; idx++) {
      ScamacIdx k;
      SCAMAC_TRY(scamac_generate_row(my_gen, my_ws, idx, SCAMAC_DEFAULT, &k, cind, val));
      SCAMAC_TRY(scamac_statistics_update(&st_thread, idx, k, cind, val));
      my_nz += k;
      if (omp_get_thread_num() == 0) {
        if ( (idx-ia) % 50000 == 0) {
          print_progress_bar((double) (idx-ia)/(ib-ia));
        }
      }
    }

    if (omp_get_thread_num() == 0) {
      print_progress_bar(1.0);
    }


    free(cind);
    free(val);
    SCAMAC_TRY(scamac_workspace_free(my_ws));

    double t2 = omp_get_wtime();

    #pragma omp barrier

    #pragma omp critical
    {
      printf("========================================\nThread number: %d [of %d threads]\n generated %"SCAMACPRIDX" rows with %"SCAMACPRIDX" non-zeros\n in %12.3f seconds.\n",
             omp_get_thread_num(),omp_get_num_threads(), ib-ia, my_nz, t2-t1);

      nrow += ib-ia;
      n_nz += my_nz;
      SCAMAC_TRY(scamac_statistics_combine(&st, &st_thread));
      elapsedt = fmax(elapsedt,t2-t1);
    }

  }

  printf("\n========================================\nGenerated a total of %"SCAMACPRIDX" rows with %"SCAMACPRIDX" non-zeros\n in %12.3f seconds [%d rows/second]\n",
         nrow, n_nz, elapsedt, (int) ((double) nrow/elapsedt));


  SCAMAC_TRY(scamac_generator_destroy(my_gen));

  char *desc;
  SCAMAC_TRY(scamac_statistics_print(&st, &desc));
  printf("\n======== matrix statistics ========\n\n");
  puts(desc);
  free(desc);

  return 0;
}
