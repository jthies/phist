#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>

/* minimal working example */

#include "scamac.h"

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
    // choose an example, just for demonstration
    SCAMAC_TRY(scamac_generator_obtain("Hubbard", &my_gen));
    SCAMAC_TRY(scamac_generator_set_int(my_gen, "n_sites", 20));
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

  ScamacIdx nrow = scamac_generator_query_nrow(my_gen);
  ScamacIdx ncol = scamac_generator_query_ncol(my_gen);

  printf("Example has %"SCAMACPRIDX" rows x %"SCAMACPRIDX" columns\n",nrow,ncol);

  SCAMAC_TRY(scamac_generator_destroy(my_gen));

  return 0;
}
