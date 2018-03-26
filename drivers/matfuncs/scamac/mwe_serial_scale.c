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

  int n_sites, n_fermions;

  //  SCAMAC_TRY(scamac_parse_argstr("Hubbard", &my_gen, &errstr));
  SCAMAC_TRY(scamac_generator_obtain("Hubbard", &my_gen));

  printf("n_sites n_fermions   nrows x ncols\n------------------------------\n");
  
  n_sites=2;
  n_fermions=n_sites/2;

  while (true) {
    
    SCAMAC_TRY(scamac_generator_set_int(my_gen, "n_sites", n_sites));
    SCAMAC_TRY(scamac_generator_set_int(my_gen, "n_fermions", n_fermions));
  
    err = scamac_generator_check(my_gen, &errstr);
    if (err) {
      break;
    }
    err = scamac_generator_finalize(my_gen);
    if (err) {
      break;
    }

    ScamacIdx nrow = scamac_generator_query_nrow(my_gen);
    ScamacIdx ncol = scamac_generator_query_ncol(my_gen);

    printf("%5d %5d       %"SCAMACPRIDX" x %"SCAMACPRIDX"\n",n_sites,n_fermions, nrow,ncol);

    n_sites++;
    n_fermions=n_sites/2;

  }

  printf("\nTerminated with error: %s\n",scamac_error_desc(err));

  SCAMAC_TRY(scamac_generator_destroy(my_gen));

  return 0;
}
