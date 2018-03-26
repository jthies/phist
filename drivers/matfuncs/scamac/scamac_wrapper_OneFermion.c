#include "scamac_internal.h"
#include "scamac_string.h"
#include "scamac_wrapper_OneFermion.h"

ScamacErrorCode scamac_wrapper_OneFermion_check(const scamac_wrapper_OneFermion_params_st * par, char ** desc) {
  ScamacErrorCode err = SCAMAC_EOK;
  scamac_string_st str;
  if (desc) {
    scamac_string_empty(&str);
  }
  // add conditions as:
  // SCAMAC_DESC_ERR(par->... <= 0,    "... <= 0");
  SCAMAC_DESC_ERR(par->n_sites<1,"n_sites<1");
  SCAMAC_DESC_WARN(par->n_sites>100,"parameter n_sites is very large");
  if (desc) {
    *desc = scamac_string_get(&str);
  }
  return err;
}

ScamacErrorCode scamac_wrapper_OneFermion_unwrap(const scamac_wrapper_OneFermion_params_st * par_in, scamac_matrix_FreeFermionChain_params_st * par_out) {
  par_out -> t = par_in -> t;
  par_out -> n_fermions = 1;
  par_out -> n_species = 1;
  par_out -> n_sites = par_in -> n_sites;
  par_out -> PBC = par_in -> PBC;
  return SCAMAC_EOK;
}


