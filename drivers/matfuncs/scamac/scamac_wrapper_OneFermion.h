/** \file
 *  \author Andreas Alvermann (<alvermann@physik.uni-greifswald.de>)
 *  \date   October 2017
 *  \brief  ...
 *  \ingroup matrix
 */

#ifndef SCAMAC_WRAPPER_ONEFERMION_H
#define SCAMAC_WRAPPER_ONEFERMION_H

#include <stdbool.h>

#include "scamac_include.h"

#include "scamac_matrix_FreeFermionChain.h"

/* >>>>
 * symmetric real
 * one fermion on a chain
 * <<<< */

typedef struct {
  /* hopping strength */
  /* = 1.0 */
  double t;
  /* number of sites */
  /* = 10 */
  int n_sites;
  /* open (false) or periodic (true) boundary conditions */
  /* = true */
  bool PBC;
} scamac_wrapper_OneFermion_params_st;


ScamacErrorCode scamac_wrapper_OneFermion_unwrap(const scamac_wrapper_OneFermion_params_st * par_in, scamac_matrix_FreeFermionChain_params_st * par_out);
ScamacErrorCode scamac_wrapper_OneFermion_check(const scamac_wrapper_OneFermion_params_st * par, char ** desc);


#endif /* SCAMAC_WRAPPER_ONEFERMION_H */
