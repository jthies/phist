/** \file
 *  \author Andreas Alvermann (<alvermann@physik.uni-greifswald.de>)
 *  \date   October 2017
 *  \brief  ...
 *  \ingroup matrix
 */

#ifndef SCAMAC_MATRIX_HARMONIC_H
#define SCAMAC_MATRIX_HARMONIC_H

#include "scamac_include.h"
#include "scamac_sparserow.h"

/* >>>>
 * symmetric real
 * quantum Harmonic oscillator
 * <<<< */

typedef struct {
  // oscillator frequency
  // = 1.0
  double omega;
  // oscillator shift
  // = 0.0
  double lambda;
  // number of bosons
  // = 100
  int n_bos;
} scamac_matrix_Harmonic_params_st;


ScamacErrorCode scamac_matrix_Harmonic_check(const scamac_matrix_Harmonic_params_st * par, char ** desc);
ScamacErrorCode scamac_matrix_Harmonic_tables_create(const scamac_matrix_Harmonic_params_st * par, void ** tab, ScamacInfo * info);
ScamacErrorCode scamac_matrix_Harmonic_generate_row(const scamac_matrix_Harmonic_params_st * par, const void * tab, void * ws, ScamacIdx irow, ScamacFlag flag, scamac_sparserow_real_st * row);

#endif /* SCAMAC_MATRIX_HARMONIC_H */
