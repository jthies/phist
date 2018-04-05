/** \file
 *  \author Andreas Alvermann (<alvermann@physik.uni-greifswald.de>)
 *  \date   October 2017
 *  \brief  ...
 *  \ingroup matrix
 */

#ifndef SCAMAC_MATRIX_TRIDIAGONAL_H
#define SCAMAC_MATRIX_TRIDIAGONAL_H

#include "scamac_include.h"
#include "scamac_sparserow.h"

/* >>>>
 * hermitian complex
 * Tridiagonal matrix
 * <<<< */

typedef struct {
// matrix dimension
  // = 100
  int n;
  // diagonal element
  // = 0.0
  double diag;
  // off-diagonal element
  // = 1.0
  double offdiag;
  // phase for off-diagonal element
  // = 0.0
  double phi;
} scamac_matrix_Tridiagonal_params_st;


ScamacErrorCode scamac_matrix_Tridiagonal_check(const scamac_matrix_Tridiagonal_params_st * par, char ** desc);
ScamacErrorCode scamac_matrix_Tridiagonal_tables_create(const scamac_matrix_Tridiagonal_params_st * par, void ** tab, ScamacInfo * info);
ScamacErrorCode scamac_matrix_Tridiagonal_generate_row(const scamac_matrix_Tridiagonal_params_st * par, const void * tab, void * ws, ScamacIdx irow, ScamacFlag flag, scamac_sparserow_cplx_st * row);

#endif /* SCAMAC_MATRIX_TRIDIAGONAL_H */
