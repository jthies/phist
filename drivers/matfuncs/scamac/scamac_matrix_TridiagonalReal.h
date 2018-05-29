/** \file
 *  \author Andreas Alvermann (<alvermann@physik.uni-greifswald.de>)
 *  \date   October 2017
 *  \brief  ...
 *  \ingroup matrix
 */

#ifndef SCAMAC_MATRIX_TRIDIAGONALREAL_H
#define SCAMAC_MATRIX_TRIDIAGONALREAL_H

#include "scamac_include.h"
#include "scamac_sparserow.h"

/* >>>>
 * general real
 * Non-symmetric tridiagonal matrix
 * <<<< */

typedef struct {
// matrix dimension
  // = 100
  int n;
  // diagonal element
  // = 0.0
  double diag;
  // off-diagonal element (below diagonal)
  // = 1.0
  double supdiag;
  // off-diagonal element (above diagonal)
  // = 1.0
  double subdiag;
} scamac_matrix_TridiagonalReal_params_st;


ScamacErrorCode scamac_matrix_TridiagonalReal_check(const scamac_matrix_TridiagonalReal_params_st * par, char ** desc);
ScamacErrorCode scamac_matrix_TridiagonalReal_tables_create(const scamac_matrix_TridiagonalReal_params_st * par, void ** tab, ScamacInfo * info);
ScamacErrorCode scamac_matrix_TridiagonalReal_generate_row(const scamac_matrix_TridiagonalReal_params_st * par, const void * tab, void * ws, ScamacIdx irow, ScamacFlag flag, scamac_sparserow_real_st * row);

#endif /* SCAMAC_MATRIX_TRIDIAGONALREAL_H */
