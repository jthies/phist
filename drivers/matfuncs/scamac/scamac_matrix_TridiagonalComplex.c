#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "scamac_internal.h"
#include "scamac_string.h"
#include "scamac_safeint.h"

#include "scamac_matrix_TridiagonalComplex.h"

ScamacErrorCode scamac_matrix_TridiagonalComplex_check(const scamac_matrix_TridiagonalComplex_params_st * par, char ** desc) {
  ScamacErrorCode err = SCAMAC_EOK;
  scamac_string_st str;
  if (desc) {
    scamac_string_empty(&str);
  }
  SCAMAC_DESC_ERR(par->n<= 0, "n <= 0");
  if (desc) {
    *desc = scamac_string_get(&str);
  }
  return err;
}

ScamacErrorCode scamac_matrix_TridiagonalComplex_tables_create(const scamac_matrix_TridiagonalComplex_params_st * par, void ** tab, ScamacInfo * info) {
  if (info) {
    info->nrow = par->n;
    info->ncol = par->n;
    info->maxnzrow = 3;
    info->maxnzcol = 3;
    info->maxnz    = scamac_safe_mult(info->nrow, info->maxnzrow);
    if (info->maxnz < 0) {
      return SCAMAC_EOVERFLOW;
    }
    info->valtype=SCAMAC_VAL_COMPLEX;
    info->symmetry=SCAMAC_GENERAL;
  }
  return SCAMAC_EOK;
}

ScamacErrorCode scamac_matrix_TridiagonalComplex_generate_row(const scamac_matrix_TridiagonalComplex_params_st * par, const void * tab, void * ws, ScamacIdx irow, ScamacFlag flag, scamac_sparserow_cplx_st * row) {

  if ( !par || !row) {
    return SCAMAC_ENULL;
  }

  if (flag & ~SCAMAC_TRANSPOSE & ~SCAMAC_CONJUGATE & ~SCAMAC_KEEPZEROS) {
    return SCAMAC_EINVALID;
  }
  bool fl_transpose = (flag & SCAMAC_TRANSPOSE) != 0;
  bool fl_conjugate = (flag & SCAMAC_CONJUGATE) != 0;


  if ( (irow<0) || (irow >= par->n) ) {
    return SCAMAC_ERANGE;
  }

  double complex dd;
  if (!fl_conjugate) {
    dd = par->diag_re + I * par->diag_im;
  } else {
    dd = par->diag_re - I * par->diag_im;
  }

  double psub, psuper;
  if (!fl_transpose) {
    psub   = par->subdiag;
    psuper = par->supdiag;
  } else {
    psub   = par->supdiag;
    psuper = par->subdiag;
  }

  if (irow == 0) {
    scamac_sparserow_cplx_add(row, dd,           0     );
    scamac_sparserow_cplx_add(row, psuper,       1     );
  } else if (0 < irow && irow < (par->n-1)) {
    scamac_sparserow_cplx_add(row, psub,         irow-1);
    scamac_sparserow_cplx_add(row, dd,           irow  );
    scamac_sparserow_cplx_add(row, psuper,       irow+1);
  } else if (irow == (par->n-1)) {
    scamac_sparserow_cplx_add(row, psub,         irow-1);
    scamac_sparserow_cplx_add(row, dd,           irow  );
  } else {
    return SCAMAC_ERANGE;
  }

  return SCAMAC_EOK;
}
