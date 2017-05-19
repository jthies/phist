#include <stdlib.h>
#include <math.h>

#include "esmac_constants.h"
#include "esmac_matrix_Harmonic.h"

esmac_matrix_Harmonic_work_t * esmac_matrix_Harmonic_alloc(const esmac_matrix_Harmonic_params_t * par, int * info) {
  esmac_matrix_Harmonic_work_t * ws = malloc(sizeof *ws);
  ws->ns = par->n_bos;
  ws->maxrowlength=3;
  if (info) {*info=0;}
  return ws;
}

void esmac_matrix_Harmonic_free(esmac_matrix_Harmonic_work_t *ws) {
  if (ws) {
    free(ws);
  }
}

/*
esmac_idx_t esmac_matrix_Harmonic_ns(const esmac_matrix_Harmonic_work_t *ws) {
  return ws->ns;
}

int esmac_matrix_Harmonic_maxrowlength(const esmac_matrix_Harmonic_work_t *ws) {
  return ws->maxrowlength;
}
*/

int esmac_matrix_Harmonic_set_info(const esmac_matrix_Harmonic_params_t * par, esmac_matrix_Harmonic_work_t *ws, 
                                    esmac_matrix_info_t *info) {
  info->nrow = par->n_bos;
  info->ncol = par->n_bos;
  info->maxrowlength = 3;
  info->valtype = ESMAC_VAL_REAL;
  info->symmetry = ESMAC_SYMMETRIC;
  return 0;
}


int esmac_matrix_Harmonic_row(const esmac_matrix_Harmonic_params_t * par, esmac_matrix_Harmonic_work_t *ws,
				 esmac_idx_t irow, esmac_idx_t *cind, double *val) {
  if (irow == 0) {
    cind[0]=0;
    cind[1]=1;
    val[0]=0.0;
    val[1]=par->lambda;
    return 2;
  } else if (0 < irow && irow < (par->n_bos-1)) {
    cind[0]=irow-1;
    cind[1]=irow;
    cind[2]=irow+1;
    val[0]=par->lambda * sqrt(1.0 * irow);
    val[1]=par->omega * irow;
    val[2]=par->lambda * sqrt(1.0 * (irow+1));
    return 3;
  } else if (irow == (par->n_bos - 1)) {
    cind[0]=irow-1;
    cind[1]=irow;
    val[0]=par->lambda * sqrt(1.0 * irow);
    val[1]=par->omega * irow;
    return 2;
  } else {
    return 0;
  }
}

