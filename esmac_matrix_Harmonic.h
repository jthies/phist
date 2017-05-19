#ifndef ESMAC_MATRIX_HARMONIC_H
#define ESMAC_MATRIX_HARMONIC_H

#include "esmac_types.h"

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
} esmac_matrix_Harmonic_params_t;


typedef struct {
  esmac_idx_t ns;
  int maxrowlength;
} esmac_matrix_Harmonic_work_t;


esmac_matrix_Harmonic_work_t * esmac_matrix_Harmonic_alloc(const esmac_matrix_Harmonic_params_t * par, int * info);
void esmac_matrix_Harmonic_free(esmac_matrix_Harmonic_work_t *ws);
esmac_idx_t esmac_matrix_Harmonic_ns(const esmac_matrix_Harmonic_work_t *ws);
//int esmac_matrix_Harmonic_maxrowlength(const esmac_matrix_Harmonic_work_t *ws);
int esmac_matrix_Harmonic_set_info(const esmac_matrix_Harmonic_params_t * par, esmac_matrix_Harmonic_work_t *ws, esmac_matrix_info_t *info);

int esmac_matrix_Harmonic_row(const esmac_matrix_Harmonic_params_t * par, esmac_matrix_Harmonic_work_t *ws,
				 esmac_idx_t irow, esmac_idx_t *cind, double *val);

#endif /* ESMAC_MATRIX_HARMONIC_H */
