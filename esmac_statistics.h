#ifndef ESMAC_STATISTICS_H
#define ESMAC_STATISTICS_H

#include <stdio.h>
#include <stdbool.h>

#include "esmac_generator.h"

typedef struct {
  // matrix data
  esmac_idx_t nr,nc;
  // pattern statistics
  esmac_idx_t n_nz;
  esmac_idx_t n_nz_left, n_nz_right;
  esmac_idx_t n_nz_row_min, n_nz_row_max;
  esmac_idx_t n_nz_row_max_left, n_nz_row_max_right; // only max - min always = 0
  esmac_idx_t n_zero_row, n_zero_diag;
  // bandwidth
  esmac_idx_t bw_left,bw_right,bw;

  // value statistics
  double vmin,vmax;
  double veps;
  double vmin_diag,vmax_diag;
  double veps_diag;

  // spectral statistics
  double gershgorin_min, gershgorin_max;
  
} esmac_matrix_statistics_t;

void esmac_print_statistics(FILE *f, const esmac_matrix_statistics_t *st);

//void esmac_compute_matrix_statistics(const esmac_generator_t * gen, esmac_matrix_statistics_t *st);
void esmac_collect_matrix_statistics(const esmac_generator_t * gen, esmac_matrix_statistics_t *st, int pattern_px, int *pattern, bool show_progress);

#endif /* ESMAC_STATISTICS_H */
