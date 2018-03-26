#ifndef SCAMAC_MATRIX_ANDERSON_H
#define SCAMAC_MATRIX_ANDERSON_H

#include "scamac_include.h"
#include "scamac_sparserow.h"
#include "scamac_rng.h"

/* >>>>
 * symmetric real
 * Anderson model of localization in 1D, 2D, 3D
 * <<<< */

typedef enum {Anderson_simple, Anderson_backforth} scamac_option_Anderson_sweep_ty;
typedef enum {Anderson_open, Anderson_periodic} scamac_option_Anderson_bc_ty;

typedef struct {
  /* dimensions of cuboid: Lx */
  /* = 5 */
  int Lx;
  /* dimensions of cuboid: Lx */
  /* = 5 */
  int Ly;
  /* dimensions of cuboid: Lx */
  /* = 5 */
  int Lz;
  /* hopping strength */
  /* = 1.0 */
  double t;
  /* random on-site potential [-ranpot, ranpot] */
  /* = 0.0 */
  double ranpot;
  // open or periodic boundary conditions
  scamac_option_Anderson_bc_ty boundary_conditions;
  /* random seed */
  /* = 1 */
  scamac_rng_seed_ty seed;
  /* mode of traversal of cuboid */
  /* = Anderson_simple */
  scamac_option_Anderson_sweep_ty sweep;
} scamac_matrix_Anderson_params_st;

typedef struct {
  // workspace variables
  scamac_ransrc_st * rng;
} scamac_matrix_Anderson_work_st;

ScamacErrorCode scamac_matrix_Anderson_check        (const scamac_matrix_Anderson_params_st * par, char ** desc);
ScamacErrorCode scamac_matrix_Anderson_tables_create(const scamac_matrix_Anderson_params_st * par, void ** tab, ScamacInfo * info);
ScamacErrorCode scamac_matrix_Anderson_work_alloc   (const scamac_matrix_Anderson_params_st * par, const void * tab, scamac_matrix_Anderson_work_st ** ws);
ScamacErrorCode scamac_matrix_Anderson_work_free    (scamac_matrix_Anderson_work_st * ws);
ScamacErrorCode scamac_matrix_Anderson_generate_row (const scamac_matrix_Anderson_params_st * par, const void * tab, scamac_matrix_Anderson_work_st * ws, ScamacIdx irow, ScamacFlag flag, scamac_sparserow_real_st * row);

#endif /* SCAMAC_MATRIX_ANDERSON_H */
