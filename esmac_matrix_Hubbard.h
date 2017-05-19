#ifndef ESMAC_MATRIX_HUBBARD_H
#define ESMAC_MATRIX_HUBBARD_H

#include "esmac_types.h"
#include "esmac_multidx.h"
#include "esmac_sparserow.h"
#include "esmac_dof_fermions.h"
#include "esmac_rng.h"

typedef struct {
  /* hopping strength */
  /* = 1.0 */
  double t;
  /* Hubbard interaction */
  /* = 0.0 */
  double U;
  /* number of sites */
  /* = 10 */
  int n_sites;
  /* number of fermions per spin orientation */
  /* = 5 */
  int n_fermions;
  // ESMAC_OBC or ESMAC_PBC
  int boundary_conditions; 
  /* random on-site potential [-ranpot, ranpot] */
  /* = 0.0 */
  double ranpot;
  /* random seed */
  /* = 1 */
  int seed;
} esmac_matrix_Hubbard_params_t;


typedef struct {
  esmac_idx_t ns;
  int maxrowlength;
  esmac_multidx_t *midx;
  esmac_sparserow_real_t * row;
  esmac_dof_fermions_t ** dof;
  esmac_rep_fermions_t ** rep, ** repinit;
  esmac_rng_t *rg;
  double * onsite;
} esmac_matrix_Hubbard_work_t;


esmac_matrix_Hubbard_work_t * esmac_matrix_Hubbard_alloc(const esmac_matrix_Hubbard_params_t * par, int * info);
void esmac_matrix_Hubbard_free(esmac_matrix_Hubbard_work_t *ws);
int esmac_matrix_Hubbard_set_info(const esmac_matrix_Hubbard_params_t * par, esmac_matrix_Hubbard_work_t *ws, esmac_matrix_info_t *info);

int esmac_matrix_Hubbard_row(const esmac_matrix_Hubbard_params_t * par, esmac_matrix_Hubbard_work_t *ws,
				 esmac_idx_t irow, esmac_idx_t *cind, double *val);

#endif /* ESMAC_MATRIX_HUBBARD_H */
