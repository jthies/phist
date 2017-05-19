#ifndef ESMAC_MATRIX_SPINCHAINXXZ_H
#define ESMAC_MATRIX_SPINCHAINXXZ_H

#include "esmac_types.h"
#include "esmac_multidx.h"
#include "esmac_sparserow.h"
#include "esmac_dof_spins.h"

typedef struct {
  /* J_x=J_y */
  /* = 1.0 */
  double Jxy;
  /* J_z */
  /* = 1.0 */
  double Jz;
  /* Bz */
  /* = 0.0 */
  double Bz;
  /* number of sites */
  /* = 10 */
  int n_sites;
  /* number of _up_ spins */
  /* = 5 */
  int n_up;
  // ESMAC_OBC or ESMAC_PBC
  int boundary_conditions; 
} esmac_matrix_SpinChainXXZ_params_t;


typedef struct {
  esmac_idx_t ns;
  int maxrowlength;
  esmac_sparserow_real_t * row;
  esmac_dof_spins_t * dof;
  esmac_rep_spins_t *rep, *repinit;
} esmac_matrix_SpinChainXXZ_work_t;


esmac_matrix_SpinChainXXZ_work_t * esmac_matrix_SpinChainXXZ_alloc(const esmac_matrix_SpinChainXXZ_params_t * par, int * info);
void esmac_matrix_SpinChainXXZ_free(esmac_matrix_SpinChainXXZ_work_t *ws);
int esmac_matrix_SpinChainXXZ_set_info(const esmac_matrix_SpinChainXXZ_params_t * par, esmac_matrix_SpinChainXXZ_work_t *ws, esmac_matrix_info_t *info);

int esmac_matrix_SpinChainXXZ_row(const esmac_matrix_SpinChainXXZ_params_t * par, esmac_matrix_SpinChainXXZ_work_t *ws,
				 esmac_idx_t irow, esmac_idx_t *cind, double *val);

#endif /* ESMAC_MATRIX_SPINCHAINXXZ_H */
