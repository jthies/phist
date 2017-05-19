#ifndef ESMAC_MATRIX_FREEBOSONCHAIN_H
#define ESMAC_MATRIX_FREEBOSONCHAIN_H

#include "esmac_types.h"
#include "esmac_dof_bosons.h"
#include "esmac_mstate.h"

typedef struct {
  /* hopping strength */
  /* = 1.0 */
  double t;
  /* number of bosonic species */
  /* = 1 */
  int n_species;
  /* number of sites */
  /* = 10 */
  int n_sites;
  /* number of bosons per species */
  /* = 5 */
  int n_bosons;
  // ESMAC_OBC or ESMAC_PBC
  int boundary_conditions; 
} esmac_matrix_FreeBosonChain_params_t;


typedef struct {
  esmac_idx_t ns;
  int maxrowlength;
  int ndof;
  esmac_dof_bosons_t ** dof;
  esmac_mstate_t *ms;
} esmac_matrix_FreeBosonChain_work_t;


esmac_matrix_FreeBosonChain_work_t * esmac_matrix_FreeBosonChain_alloc(const esmac_matrix_FreeBosonChain_params_t * par, int * info);
void esmac_matrix_FreeBosonChain_free(esmac_matrix_FreeBosonChain_work_t *ws);
int esmac_matrix_FreeBosonChain_set_info(const esmac_matrix_FreeBosonChain_params_t * par, esmac_matrix_FreeBosonChain_work_t *ws, esmac_matrix_info_t *info);

int esmac_matrix_FreeBosonChain_row(const esmac_matrix_FreeBosonChain_params_t * par, esmac_matrix_FreeBosonChain_work_t *ws,
				 esmac_idx_t irow, esmac_idx_t *cind, double *val);

#endif /* ESMAC_MATRIX_FREEBOSONCHAIN_H */
