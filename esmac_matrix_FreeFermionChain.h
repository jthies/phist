#ifndef ESMAC_MATRIX_FREEFERMIONCHAIN_H
#define ESMAC_MATRIX_FREEFERMIONCHAIN_H

#include "esmac_types.h"
#include "esmac_multidx.h"
#include "esmac_sparserow.h"
#include "esmac_dof_fermions.h"

typedef struct {
  /* hopping strength */
  /* = 1.0 */
  double t;
  /* number of fermionic species */
  /* = 1 */
  int n_species;
  /* number of sites */
  /* = 10 */
  int n_sites;
  /* number of fermions per species */
  /* = 5 */
  int n_fermions;
  // ESMAC_OBC or ESMAC_PBC
  int boundary_conditions; 
} esmac_matrix_FreeFermionChain_params_t;


typedef struct {
  esmac_idx_t ns;
  int maxrowlength;
  int ndof;
  esmac_multidx_t *midx;
  esmac_sparserow_real_t * row;
  esmac_dof_fermions_t ** dof;
  esmac_rep_fermions_t ** rep, ** repinit;
} esmac_matrix_FreeFermionChain_work_t;


esmac_matrix_FreeFermionChain_work_t * esmac_matrix_FreeFermionChain_alloc(const esmac_matrix_FreeFermionChain_params_t * par, int * info);
void esmac_matrix_FreeFermionChain_free(esmac_matrix_FreeFermionChain_work_t *ws);
int esmac_matrix_FreeFermionChain_set_info(const esmac_matrix_FreeFermionChain_params_t * par, esmac_matrix_FreeFermionChain_work_t *ws, esmac_matrix_info_t *info);

int esmac_matrix_FreeFermionChain_row(const esmac_matrix_FreeFermionChain_params_t * par, esmac_matrix_FreeFermionChain_work_t *ws,
				 esmac_idx_t irow, esmac_idx_t *cind, double *val);

#endif /* ESMAC_MATRIX_FREEFERMIONCHAIN_H */
