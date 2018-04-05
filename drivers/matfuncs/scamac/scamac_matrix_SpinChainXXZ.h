/** \file
 *  \author Andreas Alvermann (<alvermann@physik.uni-greifswald.de>)
 *  \date   October 2017
 *  \brief  ...
 *  \ingroup matrix
 */

#ifndef SCAMAC_MATRIX_SPINCHAINXXZ_H
#define SCAMAC_MATRIX_SPINCHAINXXZ_H

#include "scamac_include.h"
#include "scamac_multidx.h"
#include "scamac_sparserow.h"
#include "scamac_dof_spins.h"

/* >>>>
 * symmetric real
 * one-dimensional XXZ model
 *
 * <<<< */


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
  // SCAMAC_OBC or SCAMAC_PBC
  int boundary_conditions;
} scamac_matrix_SpinChainXXZ_params_st;

typedef struct {
  ScamacIdx ns;
  int maxnzrow;
  scamac_dof_spins_st * dof;
} scamac_matrix_SpinChainXXZ_tables_st;

typedef struct {
  scamac_rep_spins_st *rep, *repinit;
} scamac_matrix_SpinChainXXZ_work_st;


ScamacErrorCode scamac_matrix_SpinChainXXZ_check(const scamac_matrix_SpinChainXXZ_params_st * par, char ** desc);
ScamacErrorCode scamac_matrix_SpinChainXXZ_tables_create(const scamac_matrix_SpinChainXXZ_params_st * par, scamac_matrix_SpinChainXXZ_tables_st ** tab, ScamacInfo * info);
ScamacErrorCode scamac_matrix_SpinChainXXZ_tables_destroy(scamac_matrix_SpinChainXXZ_tables_st * tab);
ScamacErrorCode scamac_matrix_SpinChainXXZ_work_alloc(const scamac_matrix_SpinChainXXZ_params_st * par, const scamac_matrix_SpinChainXXZ_tables_st * tab, scamac_matrix_SpinChainXXZ_work_st ** ws);
ScamacErrorCode scamac_matrix_SpinChainXXZ_work_free(scamac_matrix_SpinChainXXZ_work_st * ws);
ScamacErrorCode scamac_matrix_SpinChainXXZ_generate_row(const scamac_matrix_SpinChainXXZ_params_st * par, const scamac_matrix_SpinChainXXZ_tables_st * tab, scamac_matrix_SpinChainXXZ_work_st * ws, ScamacIdx irow, ScamacFlag flag, scamac_sparserow_real_st * row);

#endif /* SCAMAC_MATRIX_SPINCHAINXXZ_H */
