/** \file
 *  \author Andreas Alvermann (<alvermann@physik.uni-greifswald.de>)
 *  \date   October 2017
 *  \brief  ...
 *  \ingroup matrix
 */

#ifndef SCAMAC_MATRIX_HUBBARD_H
#define SCAMAC_MATRIX_HUBBARD_H

#include "scamac_include.h"
#include "scamac_multidx.h"
#include "scamac_sparserow.h"
#include "scamac_dof_fermions.h"
#include "scamac_rng.h"

/* >>>>
 * symmetric real
 * Hubbard models are fermionic
 * solid state models
 *
 * <<<< */

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
  // SCAMAC_OBC or SCAMAC_PBC
  int boundary_conditions;
  /* random on-site potential [-ranpot, ranpot] */
  /* = 0.0 */
  double ranpot;
  /* random seed */
  /* = 1 */
  int seed;
} scamac_matrix_Hubbard_params_st;

typedef struct {
  ScamacIdx ns;
  int maxrowlength;
  scamac_multidx_st *midx;
  scamac_dof_fermions_st ** dof;
  scamac_rng_st *rg;
  double * onsite;
} scamac_matrix_Hubbard_tables_st;

typedef struct {
  scamac_rep_fermions_st ** rep, ** repinit;
} scamac_matrix_Hubbard_work_st;


ScamacErrorCode scamac_matrix_Hubbard_check(const scamac_matrix_Hubbard_params_st * par, char ** desc);
ScamacErrorCode scamac_matrix_Hubbard_tables_create(const scamac_matrix_Hubbard_params_st * par, scamac_matrix_Hubbard_tables_st ** tab, ScamacInfo * info);
ScamacErrorCode scamac_matrix_Hubbard_tables_destroy(scamac_matrix_Hubbard_tables_st * tab);
ScamacErrorCode scamac_matrix_Hubbard_work_alloc(const scamac_matrix_Hubbard_params_st * par, const scamac_matrix_Hubbard_tables_st * tab, scamac_matrix_Hubbard_work_st ** ws);
ScamacErrorCode scamac_matrix_Hubbard_work_free(scamac_matrix_Hubbard_work_st * ws);
ScamacErrorCode scamac_matrix_Hubbard_generate_row(const scamac_matrix_Hubbard_params_st * par, const scamac_matrix_Hubbard_tables_st * tab, scamac_matrix_Hubbard_work_st * ws, ScamacIdx irow, ScamacFlag flag, scamac_sparserow_real_st * row);


#endif /* SCAMAC_MATRIX_HUBBARD_H */
