/** \file
 *  \author Andreas Alvermann (<alvermann@physik.uni-greifswald.de>)
 *  \date   October 2017
 *  \brief  ...
 *  \ingroup matrix
 */

#ifndef SCAMAC_MATRIX_FREEFERMIONCHAIN_H
#define SCAMAC_MATRIX_FREEFERMIONCHAIN_H

#include "scamac_include.h"
#include "scamac_sparserow.h"
#include "scamac_multidx.h"
#include "scamac_dof_fermions.h"

/* >>>>
 * symmetric real
 * free fermions on a chain
 * <<<< */

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
  /* open (false) or periodic (true) boundary conditions */
  /* = true */
  bool PBC;
} scamac_matrix_FreeFermionChain_params_st;

typedef struct {
  ScamacIdx ns; // number of rows/columns
  ScamacIdx maxnzrow;
  int ndof;
  scamac_multidx_st *midx;
  scamac_dof_fermions_st ** dof;
} scamac_matrix_FreeFermionChain_tables_st;

typedef struct {
  int ndof;
  scamac_rep_fermions_st ** rep, ** repinit;
} scamac_matrix_FreeFermionChain_work_st;

ScamacErrorCode scamac_matrix_FreeFermionChain_check        (const scamac_matrix_FreeFermionChain_params_st * par, char ** desc);
ScamacErrorCode scamac_matrix_FreeFermionChain_tables_create(const scamac_matrix_FreeFermionChain_params_st * par, scamac_matrix_FreeFermionChain_tables_st ** tab, ScamacInfo * info);
ScamacErrorCode scamac_matrix_FreeFermionChain_tables_destroy(scamac_matrix_FreeFermionChain_tables_st * tab);
ScamacErrorCode scamac_matrix_FreeFermionChain_work_alloc   (const scamac_matrix_FreeFermionChain_params_st * par, const scamac_matrix_FreeFermionChain_tables_st * tab, scamac_matrix_FreeFermionChain_work_st ** ws);
ScamacErrorCode scamac_matrix_FreeFermionChain_work_free    (scamac_matrix_FreeFermionChain_work_st * ws);
ScamacErrorCode scamac_matrix_FreeFermionChain_generate_row (const scamac_matrix_FreeFermionChain_params_st * par, const scamac_matrix_FreeFermionChain_tables_st * tab, scamac_matrix_FreeFermionChain_work_st * ws, ScamacIdx irow, ScamacFlag flag, scamac_sparserow_real_st * row);

#endif /* SCAMAC_MATRIX_FREEFERMIONCHAIN_H */
