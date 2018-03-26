#ifndef SCAMAC_MATRIX_EXCITON_H
#define SCAMAC_MATRIX_EXCITON_H

#include "complex.h"
#include "scamac_include.h"
#include "scamac_sparserow.h"


/* >>>>
 * hermitian complex
 * Exciton on a lattice
 * <<<< */

typedef enum {Exciton_para, Exciton_ortho} scamac_option_Exciton_sym_ty;

typedef struct {
  /* spin orbit */
  /* = 128.0 */
  double so;
  /* exchange */
  /* = 666.0 */
  double ex; 
  /* mass light hole */
  /* = 0.16 */
  double mlh; 
  /* mass heavy hole */
  /* = 3.1 */
  double mhh; 
  /* mass electron */
  /* = 0.99 */
  double me;
  /* dielectric constant */ 
  /* = 6.94 */
  double eps; 
  /* eff. Coulomb length */
  /* = 1.75 */
  double lc; 
  /* momentum kx */
  /* = 0.0 */
  double kx;
  /* momentum ky */
  /* = 0.0 */
  double ky;
  /* momentum kz */
  /* = 0.0 */
  double kz;
  /* lattice constant */
  /* = 0.42696 */
  double a;
  /* cube length */
  /* = 10 */
  int L;
  /* symmetry */
  /* = Exciton_para */
  scamac_option_Exciton_sym_ty symm;
} scamac_matrix_Exciton_params_st;

typedef struct {
  ScamacIdx ns; // number of rows/columns
  ScamacIdx maxnzrow;
  // additional variables
  // ...
  double par_t1, par_t2, par_te;
  double complex eikx, eiky, eikz;
  ScamacIdx dx,dy,dz;
  double complex mat_so[3][3], mat_ex[3][3];
} scamac_matrix_Exciton_tables_st;

ScamacErrorCode scamac_matrix_Exciton_check        (const scamac_matrix_Exciton_params_st * par, char ** desc);
ScamacErrorCode scamac_matrix_Exciton_tables_create(const scamac_matrix_Exciton_params_st * par, scamac_matrix_Exciton_tables_st ** tab, ScamacInfo * info);
ScamacErrorCode scamac_matrix_Exciton_tables_destroy(scamac_matrix_Exciton_tables_st * tab);
ScamacErrorCode scamac_matrix_Exciton_generate_row (const scamac_matrix_Exciton_params_st * par, const scamac_matrix_Exciton_tables_st * tab, void * ws, ScamacIdx irow, ScamacFlag flag, scamac_sparserow_cplx_st * row);

#endif /* SCAMAC_MATRIX_EXCITON_H */
