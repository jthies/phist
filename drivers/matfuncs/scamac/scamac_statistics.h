/** \file
 *  \author Andreas Alvermann (<alvermann@physik.uni-greifswald.de>)
 *  \date   October 2017
 *  \brief  Pattern and value statistics
 *  \ingroup toolbox
 */

#ifndef SCAMAC_STATISTICS_H
#define SCAMAC_STATISTICS_H

#include <stdio.h>
#include <stdbool.h>

#include "scamac_generator.h"

typedef struct {
  // matrix data
  ScamacIdx nrow,ncol;
  // number of rows incorporated into the statistics so far
  ScamacIdx ncontributed;
  // valtype
  int valtype;

  // pattern statistics
  ScamacIdx n_nz;
  ScamacIdx n_nz_left, n_nz_right;
  ScamacIdx n_nz_row_min, n_nz_row_max;
  ScamacIdx n_nz_row_max_left, n_nz_row_max_right; // only max - min always = 0
  ScamacIdx n_zero_row, n_zero_diag;
  // bandwidth
  ScamacIdx bw_left,bw_right;

  // value statistics
  double v_min_re,v_max_re;
  double v_min_im,v_max_im;
  double v_abs;
  double v_min_re_diag,v_max_re_diag;
  double v_min_im_diag,v_max_im_diag;
  double v_abs_diag;
  // number of diagonally dominant rows
  ScamacIdx n_diag_dominant;
  // minimal |diag|-sum |offdiag|. Positive if and only if all rows are diagonally dominant
  double diag_minus_offdiag;

  // spectral statistics
  double gershgorin_min_re, gershgorin_max_re, gershgorin_min_im, gershgorin_max_im;

} scamac_matrix_statistics_st;

typedef struct {
  // matrix data
  ScamacIdx nrow, ncol;
  // pattern size
  int px,py;
  // number of rows incorporated into the pattern so far
  ScamacIdx ncontributed;
  // the pattern (array dimension px * py)
  ScamacIdx * pat;
} scamac_matrix_pattern_st;

/** \brief  collect matrix statistics and/or pattern from generator
 *  \ingroup toolbox
 *  \details A matrix statistics is a small object, and no allocation of the struct is necessary.
 *           A matrix pattern is a large object, and will be allocated in this routine.
 *           It should be free'd with scamac_pattern_free() afterwards.
 *  \todo Enable OpenMP support.
 */
ScamacErrorCode scamac_collect_statistics_and_pattern(const ScamacGenerator * gen, ScamacFlag flag, scamac_matrix_statistics_st * st, scamac_matrix_pattern_st ** pt);

/** \brief  set statistics to "empty", i.e., initialize
 *  \ingroup toolbox
 */
ScamacErrorCode scamac_statistics_empty  (scamac_matrix_statistics_st * st, ScamacIdx nrow, ScamacIdx ncol, int valtype);
/** \brief  update statistics with data from a matrix row not yet included
 *  \ingroup toolbox
 */
ScamacErrorCode scamac_statistics_update (scamac_matrix_statistics_st * st, ScamacIdx irow, ScamacIdx nzr, const ScamacIdx * cind, const double * val);
/** \brief  combine two parts of a matrix statistics for different sets of rows
 *  \ingroup toolbox
 *  \todo allow for val == NULL
 */
ScamacErrorCode scamac_statistics_combine(scamac_matrix_statistics_st * stcomb, const scamac_matrix_statistics_st * st);
/** \brief  create description of statistics for output, printing etc.
 *  \ingroup toolbox
 */
ScamacErrorCode scamac_statistics_print  (const scamac_matrix_statistics_st * st, char ** desc);


/** \brief  allocate memory for pattern
 *  \ingroup toolbox
 */
ScamacErrorCode scamac_pattern_alloc  (int px, int py, scamac_matrix_pattern_st ** pt);
/** \brief  set pattern to "empty", i.e., initialize
 *  \ingroup toolbox
 */
ScamacErrorCode scamac_pattern_empty  (scamac_matrix_pattern_st * pt, ScamacIdx nrow, ScamacIdx ncol, int valtype);
/** \brief  update pattern with data from a matrix row not yet included
 *  \ingroup toolbox
 */
ScamacErrorCode scamac_pattern_update (scamac_matrix_pattern_st * pt, ScamacIdx irow, ScamacIdx nzr, const ScamacIdx * cind);
/** \brief  combine two parts of a sparsity pattern for different sets of rows
 *  \ingroup toolbox
 */
ScamacErrorCode scamac_pattern_combine(scamac_matrix_pattern_st * ptcomb, const scamac_matrix_pattern_st * pt);
/** \brief  create ASCII pattern for output, printing etc.
 *  \note   The real plot functions are found in scamac_plot.h
 *  \ingroup toolbox
 */
ScamacErrorCode scamac_pattern_print  (const scamac_matrix_pattern_st * pt, char ** desc);
/** \brief  free memory allocated for pattern
 *  \ingroup toolbox
 */
ScamacErrorCode scamac_pattern_free   (scamac_matrix_pattern_st * pt);

#endif /* SCAMAC_STATISTICS_H */
