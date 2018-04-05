/** \file
 *  \author Andreas Alvermann (<alvermann@physik.uni-greifswald.de>)
 *  \date   October 2017
 *  \brief  Plot the sparsity pattern
 *  \ingroup toolbox
 */

#ifndef SCAMAC_PLOT_H
#define SCAMAC_PLOT_H

#include "scamac_include.h"
#include "scamac_statistics.h" // for scamac_pattern_st

// plot sparsity pattern "one pixel per entry". Use only for small matrices
ScamacErrorCode scamac_plot_pattern_onetoone(const scamac_matrix_pattern_st * pat, const char * filename);

/** \brief write pattern to file
 *  \ingroup toolbox
 */
ScamacErrorCode scamac_plot_pattern(const scamac_matrix_pattern_st * pat, int downscale, const char * filename);


#endif /* SCAMAC_PLOT_H */
