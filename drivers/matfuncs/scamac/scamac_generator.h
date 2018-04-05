/** \file
 *  \author Andreas Alvermann (<alvermann@physik.uni-greifswald.de>)
 *  \date   October 2017
 *  \brief  generic ScaMaC generator routines
 *  \ingroup library
 */

#ifndef SCAMAC_GENERATOR_H
#define SCAMAC_GENERATOR_H

#include <complex.h>

#include "scamac_include.h"

/** \brief   Check the parameters of the generator
 *  \details This routine performs some checks.
 *  \param[in] gen The generator to be checked
 *  \param[out] desc  If desc != NULL, it contains a string describing the problems with the parameters.
 *  \return Status of the check. Returns SCAMAC_EOK if check passed, SCAMAC_EFAIL otherwise.
 *  \pre gen != NULL
 *  \ingroup library
 */
ScamacErrorCode scamac_generator_check(const ScamacGenerator * gen, char ** desc);
ScamacErrorCode scamac_generator_finalize(ScamacGenerator * gen);
ScamacErrorCode scamac_generator_destroy(ScamacGenerator * gen);

ScamacErrorCode scamac_workspace_alloc(const ScamacGenerator * gen, ScamacWorkspace ** ws);
ScamacErrorCode scamac_workspace_free (ScamacWorkspace * ws);

/** \brief   Generate one row of a matrix.
 *  \details We generate a row.
 *  \param[in] gen Matrix generator, after scamac_generator_finalize()
 *  \param[inout] ws Workspace, allocated by scamac_workspace_alloc()
 *  \param[in] irow index of the row to be generated
 *  \param[in] flag flags
 *  \param[out] nzr number of non-zero elements in the row
 *  \param[out] cind column indices of non-zero elements
 *  \param[out] val values of non-zero elements
 *  \ingroup library
 */
ScamacErrorCode scamac_generate_row(const ScamacGenerator * gen, ScamacWorkspace * ws, ScamacIdx irow, ScamacFlag flag,
                                    ScamacIdx * nzr, ScamacIdx * cind, double * val);

/* more specialized versions */
ScamacErrorCode scamac_generate_row_real    (const ScamacGenerator * gen, ScamacWorkspace * ws, ScamacIdx irow, ScamacFlag flag, ScamacIdx * nzr, ScamacIdx * cind, double * val);
ScamacErrorCode scamac_generate_row_cplx    (const ScamacGenerator * gen, ScamacWorkspace * ws, ScamacIdx irow, ScamacFlag flag, ScamacIdx * nzr, ScamacIdx * cind, double complex * val);
ScamacErrorCode scamac_generate_row_int     (const ScamacGenerator * gen, ScamacWorkspace * ws, ScamacIdx irow, ScamacFlag flag, int * nzr, int * cind, double * val);
ScamacErrorCode scamac_generate_row_int_real(const ScamacGenerator * gen, ScamacWorkspace * ws, ScamacIdx irow, ScamacFlag flag, int * nzr, int * cind, double * val);
ScamacErrorCode scamac_generate_row_int_cplx(const ScamacGenerator * gen, ScamacWorkspace * ws, ScamacIdx irow, ScamacFlag flag, int * nzr, int * cind, double complex * val);

/** \brief Allocate vector cind and val, for calls to scamac_generate_row().
 * \details This function is merely a convenience.
 * \ingroup library
 */
ScamacErrorCode scamac_alloc_cind_val(const ScamacGenerator * gen, ScamacFlag flag, ScamacIdx ** cind, double ** val);

const char * scamac_generator_query_name (const ScamacGenerator * gen);

ScamacIdx scamac_generator_query_nrow    (const ScamacGenerator * gen);
ScamacIdx scamac_generator_query_ncol    (const ScamacGenerator * gen);
ScamacIdx scamac_generator_query_maxnzrow(const ScamacGenerator * gen);
ScamacIdx scamac_generator_query_maxnzcol(const ScamacGenerator * gen);
ScamacIdx scamac_generator_query_maxnz   (const ScamacGenerator * gen);
ScamacIdx scamac_generator_query_valtype (const ScamacGenerator * gen);
ScamacIdx scamac_generator_query_symmetry(const ScamacGenerator * gen);


#endif /* SCAMAC_GENERATOR_H */
