/** \file
 *  \author Andreas Alvermann (<alvermann@physik.uni-greifswald.de>)
 *  \date   October 2017
 *  \brief  Compute entire spectrum of (sparse) matrix.
 *  \warning Use only for small matrices.
 *  \ingroup toolbox
 */

#ifndef SCAMAC_SPECTRUM_H
#define SCAMAC_SPECTRUM_H

#include "scamac_include.h"
#include "scamac_sparsemat.h"

/** \brief Compute spectrum of real symmetric matrix.
 *  \ingroup toolbox
 * \param[in] sm sparse matrix
 * \param[out] spec pointer to array that contains the spectrum (allocated in routine)
 * \return Error code
 */
ScamacErrorCode scamac_spectrum_real_symmetric(const scamac_sparsemat_st * sm, double **spec);
/** \brief Compute spectrum of general symmetric matrix.
 *  \ingroup toolbox
 * \param[in] sm sparse matrix
 * \param[out] spec pointer to array that contains the spectrum (allocated in routine)
 * \return Error code
 */
ScamacErrorCode scamac_spectrum_real_general  (const scamac_sparsemat_st * sm, double **spec);
/** \brief Compute spectrum of hermitian complex matrix.
 *  \ingroup toolbox
 * \param[in] sm sparse matrix
 * \param[out] spec pointer to array that contains the spectrum (allocated in routine)
 * \return Error code
 */
ScamacErrorCode scamac_spectrum_cplx_hermitian(const scamac_sparsemat_st * sm, double **spec);
/** \brief Compute spectrum of general complex matrix.
 *  \ingroup toolbox
 * \param[in] sm sparse matrix
 * \param[out] spec pointer to array that contains the spectrum (allocated in routine)
 * \return Error code
 */
ScamacErrorCode scamac_spectrum_cplx_general  (const scamac_sparsemat_st * sm, double **spec);

#endif /* SCAMAC_SPECTRUM_H */
