/** \file
 *  \author Andreas Alvermann (<alvermann@physik.uni-greifswald.de>)
 *  \date   October 2017
 *  \brief  Compute (approximate) bounds on the spectrum
 *  \ingroup toolbox
 */

#ifndef SCAMAC_LANCZOS_H
#define SCAMAC_LANCZOS_H

#include "scamac_sparsemat.h"
#include "scamac_generator.h"

/** \brief Lanczos for extremal eigenvalues of sparse matrix
 *  \ingroup toolbox
 */
ScamacErrorCode scamac_lanczos_ev_mat(const scamac_sparsemat_st *sm, double tol, double *ev1, double *ev2, double *eps1, double *eps2);
/** \brief Lanczos for extremal eigenvalues of sparse matrix, obtained directly from a generator
 *  \ingroup toolbox
 */
ScamacErrorCode scamac_lanczos_ev(const ScamacGenerator *gen, double tol, double *ev1, double *ev2, double *eps1, double *eps2);


// the "working horses"
/** \brief Lanczos for extremal eigenvalues of real (symmetric) matrix (computational routine)
 *  \ingroup toolbox
 */
ScamacErrorCode scamac_lanczos_ev_real(const scamac_sparsemat_st *sm, double tol, double *ev1, double *ev2, double *eps1, double *eps2);
/** \brief Lanczos for extremal eigenvalues of complex (hermitian) matrix (computational routine)
 *  \ingroup toolbox
 */
ScamacErrorCode scamac_lanczos_ev_cplx(const scamac_sparsemat_st *sm, double tol, double *ev1, double *ev2, double *eps1, double *eps2);

#endif /* SCAMAC_LANCZOS_H */
