/** \file
 *  \author Andreas Alvermann (<alvermann@physik.uni-greifswald.de>)
 *  \date   October 2017
 *  \brief  Sparse matrix creation and manipulation
 *  \ingroup toolbox
 */

#ifndef SCAMAC_SPARSEMAT_H
#define SCAMAC_SPARSEMAT_H

#include <complex.h>
#include <stdbool.h>

#include "scamac_include.h"
#include "scamac_generator.h"

typedef struct {
  /// number of rows
  ScamacIdx nr;
  /// number of columns
  ScamacIdx nc;
  /// number of non-zero (i.e. stored) matrix entries
  ScamacIdx ne;
  /// maximal number of non-zeroes
  ScamacIdx nemax;
  /// row pointers (dimension NR+1)
  ScamacIdx *rptr;
  /// column indices (dimension NE)
  ScamacIdx *cind;
  /// values (dimension NE)
  double *val;
  ///
  /// type of entries
  int valtype; // SCAMAC_VAL_REAL or SCAMAC_VAL_COMPLEX
} scamac_sparsemat_st;


/*
typedef struct {
  /// number of rows
  int nr;
  /// number of columns
  int nc;
  /// number of non-zero (i.e. stored) matrix entries
  int ne;
  /// maximal number of non-zeroes
  // int nemax;
  /// row pointers (dimension NR+1)
  int *rptr;
  /// column indices (dimension NE)
  int *cind;
  /// values (dimension NE)
  double *val;
  ///
  /// type of entries
  int valtype; // SCAMAC_VAL_REAL or SCAMAC_VAL_COMPLEX
} scamac_sparsemat_st;
*/

/** \brief Allocate sparse matrix
 *  \ingroup toolbox
 */
ScamacErrorCode scamac_sparsemat_alloc(ScamacIdx nr, ScamacIdx nc, ScamacIdx ne, int valtype, scamac_sparsemat_st ** sm);
/** \brief Free allocated sparse matrix
 *  \ingroup toolbox
 */
ScamacErrorCode scamac_sparsemat_free(scamac_sparsemat_st * sm);
/** \brief Obtain sparse matrix from ScaMaC generator
 *  \ingroup toolbox
 */
ScamacErrorCode scamac_sparsemat_from_generator(const ScamacGenerator * gen, scamac_sparsemat_st ** sm);

/** \brief Sparse matrix-vector multiplication: y = alpha SM x + beta y + gamma x
 *  \ingroup toolbox
 */
ScamacErrorCode scamac_sparsemat_mvm(const scamac_sparsemat_st *sm, const double *x, double *y, double alpha, double beta, double gamma);

/** \brief Complex sparse matrix-vector multiplication: y = alpha SM x + beta y + gamma x
 *  \ingroup toolbox
 */
ScamacErrorCode scamac_sparsemat_mvm_cplx(const scamac_sparsemat_st *sm, const double complex *x, double complex *y,
    double complex alpha, double complex beta, double complex gamma);

ScamacIdx scamac_sparsemat_maxrowlength(const scamac_sparsemat_st *sm);

/** \brief query symmetry (hermiticity) of matrix
 *  \ingroup toolbox
 *  \pre assumes \em ordered, i.e., \em monotonically increasing column indices. This is true for matrices obtained from the ScaMaC generators with scamac_sparsemat_from_generator()
 *  \param[in] sm sparse matrix to be checked
 *  \param[out] symm_pattern Is the sparsity pattern symmetric?
 *  \param[out] symm_value Is the matrix symmetric (for real sm) or Hermitian (for complex sm)?
 *  \note symm_pattern==false implies symm_value==false
 *  \return Error code
 *  \todo Check for zero entries!
 */
ScamacErrorCode scamac_sparsemat_check_symmetry(const scamac_sparsemat_st * sm, bool * symm_pattern, bool * symm_value);

#endif /* SCAMAC_SPARSEMAT_H */
