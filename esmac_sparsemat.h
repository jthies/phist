#ifndef ESMAC_SPARSEMAT_H
#define ESMAC_SPARSEMAT_H

#include "esmac_types.h"
#include "esmac_generator.h"

typedef struct {
  /// number of rows
  esmac_idx_t nr;
  /// number of columns
  esmac_idx_t nc;
  /// number of non-zero (i.e. stored) matrix entries
  esmac_idx_t ne;
  /// maximal number of non-zeroes
  esmac_idx_t nemax;
  /// row pointers (dimension NR+1)
  esmac_idx_t *rptr;
  /// column indices (dimension NE)
  esmac_idx_t *cind;
  /// values (dimension NE)
  double *val;
  ///
  /// type of entries
  int valtype; // ESMAC_VAL_REAL or ESMAC_VAL_COMPLEX
} esmac_sparsemat_t;


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
  int valtype; // ESMAC_VAL_REAL or ESMAC_VAL_COMPLEX
} esmac_sparsemat_t;
*/

esmac_sparsemat_t * esmac_sparsemat_alloc(int nr, int nc, int ne);
void esmac_sparsemat_free(esmac_sparsemat_t *sm);
esmac_sparsemat_t * esmac_sparsemat_from_generator(const esmac_generator_t *gen);

// y = alpha SM x + beta y + gamma x 
int esmac_sparsemat_mvm(const esmac_sparsemat_t *sm, const double *x, double *y, double alpha, double beta, double gamma);

int esmac_sparsemat_maxrowlength(const esmac_sparsemat_t *sm);

#endif /* ESMAC_SPARSEMAT_H */
