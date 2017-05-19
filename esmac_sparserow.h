#ifndef ESMAC_SPARSEROW_H
#define ESMAC_SPARSEROW_H

#include <complex.h>

#include "esmac_constants.h"
#include "esmac_types.h"

typedef struct {
  esmac_idx_t idx;
  double val;
} esmac_idxreal_t;

typedef struct {
  int nalloc;
  int nrow;
  esmac_idxreal_t *row; 
} esmac_sparserow_real_t;

typedef struct {
  esmac_idx_t idx;
  double complex val;
} esmac_idxcplx_t;

typedef struct {
  int nalloc;
  int nrow;
  esmac_idxcplx_t *row;
} esmac_sparserow_cplx_t;


/* sparse row routines */

esmac_sparserow_real_t * esmac_sparserow_real_alloc();
esmac_sparserow_cplx_t * esmac_sparserow_cplx_alloc();

void esmac_sparserow_real_free(esmac_sparserow_real_t *row);
void esmac_sparserow_cplx_free(esmac_sparserow_cplx_t *row);

void esmac_sparserow_real_zero(esmac_sparserow_real_t *row);
void esmac_sparserow_cplx_zero(esmac_sparserow_cplx_t *row);

int esmac_sparserow_real_add(esmac_sparserow_real_t *row, double alpha, esmac_idx_t idx);
int esmac_sparserow_cplx_add(esmac_sparserow_cplx_t *row, double complex alpha, esmac_idx_t idx);

/* convert sparse row to idx/val pair of vectors. Sort and compress at the same time. Throw FATAL error if maxrowlength is exceeded */
int esmac_sparserow_real_to_idxval(esmac_sparserow_real_t *row, int maxrowlength, int keep_zeros, esmac_idx_t *idx, double *val);
int esmac_sparserow_cplx_to_idxval(esmac_sparserow_cplx_t *row, int maxrowlength, int keep_zeros, esmac_idx_t *idx, double complex *val);


#endif /* ESMAC_SPARSEROW_H */
