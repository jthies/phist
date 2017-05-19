#ifndef ESMAC_MULTIDX_H
#define ESMAC_MULTIDX_H

#include "esmac_types.h"

typedef struct {
  int n;
  /* i-th index runs from 0 to ni[i]-1 */
  esmac_idx_t *ni;
  /* product of ni[0]*ni[1]* ... */
  esmac_idx_t *niprod;
  /* number of differents index entries */
  esmac_idx_t nidx;
} esmac_multidx_t;

esmac_multidx_t * esmac_multidx_alloc(int n);
void esmac_multidx_free(esmac_multidx_t *midx);

int esmac_multidx_set(esmac_multidx_t *midx, int pos, esmac_idx_t ni);

esmac_idx_t esmac_multidx_nidx(esmac_multidx_t *midx);

esmac_idx_t esmac_multidx_decode(esmac_multidx_t *midx, int pos, esmac_idx_t idx);

esmac_idx_t esmac_multidx_upd(esmac_multidx_t *midx, int pos, esmac_idx_t idx_at_pos, esmac_idx_t idx);

#endif /* ESMAC_MULTIDX_H */
