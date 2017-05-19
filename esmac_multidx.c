#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "esmac_constants.h"

#include "esmac_multidx.h"

esmac_multidx_t * esmac_multidx_alloc(int n) {
  esmac_multidx_t * midx = malloc(sizeof *midx);
  midx->n = n;
  midx->ni = calloc(n, sizeof *(midx->ni));
  midx->niprod = calloc(n, sizeof *(midx->niprod));
  midx->nidx=0;
  return midx;
}

void esmac_multidx_free(esmac_multidx_t *midx) {
  if (midx) {
    if (midx->ni) {free(midx->ni);}
    if (midx->niprod) {free(midx->niprod);}
    free(midx);
  }
}

int esmac_multidx_set(esmac_multidx_t *midx, int pos, esmac_idx_t ni) {
  if (0 <= pos && pos < midx->n) {
    midx->ni[pos]=ni;
    // recalculate
    int i;
    midx->niprod[0]=1;
    for (i=1;i<midx->n;i++) {
      midx->niprod[i]=midx->niprod[i-1] * midx->ni[i-1];
    }
    midx->nidx=midx->niprod[midx->n-1] * midx->ni[midx->n-1];
    return 0;
  } else {
    return ESMAC_EINVAL;
  }
}

esmac_idx_t esmac_multidx_nidx(esmac_multidx_t *midx) {
  if (midx) {
    return midx->nidx;
  } else {
    return -1;
  }
}

esmac_idx_t esmac_multidx_decode(esmac_multidx_t *midx, int pos, esmac_idx_t idx) {
  if (0<=idx && idx < midx->nidx) {
    idx = idx / midx->niprod[pos];
    idx = idx % midx->ni[pos];
    return idx;
  } else {
    return -1; // signals error (or does not exist)
  }
}

esmac_idx_t esmac_multidx_upd(esmac_multidx_t *midx, int pos, esmac_idx_t idx_at_pos, esmac_idx_t idx) {
  if (0<= pos && pos < midx->n && 0<= idx_at_pos && idx_at_pos < midx->ni[pos]) {
    esmac_idx_t idxnew, idx_at_old;
    idx_at_old = (idx / midx->niprod[pos]) % midx->ni[pos];
    idxnew = idx + (idx_at_pos-idx_at_old) * midx->niprod[pos];
    return idxnew;
  } else {
    return -1;
  }
}

