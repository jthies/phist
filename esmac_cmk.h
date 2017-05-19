#ifndef ESMAC_CMK_H
#define ESMAC_CMK_H

#include "esmac_sparsemat.h"

// return Cuthill-McKee permutation
int * esmac_cmk(const esmac_sparsemat_t *sm);
// return permuted matrix
esmac_sparsemat_t * esmac_sparsemat_permute(const esmac_sparsemat_t *sm, const int *perm);

#endif /* ESMAC_CMK_H */
