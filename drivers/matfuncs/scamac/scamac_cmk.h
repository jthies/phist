#ifndef SCAMAC_CMK_H
#define SCAMAC_CMK_H

#include "scamac_sparsemat.h"

// return Cuthill-McKee permutation
int * scamac_cmk(const scamac_sparsemat_st *sm);
// return permuted matrix
scamac_sparsemat_st * scamac_sparsemat_permute(const scamac_sparsemat_st *sm, const int *perm);

#endif /* SCAMAC_CMK_H */
