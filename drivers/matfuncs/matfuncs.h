#ifndef MATFUNCS_H
#define MATFUNCS_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef PHIST_HAVE_GHOST
#include <ghost.h>
#else
typedef int64_t ghost_idx_t;
#define PRIDX "PRId64"
// just dummy values, simply ignored!
#define GHOST_SPARSEMAT_SYMM_GENERAL 0x1
#define GHOST_DT_DOUBLE 0x2
#define GHOST_DT_REAL 0x4
#define GHOST_HAVE_LONGIDX
#endif


typedef struct {

	int32_t version;
	int32_t base;
	int32_t symmetry;
	int32_t datatype;
	ghost_idx_t nrows;
	ghost_idx_t ncols;
	ghost_idx_t row_nnz;

	int32_t      hermit;
	int32_t      eig_info;
	double       eig_down;
	double       eig_top;

	} matfuncs_info_t;


#ifdef __cplusplus
extern "C" {
#endif

int SpinChainSZ(   ghost_idx_t row, ghost_idx_t *nnz, ghost_idx_t *cols, void *vals);
int crsGraphene(   ghost_idx_t row, ghost_idx_t *nnz, ghost_idx_t *cols, void *vals);

#ifdef __cplusplus
}
#endif

#endif
