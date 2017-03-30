/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef MATFUNCS_H
#define MATFUNCS_H

#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <stdint.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "phist_kernels.h"

#ifdef PHIST_HAVE_GHOST
#include <ghost.h>
#else
#define PRGIDX PRId64
#define PRLIDX PRId32
// just dummy values, simply ignored!
#define GHOST_SPARSEMAT_SYMM_GENERAL 0x1
#define GHOST_DT_DOUBLE 0x2
#define GHOST_DT_REAL 0x4
#define GHOST_IDX64_GLOBAL
#endif


typedef struct {

	int32_t version;
	int32_t base;
	int32_t symmetry;
	int32_t datatype;
	ghost_gidx nrows;
	ghost_gidx ncols;
	ghost_lidx row_nnz;

	int32_t      hermit;
	int32_t      eig_info;
	double       eig_down;
	double       eig_top;

	} matfuncs_info_t;


#ifdef __cplusplus
extern "C" {
#endif

int SpinChainSZ(   ghost_gidx row, ghost_lidx *nnz, ghost_gidx *cols, void *vals, void* data);
int crsGraphene(   ghost_gidx row, ghost_lidx *nnz, ghost_gidx *cols, void *vals, void* data);
int MATPDE_rowFunc(   ghost_gidx row, ghost_lidx *nnz, ghost_gidx *cols, void *vals, void* data);
int MATPDE3D_rowFunc(   ghost_gidx row, ghost_lidx *nnz, ghost_gidx *cols, void *vals, void* data);

#ifdef __cplusplus
}
#endif

#endif
