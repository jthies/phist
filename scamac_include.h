/** \file
 *  \author Andreas Alvermann (<alvermann@physik.uni-greifswald.de>)
 *  \date   October 2017
 *  \brief  ScaMaC data structure and macro definitions
 *  \ingroup library
 */

#ifndef SCAMAC_INCLUDE_H
#define SCAMAC_INCLUDE_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdbool.h>

#include "scamac_config.h"

// length of names in SCAMAC
#define SCAMAC_NAME_LENGTH 256

/*
 * some integer parameters, e.g., n_sites or n_fermions for the Hubbard example,
 * or the number of sites per axis for a 3D Laplace discretization grid,
 * should be "reasonably small" (even if the resulting matrix is "huge").
 * By definition, an integer n is "reasonably small" if abs(n) <= SCAMACHUGEINT.
 *
 * SCAMACHUGEINT acts mainly as a failsafe against "irresponsible" parameter choices
 * (e.g., when generating matrices in a loop).
 *
 * We set, deliberately, SCAMACHUGEINT = 2^20, but this value is not mandatory.
 * SCAMACHUGEINT should be much smaller than SCAMACIDXMAX.
 *
 * Note: The memory consumption for generator tables or workspace does not exceed
 * a reasonable multiple of SCAMACHUGEINT (here, some MByte).
 *
 */
#define SCAMACHUGEINT (1U << 20)

/* index type */
#ifdef SCAMAC_INDEX_TYPE_int64
typedef int64_t ScamacIdx;
#define SCAMACIDXMAX INT64_MAX
#define SCAMACPRIDX PRId64
#elif defined SCAMAC_INDEX_TYPE_int32
typedef int32_t ScamacIdx;
#define SCAMACIDXMAX INT32_MAX
#define SCAMACPRIDX PRId32
#else
typedef int ScamacIdx;
#error "Unknown SCAMAC_INDEX_TYPE"
#endif

/* error codes */
/* This definition must agree with the list of error names
   in scamac_desc_err() in scamac_collection.c */
typedef enum {
// routine completed successfully
  SCAMAC_EOK=0,
// routine failed for an unspecified reason
  SCAMAC_EFAIL,
// an input parameter has an unknown value (e.g., a flag)
  SCAMAC_EUNKNOWN,
  SCAMAC_ESHORTROW,
// an (integer) input parameter falls out of the valid range (e.g, is negative when a positive value is expected)
  SCAMAC_ERANGE,
// a combination of input parameters is invalid
  SCAMAC_EINVALID,
// a null pointer was supplied in the wrong place
  SCAMAC_ENULL,
// index/matrix dimension exceeds possible integer range (2^31 ~ 2E9 for int32, 2^63 ~ 9E18 for int64)
  SCAMAC_EOVERFLOW,
// memory allocation failed
  SCAMAC_EMALLOCFAIL,
// an integer parameter is too large, i.e., larger then SCAMACHUGE
  SCAMAC_EHUGEINT,
// an integer value computed in the routine is too large, i.e., larger then SCAMACHUGE
  SCAMAC_EHUGECOMP,
// a (parameter) warning. Not necessarily an error
  SCAMAC_EWARNING,
// replace
  SCAMAC_EINVAL,
// algorithm (e.g., Lanczos) not converged
  SCAMAC_ENOTCONVERGED
} ScamacErrorCode;
// the error occurred due to internal problems
#define SCAMAC_EINTERNAL (1U << 10)


// macro for error handling. Replace as you wish.
#define SCAMAC_CHKERR(n) do {if (n) {fprintf(stderr,"%s: Failed with error code: %d [%s]\n",__func__,n,scamac_desc_err(n));exit(EXIT_FAILURE);}} while (0)

#define SCAMAC_TRY(exp)  \
    do { \
      ScamacErrorCode err; \
      err = exp; \
      if (err) { \
        fprintf(stderr, "%s : line %d\n in function >%s<,\n >%s< failed with error code: %d [%s]\n",\
                __FILE__, __LINE__, __func__, #exp, err, scamac_desc_err(err)); \
        exit(EXIT_FAILURE); \
      } \
    } while (0)


/* abstract objects */
typedef struct scamac_generator_st ScamacGenerator;
typedef struct scamac_workspace_st ScamacWorkspace;

typedef struct scamac_info_st ScamacInfo;

// flags for scamac_info_t
typedef int ScamacFlag;
#define SCAMAC_NONE 0
// valtype
#define SCAMAC_VAL_REAL 1
#define SCAMAC_VAL_COMPLEX 2
// symmetry
#define SCAMAC_GENERAL 1
#define SCAMAC_SYMMETRIC 2
#define SCAMAC_HERMITIAN 3

// flags for generate row
#define SCAMAC_DEFAULT 0
#define SCAMAC_TRANSPOSE (1U << 0)
#define SCAMAC_CONJUGATE (1U << 1)
#define SCAMAC_CONJUGATETRANSPOSE (SCAMAC_TRANSPOSE | SCAMAC_CONJUGATE)
#define SCAMAC_KEEPZEROS (1U << 2)
//#define SCAMAC_INT32 (1U << 3)
//#define SCAMAC_SINGLE (1U << 4)

// boundary conditions
#define SCAMAC_OBC 1
#define SCAMAC_PBC 2

// parameter types
//typedef int ScamacPar;
typedef int ScamacPar;
#define SCAMAC_PAR_NONE -1
#define SCAMAC_PAR_INT 1
#define SCAMAC_PAR_DOUBLE 2
#define SCAMAC_PAR_BOOL 3

#endif /* SCAMAC_INCLUDE_H */
