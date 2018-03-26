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
/* we support up to 2^6 = 63 error codes + 1 "EOK"
 * The error codes ENULL ... ECORRUPTED can be OR'ed with the position of the problematic input value,
 * that is if the i-th parameter is problematic, with 1 <= i <= 31, functions return, e.g,
 * SCAMAC_EINVALID | (i << SCAMAC_ESHIFT) or SCAMAC_EINVALID | (i << SCAMAC_ESHIFT) | SCAMAC_EINTERNAL
 * The shift is given by SCAMAC_ESHIFT = 7
 */
/* This definition must agree with the list of error names
   in scamac_error_desc() in scamac_collection.c */
#define SCAMAC_ESHIFT 6
#define SCAMAC_EMASK ((1U << SCAMAC_ESHIFT) - 1)
typedef enum {
// +++ basic status
// routine completed successfully
  SCAMAC_EOK=0,
// routine failed for an unspecified reason
  SCAMAC_EFAIL,
// +++
// +++ problems with specific function parameters (on input)
// a null pointer was supplied in the wrong place
  SCAMAC_ENULL,
// an input parameter, or the combination of all parameters, is invalid
  SCAMAC_EINVALID,
// an (integer or double) input parameter is outside of the valid range (e.g, is negative when a positive value is expected)
  SCAMAC_ERANGE,
// an object passed to a routine, probably constructed elsewhere, is corrupted. Can signal an internal error, 
// or improper use of ScaMaC routines & algorithms.
// example: If a complex matrix is passed to a routine that expects a real matrix, EINVALID is returned.
// But if the valtype of the matrix is neither real nor complex, ECORRUPTED is returned, because that's not a valid matrix at all.
  SCAMAC_ECORRUPTED,
// +++
// +++
// the requested computation is outside of the scope of the function (e.g., Lanczos for non-symmetric matrices)
  SCAMAC_ESCOPE,
// the combination of all input parameters is invalid (but each individual parameter may be valid)
  SCAMAC_EINPUT,
// +++
// +++ problems during execution or on return
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
// algorithm (e.g., Lanczos) not converged
  SCAMAC_ENOTCONVERGED,
// row to short
  SCAMAC_ESHORTROW,
// the error originates from an internal call to a ScaMaC routine. This is a flag that can be combined with the other error codes,
// as in SCAMAC_EFAIL | SCAMAC_EINTERNAL
  SCAMAC_EINTERNAL = (1U << (2*SCAMAC_ESHIFT))
} ScamacErrorCode;

// discard warnings (SCAMAC_EWARNING)
#define SCAMAC_DISWARN(err) \
( ((err & SCAMAC_EMASK) == SCAMAC_EWARNING) ? SCAMAC_EOK : err)

// macros for error handling. Replace as you wish.
// Calls exit(EXIT_FAILURE) on abort.
#define SCAMAC_CHKERR(err)						\
  do { \
    if (err) {								\
      fprintf(stderr, "\n*** ABORT ***\n%s\n\nfunction >%s< at line %d failed with\n\n", \
	      __FILE__, __func__,  __LINE__);			\
      fprintf(stderr,"%s%s\n",scamac_error_desc(err),scamac_error_dpar(err)); \
      fprintf(stderr,"*************\n"); \
      exit(EXIT_FAILURE);						\
    }									\
  } while (0)

#define SCAMAC_TRY(exp)							\
  do {									\
    ScamacErrorCode err;						\
    err = exp;								\
    if (err) {								\
      fprintf(stderr, "\n*** ABORT ***\n%s\n\nin function >%s< at line %d\n\n>%s< failed with\n\n", \
	      __FILE__, __func__,  __LINE__, #exp);			\
      fprintf(stderr,"%s%s\n",scamac_error_desc(err),scamac_error_dpar(err)); \
      fprintf(stderr,"*************\n"); \
      exit(EXIT_FAILURE);						\
    }									\
  } while (0)

// macros for error handling with MPI. Replace as you wish.
// Calls MPI_Abort on abort.
#define SCAMAC_CHKERR_MPI(err)						\
  do { \
    if (err) {								\
      int mpi_rank; \
      MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank); \
      fprintf(stderr, "\n*** ABORT *** Process %d\n%s\n\nfunction >%s< at line %d failed with\n\n", \
	      mpi_rank, __FILE__, __func__,  __LINE__);			\
      fprintf(stderr,"%s%s\n",scamac_error_desc(err),scamac_error_dpar(err)); \
      fprintf(stderr,"*************\n"); \
      MPI_Abort(MPI_COMM_WORLD, 1);						\
    }									\
  } while (0)

#define SCAMAC_TRY_MPI(exp)							\
  do {									\
    ScamacErrorCode err;						\
    err = exp;								\
    if (err) {								\
      int mpi_rank; \
      MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank); \
      fprintf(stderr, "\n*** ABORT *** Process %d\n%s\n\nin function >%s< at line %d\n\n>%s< failed with\n\n", \
	      mpi_rank, __FILE__, __func__,  __LINE__, #exp);			\
      fprintf(stderr,"%s%s\n",scamac_error_desc(err),scamac_error_dpar(err)); \
      fprintf(stderr,"*************\n"); \
      MPI_Abort(MPI_COMM_WORLD, 1);						\
    }									\
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
#define SCAMAC_PAR_RNGSEED 4
#define SCAMAC_PAR_OPTION 5

#endif /* SCAMAC_INCLUDE_H */
