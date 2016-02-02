#ifndef PHIST_DRIVER_UTILS_H
#define PHIST_DRIVER_UTILS_H

#include "phist_config.h"

#ifndef _ST_
#error "this file should be included after a phist_gen_X header"
#endif

#ifdef MAX
  #undef MAX
  #endif
  #define MAX(a,b) ((a)<(b)?(b):(a))
#ifdef MIN
#undef MIN
#endif
#define MIN(a,b) ((b)<(a)?(b):(a))

typedef _ST_ ST;
typedef _MT_ MT;
typedef TYPE(mvec_ptr) mvec_ptr_t;
typedef TYPE(const_mvec_ptr) const_mvec_ptr_t;

typedef TYPE(sdMat_ptr) sdMat_ptr_t;
typedef TYPE(const_sdMat_ptr) const_sdMat_ptr_t;

typedef TYPE(sparseMat_ptr) sparseMat_ptr_t;
typedef TYPE(const_sparseMat_ptr) const_sparseMat_ptr_t;

#ifdef PHIST_OPERATOR_H
typedef TYPE(linearOp_ptr) op_ptr_t;
typedef TYPE(const_linearOp_ptr) const_op_ptr_t;
#endif

#ifdef __cplusplus
extern "C" {
#endif

//! read matrix from some supported file format

//! auto-detects the file type by looking at the file extension
//! and calls the appropriate kernel routine to read the matrix
void SUBR(sparseMat_read)(TYPE(sparseMat_ptr)* A, const_comm_ptr_t comm, 
        char* filename, int* iflag);

//! quick matrix generation/input routine

//! generate a test matrix described by a string. Currently this
//! is only implemented in real double precision (D), in which case
//! we support e.g. "graphene<L>" (an L x L graphene problem) or "anderson<L>"
//! (an L x L x L Anderson model problem with periodic BC). 
//!
//! Another class of test problems is available if the ESSEX-Physics library is
//! available. If problem starts with "BAPP-", the remainder of the string is
//! passed to the test case generation functions in essex/physics/bapps/.
//!
//! Any other string
//! is assumed to be a matrix filename and passed to sparseMat_read, which recognizes
//! files ending on '.mm', '.bin' or '.rua'.
//!
//! Example: phist_Dcreate_matrix(&A, comm, "anderson42",&iflag) will generate the Anderson
//! model problem on a 42^3 grid, whereas the string "anderson42.mm" will be interpreted
//! as a filename. The special string "usage" causes a usage message to be printed, which
//! contains the supported matrix file formats and the problems implemented at the moment.
void SUBR(create_matrix)(TYPE(sparseMat_ptr)* mat, const_comm_ptr_t comm,
        const char* problem, int* iflag);

//! For testing linear solvers, generates an 'exact solution' sol and right-hand side rhs
//! for some matrix creted by create_matrix. For most test cases, this will be some random
//! sol vector and the rhs is computed by rhs=A*sol, but for the cases stemming from PDEs
//! (BENCH3D-A*,B*) we presribe an analytical solution and generate the correct F for it.
//! The mvecs sol and rhs must be created beforehand and may have an arbitrary number of 
//! columns.
void SUBR(create_sol_and_rhs)(const char* problem, TYPE(const_sparseMat_ptr) A,
                        TYPE(mvec_ptr) sol, TYPE(mvec_ptr) rhs, int* iflag);

int phist_sizeof_lidx();
int phist_sizeof_gidx();

#ifdef __cplusplus
} //extern "C"
#endif

#endif
