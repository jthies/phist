#ifndef PHIST_TYPEDEFS_H
#define PHIST_TYPEDEFS_H

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <magma.h>
#include "phist_macros.h"

//! complex data types
#ifdef __cplusplus
#include <magma_operators.h>
#include <complex>
typedef std::complex<float> phist_s_complex;
typedef std::complex<double> phist_d_complex;
#else
#include <complex.h>
typedef  complex float phist_s_complex;
typedef  complex double phist_d_complex;
//typedef magmaFloatComplex phist_s_complex;
//typedef magmaDoubleComplex phist_d_complex;
#endif


//! type of node-local indices. 
typedef magma_int_t phist_lidx;

//! type of global indices
typedef magma_int_t phist_gidx;

#define PRlidx "d"
#define PRgidx "d"

#include "phist_void_aliases.h"

#endif
