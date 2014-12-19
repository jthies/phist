#ifndef PHIST_TYPEDEFS_H
#define PHIST_TYPEDEFS_H

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include <inttypes.h>

//! complex data types
#ifdef __cplusplus
#include <complex>
typedef std::complex<float> s_complex_t;
typedef std::complex<double> d_complex_t;
#else
#include <complex.h>
typedef  complex float s_complex_t;
typedef  complex double d_complex_t;
#endif

//! type of node-local indices
typedef int32_t lidx_t;

//! type of global indices
typedef int64_t gidx_t;

#define PRlidx "d"
#define PRgidx "ld"

#include "phist_void_aliases.h"

#endif