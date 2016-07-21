#ifndef PHIST_TYPEDEFS_H
#define PHIST_TYPEDEFS_H

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#ifndef DOXYGEN
#ifdef __cplusplus
//#include <cinttypes>
#include <complex>
#include <cstddef>
#else
//#include <inttypes.h>
#include <complex.h>
#include <stddef.h>
#endif
#endif

#ifdef __cplusplus
//! single precision complex type
typedef std::complex<float> phist_s_complex;
typedef std::complex<double> phist_d_complex;
//! type of global indices
typedef std::ptrdiff_t phist_gidx;
#else
typedef float complex phist_s_complex;
typedef double complex phist_d_complex;
//! type of global indices
#ifdef PHIST_FORCE_INT_GIDX
typedef int phist_gidx;
#define PRgidx "d"
#else
typedef ptrdiff_t phist_gidx;
#define PRgidx "ld"
#endif
#endif

// we want ptrdiff_t (aka long long int on 64 bit systems) as local index,
// but a bug in Trilinos prevents us from using it right now. until then,
// we use int as local index type

//! type of node-local indices
typedef int phist_lidx;
#define PRlidx "d"

#include "phist_void_aliases.h"

#endif
