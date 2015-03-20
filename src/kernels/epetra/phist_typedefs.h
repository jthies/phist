#ifndef PHIST_TYPEDEFS_H
#define PHIST_TYPEDEFS_H

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#ifdef __cplusplus

#ifndef DOXYGEN
//#include <cinttypes>
#include <complex>
//! single precision complex type
typedef std::complex<float> s_complex_t;
typedef std::complex<double> d_complex_t;
#endif

#else

#ifndef DOXYGEN
//#include <inttypes.h>
#include <complex.h>
typedef float complex s_complex_t;
typedef double complex d_complex_t;
#endif

#endif

#include "Epetra_config.h"

//! type of node-local indices
typedef int lidx_t;

//! type of global indices
#ifdef EPETRA_NO_64BIT_GLOBAL_INDICES
typedef int gidx_t;
#else
typedef long long int gidx_t;
#endif

#define PRlidx "d"
#define PRgidx "d"

#include "phist_void_aliases.h"

#endif
