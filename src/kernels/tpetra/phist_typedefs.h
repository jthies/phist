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
#include "Kokkos_DefaultNode.hpp"
#include "Kokkos_OpenMPNode.hpp"
#else
//#include <inttypes.h>
#include <complex.h>
#include <stddef.h>
#endif
#endif

#ifdef __cplusplus
//! TODO - do we want that here or do we want to give
//!        the user a choice?
typedef Kokkos::DefaultNode::DefaultNodeType node_t;
//typedef Kokkos::OpenMPNode node_t;
//! single precision complex type
typedef std::complex<float> s_complex_t;
typedef std::complex<double> d_complex_t;
//! type of global indices
typedef std::ptrdiff_t gidx_t;
#else
typedef float complex s_complex_t;
typedef double complex d_complex_t;
//! type of global indices
typedef ptrdiff_t gidx_t;
#endif

// we want ptrdiff_t (aka long long int on 64 bit systems) as local index,
// but a bug in Trilinos prevents us from using it right now. until then,
// we use int as local index type

//! type of node-local indices
typedef int lidx_t;

#define PRlidx "d"
#define PRgidx "ld"

#include "phist_void_aliases.h"

#endif
