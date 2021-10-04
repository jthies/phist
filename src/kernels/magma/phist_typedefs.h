/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef PHIST_TYPEDEFS_H
#define PHIST_TYPEDEFS_H

#include "phist_config.h"

#ifndef DOXYGEN

/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <magma.h>
#endif /* DOXYGEN */

#include "phist_macros.h"

//! complex data types
#ifdef __cplusplus
#ifndef DOXYGEN
#include <magma_operators.h>
#include <complex>
#endif /* DOXYGEN */
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

#ifndef DOXYGEN
#include "phist_void_aliases.h"
#endif //DOXYGEN

#endif
