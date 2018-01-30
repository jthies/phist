/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
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

#include <petsc.h>

#endif

//! complex data types
#ifdef __cplusplus
#ifndef DOXYGEN
#include <complex>
#endif
typedef std::complex<float> phist_s_complex;
typedef std::complex<double> phist_d_complex;
#else
#ifndef DOXYGEN
#include <complex.h>
#endif
typedef  complex float phist_s_complex;
typedef  complex double phist_d_complex;
#endif


//! type of node-local indices. 
typedef PetscInt phist_lidx;

//! type of global indices
typedef PetscInt phist_gidx;

#define PRlidx "d"
#define PRgidx "d"

#include "phist_void_aliases.h"

#endif
