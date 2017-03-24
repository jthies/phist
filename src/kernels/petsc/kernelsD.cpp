/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
/*! \file petsc/kernels.cpp
 * wraps implementation of petsc kernel routines
 * \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
 * \author "Jonas Thies <Jonas.Thies@DLR.de>
 *
*/

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <petscsys.h>
#include <petsc/private/petscimpl.h>

#include <fstream>
#include <string>

#include "phist_macros.h"
#include "phist_kernel_perfmodels.hpp"
#include "../phist_kernels.h"
#include "phist_lapack.h"

#include "../common/default_context.h"

#include "phist_typedefs.h"
#include "typedefs.hpp"
#include "phist_ScalarTraits.hpp"


#include "phist_gen_d.h"
#if defined(PETSC_USE_REAL_DOUBLE) && !defined(PETSC_USE_COMPLEX)
#include "kernels_def.hpp"
#include "../common/default_context_def.hpp"
#include "../common/kernels_no_gpu.cpp"
#include "../common/kernels_no_carp.cpp"
#include "../common/kernels_no_inplace_VC.cpp"
#include "../common/kernels_no_VC_add_WD.cpp"
#include "../common/kernels_no_fused.cpp"
#else
#include "../common/kernels_no_impl.cpp"
#endif
