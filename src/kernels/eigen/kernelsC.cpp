/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
/*! \file Eigen/kernelsC.cpp
 * wraps implementation of Eigen kernel routines
 * \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
 * \author "Jonas Thies <j.thies@tudelft.nl>
 *
*/

#include "phist_config.h"

#if defined(PHIST_HAVE_SP)&&defined(PHIST_HAVE_CMPLX)

/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <fstream>
#include <string>
#include <iostream>

#include "phist_macros.h"
#include "../phist_mpi_kernels.h"
#include "phist_kernel_perfmodels.hpp"
#include "../phist_kernels.h"
#include "phist_lapack.h"

#include "../common/default_context.h"

#include "phist_typedefs.h"
#include "typedefs.hpp"
#include "phist_ScalarTraits.hpp"


#include "phist_gen_c.h"
#include "kernels_def.hpp"
#include "../common/default_context_def.hpp"
#include "../common/default_mvec_get_data_def.hpp"
#include "../common/kernels_no_gpu.cpp"
#include "../common/kernels_no_carp.cpp"
#include "../common/kernels_no_inplace_VC.cpp"
#include "../common/kernels_no_VC_add_WD.cpp"
#include "../common/kernels_no_fused.cpp"

#endif
