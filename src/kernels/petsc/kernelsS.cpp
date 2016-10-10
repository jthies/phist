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

#include "phist_typedefs.h"
#include "typedefs.hpp"
#include "phist_ScalarTraits.hpp"


#include "phist_gen_s.h"
#if defined(PETSC_USE_REAL_SINGLE) && !defined(PETSC_USE_COMPLEX)
#include "kernels_def.hpp"
#include "../common/kernels_no_gpu.cpp"
#include "../common/kernels_no_carp.cpp"
#include "../common/kernels_no_inplace_VC.cpp"
#include "../common/kernels_no_VC_add_WD.cpp"
#include "../common/kernels_no_fused.cpp"
#else
#include "../common/kernels_no_impl.cpp"
#endif
