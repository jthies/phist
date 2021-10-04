/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/

#include "phist_config.h"

#if defined(PHIST_HAVE_SP)&&defined(PHIST_HAVE_CMPLX)

/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "phist_macros.h"
#include "phist_typedefs.h"
#include "phist_tpetra_typedefs.hpp"
#include "../phist_kernels.h"
#include "phist_trilinos_macros.hpp"
#include "phist_ScalarTraits.hpp"
#include "../common/default_context.h"

#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_DefaultComm.hpp"
#ifdef PHIST_HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"
#endif
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Tpetra_MatrixIO.hpp"

#ifdef PHIST_HAVE_BELOS
#include "Belos_config.h"
# ifdef HAVE_BELOS_TSQR
# include "Tpetra_TsqrAdaptor.hpp"
# include "BelosTpetraAdapter.hpp"
# include "BelosTsqrOrthoManager.hpp"
# endif
#endif

#include "phist_kernel_perfmodels.hpp"

#ifdef PHIST_HAVE_LIKWID
#include <likwid.h>
#endif

#include <fstream>
#include <sstream>


using namespace phist::tpetra;

#include "Teuchos_config.h"
#include "TpetraCore_config.h"
#include "phist_gen_c.h"
# if defined(HAVE_TEUCHOS_FLOAT)&&defined(HAVE_TEUCHOS_COMPLEX)&&defined(HAVE_TPETRA_INST_COMPLEX_FLOAT)
# include "kernels_def.hpp"
# include "carp_def.hpp"
# include "../common/default_mvec_get_data_def.hpp"
# else
# warning "Your phist_config.h defines PHIST_HAVE_COMPLEX and PHIST_HAVE_SP, but your Trilinos installation does not \
support the float complex data type, so all 'C' kernels will return with *iflag=-99 (PHIST_NOT_IMPLEMENTED). In order \
to change this check Teuchos_config.h for HAVE_TEUCHOS_FLOAT and HAVE_TEUCHOS_COMPLEX, and TpetraCore_config.h for \
HAVE_TPETRA_INST_COMPLEX_FLOAT, and install Trilinos such that they are defined. To get rid of this warning, set the \
CMake variables PHIST_ENABLE_SP=OFF and PHISt_ENABLE_COMPLEX=OFF."
# include "../common/kernels_no_impl.cpp"
# endif
#endif
