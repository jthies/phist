/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/

#include "phist_config.h"
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
#include "Tpetra_DefaultPlatform.hpp"

#ifdef PHIST_HAVE_BELOS
#include "Belos_config.h"
# ifdef HAVE_BELOS_TSQR
# include "Tpetra_TsqrAdaptor.hpp"
# include "BelosTpetraAdapter.hpp"
# include "BelosTsqrOrthoManager.hpp"
# endif
#endif

#ifdef PHIST_HAVE_LIKWID
#include <likwid.h>
#endif

#include <fstream>
#include <sstream>


using namespace phist::tpetra;

#include "TpetraCore_config.h"

#if !defined(HAVE_TPETRA_INST_DOUBLE)
#warning "Your TpetraCore_config.h does not define HAVE_TPETRA_INST_DOUBLE, this may cause problems when compiling \
phist. Presently we will try to compile the 'D' kernels anyway, if your Trilinos installation is correct and you run \
into problems compiling phist, please contact the developers because this case is untested."
#endif

#include "phist_gen_d.h"
#include "kernels_def.hpp"
#include "carp_def.hpp"
