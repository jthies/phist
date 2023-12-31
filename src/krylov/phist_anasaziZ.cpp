/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_config.h"

#if defined(PHIST_HAVE_CMPLX)

/* check wether we have to disable Anasazi because the Trilinos installation doesn't support this data type */
#if defined(PHIST_KERNEL_LIB_TPETRA)
#include "TpetraCore_config.h"
# ifndef HAVE_TPETRA_INST_COMPLEX_DOUBLE
# undef PHIST_HAVE_ANASAZI
# endif
#elif defined(PHIST_HAVE_TEUCHOS)
# include "Teuchos_config.h"
# ifndef HAVE_TEUCHOS_COMPLEX
# undef PHIST_HAVE_ANASAZI
# endif
#endif

#include "phist_tools.h"
#include "phist_MemOwner.hpp"
#include "phist_anasazi.h"

#ifdef PHIST_HAVE_ANASAZI

#include "phist_types.hpp"
#include "phist_kernels.hpp"
#include "phist_core.hpp"
#include "phist_ScalarTraits.hpp"

// we include the Belos adaptors alongside the Anasazi adapters because
// the TraceMinDavidsonSolMgr requires them.
#include "Belos_PhistAdapter.hpp"
#include "Anasazi_PhistAdapter.hpp"
#include "phist_AnasaziOperatorTraits.hpp"

// Trilinos stuff
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_FancyOStream.hpp"

#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziSolverManager.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
# include "AnasaziTraceMinDavidsonSolMgr.hpp"
# include "AnasaziLOBPCGSolMgr.hpp"
#endif

#include "phist_gen_z.h"
#include "phist_anasazi_def.hpp"

#endif
