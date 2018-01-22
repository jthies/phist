/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_config.h"

#ifdef PHIST_KERNEL_LIB_EPETRA
#undef PHIST_HAVE_ANASAZI
#elif defined(PHIST_KERNEL_LIB_TPETRA)
#include "TpetraCore_config.h"
# ifndef HAVE_TPETRA_INST_FLOAT
# undef PHIST_HAVE_ANASAZI
# endif
#endif
#ifdef PHIST_HAVE_SP

#include "phist_tools.h"
#include "phist_kernels.h"
#include "phist_operator.h"
#include "phist_ScalarTraits.hpp"
#include "phist_anasazi.h"

#ifdef PHIST_HAVE_ANASAZI

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
#endif

#include "phist_gen_s.h"
#include "phist_anasazi_def.hpp"
#endif
