/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_config.h"

#include "phist_tools.h"
#include "phist_kernels.h"
#include "phist_operator.h"
#include "phist_ScalarTraits.hpp"
#include "phist_anasazi.h"

#ifdef PHIST_HAVE_ANASAZI

// this solver is to unstable right now, it doesn't compile with Trilinos 12.x if TSQR is enabled
//#define PHIST_HAVE_ANASAZI_TRACEMIN_DAVIDSON

#include "phist_rcp_helpers.hpp"
#include "phist_AnasaziOperatorTraits.hpp"

// Trilinos stuff
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_FancyOStream.hpp"

# ifdef PHIST_KERNEL_LIB_GHOST
#  include "ghost.h"
#  include "Belos_GhostAdapter.hpp"
#  include "Anasazi_GhostAdapter.hpp"
//#  include "Ghost_TsqrAdapter.hpp"
# elif defined(PHIST_KERNEL_LIB_EPETRA)
#  include "Epetra_MultiVector.h"
#  include "BelosEpetraAdapter.hpp"
#  include "AnasaziEpetraAdapter.hpp"
# elif defined(PHIST_KERNEL_LIB_TPETRA)
#  include "Tpetra_MultiVector.hpp"
#  include "BelosTpetraAdapter.hpp"
#  include "AnasaziTpetraAdapter.hpp"
#  include "phist_tpetra_typedefs.hpp"
# else
#  warning "Anasazi not supported for this kernel lib"
#  undef PHIST_HAVE_ANASAZI
# endif
#endif

// Block Krylov-Schur solver manager from the Anasazi package
#ifdef PHIST_HAVE_ANASAZI
// adaptation of a basic ortho class from Trili 11.12 to avoid 
// col-wise norm calculations
//#include "phist_AnasaziMatOrthoManager.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziSolverManager.hpp"
/* use our own adaptation of this file from Trilinos 11.12.1 because 
 * the original did not support TSQR
 */
/*#include "phist_AnasaziBlockKrylovSchurSolMgr.hpp"*/
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
# ifdef PHIST_HAVE_ANASAZI_TRACEMIN_DAVIDSON
# include "AnasaziTraceMinDavidsonSolMgr.hpp"
# endif
#endif

#include "phist_gen_d.h"
#include "phist_anasazi_def.hpp"

