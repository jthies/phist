/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_config.h"

#if defined(PHIST_HAVE_SP)&&defined(PHIST_HAVE_CMPLX)

#ifdef PHIST_HAVE_ANASAZI

#include "phist_tools.h"
#include "phist_kernels.h"
#include "phist_operator.h"
#include "phist_ScalarTraits.hpp"
#include "phist_anasazi.h"

#include "phist_AnasaziOperatorTraits.hpp"

// Trilinos stuff
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_FancyOStream.hpp"

// we include the Belos adaptors alongside the Anasazi adapters because
// the TraceMinDavidsonSolMgr requires them.

// some hacks to prevent TSQR orthomanager and some internal saddlepoint vector type
// to clash in Trilinos 12.2.1
#include "Belos_config.h"

#include "Belos_PhistAdapter.hpp"
#include "Anasazi_PhistAdapter.hpp"

#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziSolverManager.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziTraceMinDavidsonSolMgr.hpp"
#endif

#include "phist_gen_c.h"
#include "phist_anasazi_def.hpp"
#endif
