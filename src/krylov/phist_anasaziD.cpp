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
#include "phist_MemOwner.hpp"
#include "phist_anasazi.h"

#ifdef PHIST_HAVE_ANASAZI


#include "phist_types.hpp"
#include "phist_kernels.hpp"
#include "phist_core.hpp"
#include "phist_ScalarTraits.hpp"
// glue code
#include "Belos_PhistAdapter.hpp"
#include "Anasazi_PhistAdapter.hpp"
#include "phist_AnasaziOperatorTraits.hpp"

// Trilinos stuff
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_FancyOStream.hpp"

#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziSVQBOrthoManager.hpp"
#include "AnasaziSolverManager.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziTraceMinDavidsonSolMgr.hpp"
#include "AnasaziLOBPCGSolMgr.hpp"
#endif

#include "phist_gen_d.h"
#include "phist_anasazi_def.hpp"

