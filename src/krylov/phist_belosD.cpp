/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_config.h"

#include "phist_tools.h"
#include "phist_belos.h"
#include "phist_MemOwner.hpp"

#ifdef PHIST_HAVE_BELOS

#include "phist_kernels.hpp"
#include "phist_core.hpp"
#include "phist_ScalarTraits.hpp"

// Trilinos stuff
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_FancyOStream.hpp"

#include "phist_BelosMV.hpp"
#include "Belos_PhistAdapter.hpp"
#include "phist_BelosOperatorTraits.hpp"

#include "BelosICGSOrthoManager.hpp"
#include "BelosSolverManager.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosPseudoBlockCGSolMgr.hpp"
#endif

#include "phist_gen_d.h"
#include "phist_belos_def.hpp"
