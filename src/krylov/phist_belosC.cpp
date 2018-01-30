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

#if defined(PHIST_KERNEL_LIB_TPETRA)
#include "TpetraCore_config.h"
# ifndef HAVE_TPETRA_INST_COMPLEX_FLOAT
# undef PHIST_HAVE_BELOS
# endif
#endif

#include "phist_tools.h"
#include "phist_kernels.h"
#include "phist_operator.h"
#include "phist_ScalarTraits.hpp"
#include "phist_belos.h"

#ifdef PHIST_HAVE_BELOS

// Trilinos stuff
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_FancyOStream.hpp"

#include "Belos_PhistAdapter.hpp"
#include "phist_BelosOperatorTraits.hpp"

// we provide our own ICGS/CholQR ortho manager based on the 'orthog' routine.
// Since there is no easy mechanism to extend Belos with new OrthoManagers,
// we do this by specializing their class for our own MV and OP types.
#include "BelosICGSOrthoManager.hpp"
#include "phist_BelosICGSOrthoManager.hpp"

#include "BelosSolverManager.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosPseudoBlockCGSolMgr.hpp"
#endif


#include "phist_gen_c.h"
#include "phist_belos_def.hpp"
#endif
