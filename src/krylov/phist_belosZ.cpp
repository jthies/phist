/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_config.h"

#if defined(PHIST_HAVE_CMPLX)

/* check wether we have to disable Anasazi because the Trilinos installation doesn't support this data type */
#if defined(PHIST_KERNEL_LIB_TPETRA)
#include "TpetraCore_config.h"
# ifndef HAVE_TPETRA_INST_COMPLEX_DOUBLE
# undef PHIST_HAVE_BELOS
# endif
#elif defined(PHIST_HAVE_TEUCHOS)
# include "Teuchos_config.h"
# ifndef HAVE_TEUCHOS_COMPLEX
# undef PHIST_HAVE_BELOS
# endif
#endif

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

#include "phist_gen_z.h"
#include "phist_belos_def.hpp"

#endif
