#include <cstdlib>
#include <cstring>
#include <cmath>

#include "phist_belos.h"

#include "phist_ScalarTraits.hpp"
#include "phist_rcp_helpers.hpp"
#include "phist_BelosOperatorTraits.hpp"

// Trilinos stuff
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_FancyOStream.hpp"

// GMRES solver manager from the Belos package
#include "BelosSolverManager.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosPseudoBlockCGSolMgr.hpp"

#ifdef PHIST_KERNEL_LIB_GHOST
#include "ghost.h"
#include "Belos_GhostAdapter.hpp"
#elif defined(PHIST_KERNEL_LIB_EPETRA)
#include "Epetra_MultiVector.h"
#include "BelosEpetraAdapter.hpp"
#elif defined(PHIST_KERNEL_LIB_TPETRA)
#include "Tpetra_MultiVector.hpp"
#include "BelosTpetraAdapter.hpp"
#elif defined(PHIST_KERNEL_LIB_FORTRAN)
#warning "belos not supported with fortran kernels lib"
#define NO_BELOS_IMPLEMENTATION
#else
#error "belos only supports ghost, epetra and tpetra right now"
#endif

#include "phist_gen_d.h"
#include "phist_belos_def.hpp"

#ifdef PHIST_KERNEL_LIB_EPETRA
#define NO_BELOS_IMPLEMENTATION
#endif

#include "phist_gen_s.h"
#include "phist_belos_def.hpp"
#include "phist_gen_c.h"
#include "phist_belos_def.hpp"
#include "phist_gen_z.h"
#include "phist_belos_def.hpp"

