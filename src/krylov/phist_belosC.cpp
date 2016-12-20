#include "phist_config.h"

#ifdef PHIST_HAVE_SP

#include "phist_tools.h"
#include "phist_kernels.h"

#include <cstdlib>
#include <cstring>
#include <cmath>

#include "phist_belos.h"

#ifdef PHIST_HAVE_BELOS

// Trilinos stuff
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_FancyOStream.hpp"

#include "phist_rcp_helpers.hpp"

# ifdef PHIST_KERNEL_LIB_GHOST
#  include "ghost.h"
#  include "Belos_GhostAdapter.hpp"
# elif defined(PHIST_KERNEL_LIB_EPETRA)
#  include "Epetra_MultiVector.h"
#  include "BelosEpetraAdapter.hpp"
# elif defined(PHIST_KERNEL_LIB_TPETRA)
#  include "Tpetra_MultiVector.hpp"
#  include "BelosTpetraAdapter.hpp"
#  include "phist_tpetra_typedefs.hpp"
# else
// use general phist interface to Belos (may not be complete)
#  warning "Belos not supported for this kernel lib"
#  undef PHIST_HAVE_BELOS
# endif

#include "phist_BelosOperatorTraits.hpp"
#endif


// GMRES solver manager from the Belos package
#ifdef PHIST_HAVE_BELOS
#include "BelosSolverManager.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosPseudoBlockCGSolMgr.hpp"
#endif

#ifdef PHIST_KERNEL_LIB_EPETRA
#undef PHIST_HAVE_BELOS
#endif

#include "phist_gen_c.h"
#include "phist_belos_def.hpp"
#endif
