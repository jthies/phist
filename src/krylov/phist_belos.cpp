#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <cstdlib>
#include <cstring>
#include <cmath>

#include "phist_belos.h"

#include "phist_ScalarTraits.hpp"

#ifdef PHIST_HAVE_BELOS

#include "phist_rcp_helpers.hpp"

// Trilinos stuff
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_FancyOStream.hpp"

# ifdef PHIST_KERNEL_LIB_GHOST
#  include "ghost.h"
#  include "Belos_GhostAdapter.hpp"
# elif defined(PHIST_KERNEL_LIB_EPETRA)
#  include "Epetra_MultiVector.h"
#  include "BelosEpetraAdapter.hpp"
# elif defined(PHIST_KERNEL_LIB_TPETRA)
#  include "Tpetra_MultiVector.hpp"
#  include "BelosTpetraAdapter.hpp"
# else
#  warning "belos only supported with ghost, epetra and tpetra right now"
#  undef PHIST_HAVE_BELOS
# endif
#endif

#ifdef PHIST_HAVE_BELOS
#include "phist_BelosOperatorTraits.hpp"

// GMRES solver manager from the Belos package
#include "BelosSolverManager.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosPseudoBlockCGSolMgr.hpp"
#endif

#include "phist_gen_d.h"
#include "phist_belos_def.hpp"

#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_belos_def.hpp"
#include "phist_gen_c.h"
#include "phist_belos_def.hpp"
#endif
#include "phist_gen_z.h"
#include "phist_belos_def.hpp"

