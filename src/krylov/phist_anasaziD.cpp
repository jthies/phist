#include "phist_config.h"

/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <cstdlib>
#include <cstring>
#include <cmath>

#include "phist_anasazi.h"
#include "phist_kernels.h"
#include "phist_ScalarTraits.hpp"

#ifdef PHIST_HAVE_ANASAZI

#include "phist_rcp_helpers.hpp"
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
#ifdef BELOS_HAVE_TSQR
#warning "a bug in Trilinos 12.2.1 prevents compiling TraceMinDavidson if TSQR is enabled in Belos"
#endif

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
# ifndef OLD_TRILINOS
# include "AnasaziTraceMinDavidsonSolMgr.hpp"
# endif
#endif

#include "phist_gen_d.h"
#include "phist_anasazi_def.hpp"

