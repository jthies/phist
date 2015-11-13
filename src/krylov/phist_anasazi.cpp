#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <cstdlib>
#include <cstring>
#include <cmath>

#include "phist_anasazi.h"

#include "phist_ScalarTraits.hpp"

#ifdef PHIST_HAVE_ANASAZI

#include "phist_rcp_helpers.hpp"
#include "phist_AnasaziOperatorTraits.hpp"

// Trilinos stuff
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_FancyOStream.hpp"

# ifdef PHIST_KERNEL_LIB_GHOST
#  include "ghost.h"
#  include "Anasazi_GhostAdapter.hpp"
# elif defined(PHIST_KERNEL_LIB_EPETRA)
#  include "Epetra_MultiVector.h"
#  include "AnasaziEpetraAdapter.hpp"
# elif defined(PHIST_KERNEL_LIB_TPETRA)
#  include "Tpetra_MultiVector.hpp"
#  include "AnasaziTpetraAdapter.hpp"
# else
// use general phist/anasazi interface
#include "phist_MultiVector.hpp"
#include "phist_AnasaziAdapter.hpp"
# endif
#endif

// Block Krylov-Schur solver manager from the Anasazi package
#ifdef PHIST_HAVE_ANASAZI
// adaptation of a basic ortho class from Trili 11.12 to avoid 
// col-wise norm calculations
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziSolverManager.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#endif

#include "phist_gen_d.h"
#include "phist_anasazi_def.hpp"

#ifdef PHIST_KERNEL_LIB_EPETRA
#undef PHIST_HAVE_ANASAZI
#endif
#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_anasazi_def.hpp"
#include "phist_gen_c.h"
#include "phist_anasazi_def.hpp"
#endif
#include "phist_gen_z.h"
#include "phist_anasazi_def.hpp"

