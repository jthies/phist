#include <cstdlib>
#include <cstring>
#include <cmath>

#include "phist_bgmres.h"

#include "ghost.h"

#include "phist_ScalarTraits.hpp"
#include "phist_BelosOperatorTraits.hpp"

// Trilinos stuff
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_FancyOStream.hpp"

// GMRES solver manager from the Belos package
#include "BelosBlockGmresSolMgr.hpp"

#ifdef PHIST_KERNEL_LIB_GHOST
#include "BelosGhostAdapter.hpp"
#elif defined(PHIST_KERNEL_LIB_TPETRA)
// classes for matrices and vectors
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "BelosTpetraAdapter.hpp"
#else
//TODO - get Epetra in the boat again
#error "bgmres only supports ghost and tpetra right now"
#endif

#include "phist_gen_s.h"
#include "phist_bgmres_def.hpp"
#include "phist_gen_d.h"
#include "phist_bgmres_def.hpp"
#include "phist_gen_c.h"
#include "phist_bgmres_def.hpp"
#include "phist_gen_z.h"
#include "phist_bgmres_def.hpp"

