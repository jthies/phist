#include <cstdlib>
#include <cstring>
#include <cmath>

#include "phist_bgmres.h"

#include "ghost.h"

#include "phist_ScalarTraits.hpp"
#include "phist_rcp_helpers.hpp"
#include "phist_BelosOperatorTraits.hpp"

// Trilinos stuff
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_FancyOStream.hpp"

// GMRES solver manager from the Belos package
#include "BelosBlockGmresSolMgr.hpp"

#ifdef PHIST_KERNEL_LIB_GHOST
#include "Belos_GhostAdapter.hpp"
#elif defined(PHIST_KERNEL_LIB_EPETRA)
#include "Epetra_MultiVector.h"
#include "BelosEpetraAdapter.hpp"
#elif defined(PHIST_KERNEL_LIB_TPETRA)
#include "Tpetra_MultiVector.hpp"
#include "BelosTpetraAdapter.hpp"
#else
#error "bgmres only supports ghost, epetra and tpetra right now"
#endif

#include "phist_gen_d.h"
#include "phist_bgmres_def.hpp"

// disallow using Belos with other than double type for Epetra,
// that makes no sense because Epetra has only double as type.
#ifndef PHIST_KERNEL_LIB_EPETRA
#include "phist_gen_s.h"
#include "phist_bgmres_def.hpp"
#include "phist_gen_c.h"
#include "phist_bgmres_def.hpp"
#include "phist_gen_z.h"
#include "phist_bgmres_def.hpp"
#endif
