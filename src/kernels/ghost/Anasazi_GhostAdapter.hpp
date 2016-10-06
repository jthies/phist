#ifndef ANASAZI_GHOST_ADAPTER_HPP
#define ANASAZI_GHOST_ADAPTER_HPP

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#if defined(PHIST_HAVE_ANASAZI)&&defined(PHIST_HAVE_BELOS)

#include <ghost.h>
#include "phist_typedefs.h"
#include "phist_ScalarTraits.hpp"
#include "./typedefs.hpp"
#include "phist_tools.h"
#include "phist_macros.h"

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_Array.hpp>

#include <AnasaziConfigDefs.hpp>
#include <AnasaziTypes.hpp>
#include <AnasaziMultiVecTraits.hpp>

#include "phist_GhostMV.hpp"
#include "phist_rcp_helpers.hpp"

#ifdef HAVE_ANASAZI_TSQR
#  include <Ghost_TsqrAdaptor.hpp>
#endif // HAVE_ANASAZI_TSQR

#include "Belos_GhostAdapter.hpp"

using namespace phist;

namespace Anasazi 
{

  template<typename ST>
  class MultiVecTraits<ST,GhostMV>: public Belos::MultiVecTraits<ST,GhostMV>
  {
  };
}// namespace Anasazi


#endif
#endif
