#ifndef ANASAZI_PHIST_ADAPTER_HPP
#define ANASAZI_PHIST_ADAPTER_HPP

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#if defined(PHIST_HAVE_ANASAZI)&&defined(PHIST_HAVE_BELOS)

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

#include "phist_MultiVector.hpp"
#include "phist_rcp_helpers.hpp"

#include "phist_BelosAdapter.hpp"

using namespace phist;

namespace Anasazi
{

  template<typename ST>
  class MultiVecTraits<ST,MultiVector<ST> >
        : public Belos::MultiVecTraits<ST,MultiVector<ST> >
  {
  };
}// namespace Anasazi


#endif
#endif
