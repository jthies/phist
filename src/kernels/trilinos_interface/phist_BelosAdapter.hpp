#ifndef PHIST_BELOS_ADAPTER_HPP
#define PHIST_BELOS_ADAPTER_HPP

#include "phist_config.h"

#ifdef PHIST_HAVE_BELOS
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif


#include "phist_typedefs.h"
#include "phist_ScalarTraits.hpp"
#include "phist_tools.h"
#include "phist_macros.h"
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_Array.hpp>

#include <BelosConfigDefs.hpp>
#include <BelosTypes.hpp>
#include <BelosMultiVecTraits.hpp>

#include "phist_MultiVector.hpp"
#include "phist_rcp_helpers.hpp"

#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_BelosAdapter_def.hpp"

#include "phist_gen_c.h"
#include "phist_BelosAdapter_def.hpp"
#endif

#include "phist_gen_d.h"
#include "phist_BelosAdapter_def.hpp"

#include "phist_gen_z.h"
#include "phist_BelosAdapter_def.hpp"
#endif
#endif
