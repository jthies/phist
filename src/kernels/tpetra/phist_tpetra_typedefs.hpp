#ifndef KERNELS_TPETRA_TYPEDEFS_HPP
#define KERNELS_TPETRA_TYPEDEFS_HPP

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "phist_typedefs.h"

#include <Tpetra_Map_decl.hpp>
#include "Teuchos_RCP.hpp"
#include "Teuchos_Comm.hpp"

namespace phist
{
namespace tpetra
{
using map_type = Tpetra::Map<phist_lidx, phist_gidx>;
using comm_type =  Teuchos::Comm<int>;


} //namespace tpetra
} //namespace phist
#endif