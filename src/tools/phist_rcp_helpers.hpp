#ifndef PHIST_RCP_HELPERS_HPP
#define PHIST_RCP_HELPERS_HPP

#include <mpi.h>
#include "Teuchos_RCP.hpp"

#ifdef PHIST_KERNEL_LIB_GHOST
#include "ghost.h"
#include "ghost/phist_GhostMV.hpp"
#endif

namespace phist 
  {
#ifdef PHIST_KERNEL_LIB_GHOST

  //! rcp for ghost_vec_t, includes creating the GhostMV wrapper
  Teuchos::RCP<GhostMV> rcp(ghost_vec_t* rawPtr, bool ownMem=true);

  //! specialization for const ghost_vec_t
  Teuchos::RCP<const GhostMV> rcp(const ghost_vec_t* rawPtr, bool ownMem=true);

#else
  //! standard rcp() function from Teuchos
  using Teuchos::rcp;
#endif

  //! get mvec pointer from reference
  template<typename MV>
  void* ref2ptr(MV& V)
    {
    return &V;
    }

  //! get mvec pointer from reference
  template<typename MV>
  const void* ref2ptr(const MV& V)
    {
    return &V;
    }

  } // namespace phist
#endif
