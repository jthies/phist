#include "phist_config.h"
#include "phist_rcp_helpers.hpp"


namespace phist 
  {
#ifdef PHIST_KERNEL_LIB_GHOST

  // rcp for ghost_densemat_t, includes creating the GhostMV wrapper
  Teuchos::RCP<GhostMV> rcp(ghost_densemat_t* rawPtr, bool ownMem)
    {
    // if the RCP should own the memory of the vector, it owns the wrapper as well.
    // If it does not own the memory of the vector, it may still destroy the wrapper
    // if it is no longer needed.
    return Teuchos::rcp(new GhostMV(rawPtr,ownMem),true);
    }

  // specialization for const ghost_densemat_t
  Teuchos::RCP<const GhostMV> rcp(const ghost_densemat_t* rawPtr, bool ownMem)
    {
    // if the RCP should own the memory of the vector, it owns the wrapper as well.
    // If it does not own the memory of the vector, it may still destroy the wrapper
    // if it is no longer needed.
    return Teuchos::rcp(new GhostMV(const_cast<ghost_densemat_t*>(rawPtr),ownMem),true);
    }

  //!\name ref2ptr specialization for Ghost: input GhostMV, output ghost_densemat_t*
  //@{

  // get mvec pointer from reference (implementation for GhostMV that returns
  // ghost_densemat_t*)
  template<>
  void* ref2ptr(GhostMV& V)
    {
    return (void*)V.get();
    }

  // get const mvec pointer from const reference (implementation for GhostMV that returns
  // const ghost_densemat_t*)
  template<>
  const void* ref2ptr(const GhostMV& V)
    {
    return (const void*)V.get();
    }
  //@}

#endif
  } // namespace phist
