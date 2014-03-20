#ifndef PHIST_RCP_HELPERS_HPP
#define PHIST_RCP_HELPERS_HPP

#include <mpi.h>
#include "Teuchos_RCP.hpp"

#ifdef PHIST_KERNEL_LIB_GHOST
#include "ghost.h"
#include "ghost/phist_GhostMV.hpp"
#endif

//! in this file we define various helper functions to cope with
//! the problem that ghost vectors (classical C structs) do not 
//! have a destructor and so the Trilinos 'reference counting pointer'
//! (RCP) concept does not work properly. In order for memory management
//! to work correctly, the phist::rcp function should be used instead
//| of Teuchos::rcp to create an RCP from a raw pointer (for C++ vector
//! classes it will just default to Teuchos::rcp). When getting a reference
//! to an object the ref2ptr function should be used (rather than the ampersand &)
//! to get the raw pointer. This is because our GhostMV wrapper is not a ghost_densemat_t
//! but will return one by the get() function. Again, for other vector classes than
//! ghost_densemat_t ref2ptr(V) just defaults to &V.

namespace phist 
  {
#ifdef PHIST_KERNEL_LIB_GHOST


  //! rcp for ghost_densemat_t, includes creating the GhostMV wrapper
  Teuchos::RCP<GhostMV> rcp(ghost_densemat_t* rawPtr, bool ownMem=true);

  //! specialization for const ghost_densemat_t
  Teuchos::RCP<const GhostMV> rcp(const ghost_densemat_t* rawPtr, bool ownMem=true);

#elif defined(PHIST_KERNEL_LIB_FORTRAN)

  //! rcp for fortran mvec, includes creating the FortranMV wrapper
  //Teuchos::RCP<FortranMV> rcp(void* rawPtr, bool ownMem=true);

  ////! specialization for const FortranMV
  //Teuchos::RCP<const FortranMV> rcp(const void* rawPtr, bool ownMem=true);


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

#ifdef PHIST_KERNEL_LIB_GHOST
  
  template<>
  void* ref2ptr(GhostMV& V);

  template<>
  const void* ref2ptr(const GhostMV& V);

#elif defined(PHIST_KERNEL_LIB_FORTRAN)

  //template<>
  //void* ref2ptr(FortranMV& V);

  //template<>
  //const void* ref2ptr(const FortranMV& V);
#endif
  } // namespace phist
#endif
