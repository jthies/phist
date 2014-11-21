#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "phist_rcp_helpers.hpp"


namespace phist
  {
#ifdef PHIST_HAVE_TEUCHOS

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

  //@}


#endif

// special functions for the general phist/belos interface,
// to avoid the problem that all
// our pointers are void
#ifdef HAVE_SP
  //
  Teuchos::RCP<MultiVector<float> > Srcp(Smvec_ptr_t rawPtr, bool ownMem)
    {
    // if the RCP should own the memory of the vector, it owns the wrapper as well.
    // If it does not own the memory of the vector, it may still destroy the wrapper
    // if it is no longer needed.
    return Teuchos::rcp(new MultiVector<float>(rawPtr,ownMem),true);
    }

  //
  Teuchos::RCP<const MultiVector<float> > Srcp(Sconst_mvec_ptr_t rawPtr, bool ownMem)
    {
    return Teuchos::rcp(new MultiVector<float>(const_cast<Smvec_ptr_t>(rawPtr),ownMem),true);
    }

  //
  Teuchos::RCP<MultiVector<s_complex_t> > Crcp(Cmvec_ptr_t rawPtr, bool ownMem)
  {
    // if the RCP should own the memory of the vector, it owns the wrapper as well.
    // If it does not own the memory of the vector, it may still destroy the wrapper
    // if it is no longer needed.
    return Teuchos::rcp(new MultiVector<s_complex_t>(rawPtr,ownMem),true);
  }

  // specialization for const ghost_densemat_t
  Teuchos::RCP<const MultiVector<s_complex_t> > Crcp(Cconst_mvec_ptr_t rawPtr, bool ownMem)
  {
    return Teuchos::rcp(new MultiVector<s_complex_t>(const_cast<Cmvec_ptr_t>(rawPtr),ownMem),true);
  }

#endif


  //
  Teuchos::RCP<MultiVector<double> > Drcp(Dmvec_ptr_t rawPtr, bool ownMem)
  {
    // if the RCP should own the memory of the vector, it owns the wrapper as well.
    // If it does not own the memory of the vector, it may still destroy the wrapper
    // if it is no longer needed.
    return Teuchos::rcp(new MultiVector<double>(rawPtr,ownMem),true);
  }

  //
  Teuchos::RCP<const MultiVector<double> > Drcp(Dconst_mvec_ptr_t rawPtr, bool ownMem)
  {
    return Teuchos::rcp(new MultiVector<double>(const_cast<Dmvec_ptr_t>(rawPtr),ownMem),true);
  }

  //
  Teuchos::RCP<MultiVector<d_complex_t> > Zrcp(Zmvec_ptr_t rawPtr, bool ownMem)
  {
    // if the RCP should own the memory of the vector, it owns the wrapper as well.
    // If it does not own the memory of the vector, it may still destroy the wrapper
    // if it is no longer needed.
    return Teuchos::rcp(new MultiVector<d_complex_t>(rawPtr,ownMem),true);
  }

  // specialization for const ghost_densemat_t
  Teuchos::RCP<const MultiVector<d_complex_t> > Zrcp(Zconst_mvec_ptr_t rawPtr, bool ownMem)
  {
    return Teuchos::rcp(new MultiVector<d_complex_t>(const_cast<Zmvec_ptr_t>(rawPtr),ownMem),true);
  }

#endif
  } // namespace phist
