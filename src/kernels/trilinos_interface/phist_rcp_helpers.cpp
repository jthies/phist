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
  // rcp for ghost_densemat, includes creating the GhostMV wrapper
  Teuchos::RCP<GhostMV> rcp(ghost_densemat* rawPtr, bool ownMem)
    {
    // if the RCP should own the memory of the vector, it owns the wrapper as well.
    // If it does not own the memory of the vector, it may still destroy the wrapper
    // if it is no longer needed.
    return Teuchos::rcp(new GhostMV(rawPtr,ownMem),true);
    }

  // specialization for const ghost_densemat
  Teuchos::RCP<const GhostMV> rcp(const ghost_densemat* rawPtr, bool ownMem)
    {
    // if the RCP should own the memory of the vector, it owns the wrapper as well.
    // If it does not own the memory of the vector, it may still destroy the wrapper
    // if it is no longer needed.
    return Teuchos::rcp(new GhostMV(const_cast<ghost_densemat*>(rawPtr),ownMem),true);
    }

  //!\name ref2ptr specialization for Ghost: input GhostMV, output ghost_densemat*
  //@{

  // get mvec pointer from reference (implementation for GhostMV that returns
  // ghost_densemat*)
  template<>
  void* ref2ptr(GhostMV& V)
    {
    return (void*)V.get();
    }

  //!
  template<>
  void const* ref2ptr(GhostMV const& V)
    {
    return (void const*)V.get();
    }

  //@}


#elif !defined(PHIST_KERNEL_LIB_EPETRA) && !defined(PHIST_KERNEL_LIB_TPETRA)

// special functions for the general phist/belos interface,
// to avoid the problem that all
// our pointers are void
#ifdef HAVE_SP
  //
  Teuchos::RCP<MultiVector<float> > Srcp(Smvec_ptr rawPtr, bool ownMem)
    {
    // if the RCP should own the memory of the vector, it owns the wrapper as well.
    // If it does not own the memory of the vector, it may still destroy the wrapper
    // if it is no longer needed.
    return Teuchos::rcp(new MultiVector<float>(rawPtr,ownMem),true);
    }

  //
  Teuchos::RCP<const MultiVector<float> > Srcp(Sconst_mvec_ptr rawPtr, bool ownMem)
    {
    return Teuchos::rcp(new MultiVector<float>(const_cast<Smvec_ptr>(rawPtr),ownMem),true);
    }

  //
  Teuchos::RCP<MultiVector<phist_s_complex> > Crcp(Cmvec_ptr rawPtr, bool ownMem)
  {
    // if the RCP should own the memory of the vector, it owns the wrapper as well.
    // If it does not own the memory of the vector, it may still destroy the wrapper
    // if it is no longer needed.
    return Teuchos::rcp(new MultiVector<phist_s_complex>(rawPtr,ownMem),true);
  }

  // specialization for const ghost_densemat
  Teuchos::RCP<const MultiVector<phist_s_complex> > Crcp(Cconst_mvec_ptr rawPtr, bool ownMem)
  {
    return Teuchos::rcp(new MultiVector<phist_s_complex>(const_cast<Cmvec_ptr>(rawPtr),ownMem),true);
  }

#endif


  //
  Teuchos::RCP<MultiVector<double> > Drcp(Dmvec_ptr rawPtr, bool ownMem)
  {
    // if the RCP should own the memory of the vector, it owns the wrapper as well.
    // If it does not own the memory of the vector, it may still destroy the wrapper
    // if it is no longer needed.
    return Teuchos::rcp(new MultiVector<double>(rawPtr,ownMem),true);
  }

  //
  Teuchos::RCP<const MultiVector<double> > Drcp(Dconst_mvec_ptr rawPtr, bool ownMem)
  {
    return Teuchos::rcp(new MultiVector<double>(const_cast<Dmvec_ptr>(rawPtr),ownMem),true);
  }

  //
  Teuchos::RCP<MultiVector<phist_d_complex> > Zrcp(Zmvec_ptr rawPtr, bool ownMem)
  {
    // if the RCP should own the memory of the vector, it owns the wrapper as well.
    // If it does not own the memory of the vector, it may still destroy the wrapper
    // if it is no longer needed.
    return Teuchos::rcp(new MultiVector<phist_d_complex>(rawPtr,ownMem),true);
  }

  // specialization for const ghost_densemat
  Teuchos::RCP<const MultiVector<phist_d_complex> > Zrcp(Zconst_mvec_ptr rawPtr, bool ownMem)
  {
    return Teuchos::rcp(new MultiVector<phist_d_complex>(const_cast<Zmvec_ptr>(rawPtr),ownMem),true);
  }

#endif
#endif
  } // namespace phist
