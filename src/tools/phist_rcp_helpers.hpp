#ifndef PHIST_RCP_HELPERS_HPP
#define PHIST_RCP_HELPERS_HPP

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
// we only need RCP's if we want to interact
// with the iterative solvers in Belos
#ifdef PHIST_HAVE_TEUCHOS
#include "Teuchos_RCP.hpp"

#ifdef PHIST_KERNEL_LIB_GHOST
#include "ghost.h"
#include "ghost/phist_GhostMV.hpp"
#endif

#include "trilinos_interface/phist_MultiVector.hpp"
#include "phist_void_aliases.h"

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

#elif defined(PHIST_KERNEL_LIB_TPETRA)||defined(PHIST_KERNEL_LIB_EPETRA)
// for any C++ based libraries we can just use Teuchos::rcp

  //! standard rcp() function from Teuchos
  using Teuchos::rcp;
#endif

#ifdef PHIST_HAVE_SP
  //! rcp for Smvec_t
  Teuchos::RCP<MultiVector<float>> Srcp(Smvec_ptr_t rawPtr, 
        bool ownMem=true);

  //! rcp for const Smvec_t
  Teuchos::RCP<const MultiVector<float> > Srcp(Sconst_mvec_ptr_t rawPtr, 
        bool ownMem=true);

  //! rcp for Cmvec_t
  Teuchos::RCP<MultiVector<s_complex_t> > Crcp(Cmvec_ptr_t rawPtr, 
        bool ownMem=true);

  //! rcp for const Cmvec_t
  Teuchos::RCP<const MultiVector<s_complex_t> > Crcp(Cconst_mvec_ptr_t rawPtr, 
        bool ownMem=true);
#endif

  //! rcp for Dmvec_t
  Teuchos::RCP<MultiVector<double> > Drcp(Dmvec_ptr_t rawPtr, 
        bool ownMem=true);

  //! rcp for const Dmvec_t
  Teuchos::RCP<const MultiVector<double> > Drcp(Dconst_mvec_ptr_t rawPtr, 
        bool ownMem=true);

  //! rcp for Zmvec_t
  Teuchos::RCP<MultiVector<d_complex_t> > Zrcp(Zmvec_ptr_t rawPtr, 
        bool ownMem=true);

  //! rcp for const Zmvec_t
  Teuchos::RCP<const MultiVector<d_complex_t> > Zrcp(Zconst_mvec_ptr_t rawPtr, 
        bool ownMem=true);


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
#endif

  template<typename ST>
  void* ref2ptr(MultiVector<ST>& V);
  
  template<typename ST>
  const void* ref2ptr(const MultiVector<ST>& V);

  } // namespace phist
#endif /* PHIST_HAVE_TEUCHOS */
#endif
