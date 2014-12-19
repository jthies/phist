#ifndef KERNELS_GHOST_TYPEDEFS_HPP
#define KERNELS_GHOST_TYPEDEFS_HPP

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#ifdef PHIST_HAVE_TEUCHOS
#include "Teuchos_RCP.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_DataAccess.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#endif

#include "ghost.h"

#include "phist_typedefs.h"
#include "phist_ScalarTraits.hpp"

/* Kokkos is the Trilinos package for node-optimized kernel routines.
   we use it in conjunction with TSQR, the robust and communication-efficient
   'tall skinny matrix' QR factoriztion.
   
   While Kokkos supports a large number of node types (serial, CUDA, OpenMP, Threadpool 
   (aka TPI), Intel TBB...),
   TSQR only supports TBB and TPI (as of version 11.12.2).
*/
#ifdef PHIST_HAVE_KOKKOS
#include "KokkosClassic_config.h"
/*typedef Kokkos::DefaultNode::DefaultNodeType node_t;*/
# if defined(HAVE_KOKKOSCLASSIC_THREADPOOL)
# include "Kokkos_TPINode.hpp"
typedef Kokkos::TPINode node_t;
# elif defined(HAVE_KOKKOSCLASSIC_TBB)
# include "Kokkos_TBBNode.hpp"
typedef Kokkos::TbbNode node_t;
# else
# warning "Your Trilinos installation does not support ThreadPool or Intel TBB, so we use the serial node for the TSQR interface!"
# include "Kokkos_SerialNode.hpp"
typedef Kokkos::SerialNode node_t;
# endif
#endif
template <typename ST>
class Traits
  {

public:

#ifdef PHIST_HAVE_TEUCHOS
  //! serial dense matrix from Teuchos, we need this for e.g. the BLAS interface.
  //! Note: the index type *must* be int here, not int64_t, so we decided to have
  //! phist local indices ints, even if ghost uses int64_t.
  typedef Teuchos::SerialDenseMatrix<int,ST> Teuchos_sdMat_t;

#endif

  //! multi vectors
  typedef ghost_densemat_t mvec_t;

  //! serial dense matrix - just a multivector with a serial map.
  typedef ghost_densemat_t sdMat_t;

  //! CRS matrices
  typedef ghost_sparsemat_t crsMat_t;

#ifdef PHIST_HAVE_TEUCHOS
  //! create a Teuchos' view of a local mvec/sdMat
  static Teuchos::RCP<const Teuchos_sdMat_t> CreateTeuchosView(Teuchos::RCP<const sdMat_t> M, int* iflag)
    {
    *iflag=0;
    lidx_t stride = M->traits.nrowspadded;
    lidx_t nrows = M->traits.nrows;
    lidx_t ncols = M->traits.ncols;

    if (M->traits.datatype != phist::ScalarTraits<ST>::ghost_dt)
      {
      *iflag=-1;
      return Teuchos::null;
      }
    
    const ST *M_val = (const ST*)M->val;
    Teuchos::RCP<const Teuchos_sdMat_t> M_view
                  = Teuchos::rcp(new Teuchos_sdMat_t(Teuchos::View,M_val,stride,nrows,ncols));
    return M_view;     
    }

  //! create a non-const Teuchos' view of a local mvec/sdMat
  static Teuchos::RCP<Teuchos_sdMat_t> CreateTeuchosViewNonConst(Teuchos::RCP<sdMat_t> M, int* iflag)
    {
    *iflag=0;
    
    if (M->traits.flags & GHOST_DENSEMAT_SCATTERED)
      {
      *iflag=-1;
      PHIST_OUT(PHIST_ERROR,"to create a Teuchos view we need constant stride in ghost_vec_t");
      return Teuchos::null;
      }

    int stride = M->traits.nrowspadded;
    int nrows = M->traits.nrows;
    int ncols = M->traits.ncols;
    
    if (M->traits.datatype != phist::ScalarTraits<ST>::ghost_dt)
      {
      *iflag=-2;
      PHIST_OUT(PHIST_ERROR,"incorrect data type in ghost_vec_t");
      return Teuchos::null;
      }

    ST *M_val;
    ghost_densemat_valptr(M.get(),(void**)(&M_val));
    
    Teuchos::RCP<Teuchos_sdMat_t> M_view
                  = Teuchos::rcp(new Teuchos_sdMat_t(Teuchos::View,M_val,stride,nrows,ncols));
    return M_view;                  
    }
#endif
  };
  
#endif
