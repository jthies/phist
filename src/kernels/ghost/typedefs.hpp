#ifndef KERNELS_GHOST_TYPEDEFS_HPP
#define KERNELS_GHOST_TYPEDEFS_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_DataAccess.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include "Kokkos_DefaultNode.hpp"
#include "Kokkos_DefaultKernels.hpp"

#include "ghost.h"

#include "phist_typedefs.h"
#include "phist_ScalarTraits.hpp"

typedef Kokkos::DefaultNode::DefaultNodeType node_t;

template <typename ST>
class Traits
  {

public:
  
  //!
  typedef typename Kokkos::DefaultKernels<ST,lidx_t,node_t>::SparseOps localOps_t;

  //! multi vectors
  typedef ghost_vec_t mvec_t;

  //! serial dense matrix - just a multivector with a serial map.
  typedef ghost_vec_t sdMat_t;

  //! serial dense matrix from Teuchos, we need this for e.g. the BLAS interface.
  typedef Teuchos::SerialDenseMatrix<lidx_t,ST> Teuchos_sdMat_t;

  //! CRS matrices
  typedef ghost_mat_t crsMat_t;

  //! create a Teuchos' view of a local mvec/sdMat
  static Teuchos::RCP<const Teuchos_sdMat_t> CreateTeuchosView(Teuchos::RCP<const sdMat_t> M, int* ierr)
    {
    *ierr=0;
    lidx_t stride = M->traits->nrowspadded;
    lidx_t nrows = M->traits->nrows;
    lidx_t ncols = M->traits->nvecs;

    if (M->traits->datatype != phist::ScalarTraits<ST>::ghost_dt)
      {
      *ierr=-1;
      return Teuchos::null;
      }
    
    const ST *M_val = (const ST*)M->val;
    Teuchos::RCP<const Teuchos_sdMat_t> M_view
                  = Teuchos::rcp(new Teuchos_sdMat_t(Teuchos::View,M_val,stride,nrows,ncols));
    return M_view;     
    }

  //! create a non-const Teuchos' view of a local mvec/sdMat
  static Teuchos::RCP<Teuchos_sdMat_t> CreateTeuchosViewNonConst(Teuchos::RCP<sdMat_t> M, int* ierr)
    {
    *ierr=0;
    
    if (M->traits->flags & GHOST_VEC_SCATTERED)
      {
      *ierr=-1;
      PHIST_OUT(PHIST_ERROR,"to create a Teuchos view we need constant stride in ghost_vec_t");
      return Teuchos::null;
      }

    int stride = M->traits->nrowspadded;
    int nrows = M->traits->nrows;
    int ncols = M->traits->nvecs;
    
    if (M->traits->datatype != phist::ScalarTraits<ST>::ghost_dt)
      {
      *ierr=-2;
      PHIST_OUT(PHIST_ERROR,"incorrect data type in ghost_vec_t");
      return Teuchos::null;
      }

    ST *M_val = (ST*)M->val[0];
    
    Teuchos::RCP<Teuchos_sdMat_t> M_view
                  = Teuchos::rcp(new Teuchos_sdMat_t(Teuchos::View,M_val,stride,nrows,ncols));
    return M_view;                  
    }

  };
  
#endif
