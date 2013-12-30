#ifndef KERNELS_FORTRAN_TYPEDEFS_HPP
#define KERNELS_FORTRAN_TYPEDEFS_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_DataAccess.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include "Kokkos_DefaultNode.hpp"
#include "Kokkos_DefaultKernels.hpp"

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
  typedef void mvec_t;

  //! serial dense matrix
  typedef void sdMat_t;

  //! serial dense matrix from Teuchos, we need this for e.g. the BLAS interface.
  typedef Teuchos::SerialDenseMatrix<lidx_t,ST> Teuchos_sdMat_t;

  //! CRS matrices
  typedef void crsMat_t;

  //! create a Teuchos' view of a local sdMat
  //static Teuchos::RCP<const Teuchos_sdMat_t> CreateTeuchosView(Teuchos::RCP<const sdMat_t> M, int* ierr)
  //{
    //*ierr=0;

    //const ST *M_val = NULL;
    //lidx_t lda;
    //int nrows, ncols;
    //PHIST_CHK_IERR(SUBR(sdMat_extract_view)(*M,&M_val,lda,ierr),*ierr);
    //PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(*M,&nrows,ierr),*ierr);
    //PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(*M,&ncols,ierr),*ierr);

    //Teuchos::RCP<const Teuchos_sdMat_t> M_view
      //= Teuchos::rcp(new Teuchos_sdMat_t(Teuchos::View,M_val,lda,nrows,ncols));
    //return M_view;
  //}

  //! create a non-const Teuchos' view of a local mvec/sdMat
  //static Teuchos::RCP<Teuchos_sdMat_t> CreateTeuchosViewNonConst(Teuchos::RCP<sdMat_t> M, int* ierr)
  //{
    //*ierr=0;

    //ST *M_val = NULL;
    //lidx_t lda;
    //int nrows, ncols;
    //PHIST_CHK_IERR(SUBR(sdMat_extract_view)(*M,&M_val,lda,ierr),*ierr);
    //PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(*M,&nrows,ierr),*ierr);
    //PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(*M,&ncols,ierr),*ierr);

    //Teuchos::RCP<Teuchos_sdMat_t> M_view
      //= Teuchos::rcp(new Teuchos_sdMat_t(Teuchos::View,M_val,lda,nrows,ncols));
    //return M_view;
  //}

};
  
#endif
