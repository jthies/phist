/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef KERNELS_BUILTIN_TYPEDEFS_HPP
#define KERNELS_BUILTIN_TYPEDEFS_HPP

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
#ifdef PHIST_HAVE_KOKKOS__disabled
#include "Kokkos_DefaultNode.hpp"
#include "Kokkos_DefaultKernels.hpp"

typedef Kokkos::DefaultNode::DefaultNodeType node_type;
#endif

#include "phist_typedefs.h"
#include "phist_ScalarTraits.hpp"


template <typename ST>
class Traits
{

public:
  
  //! multi vectors
  typedef void mvec_t;

  //! serial dense matrix
  typedef void sdMat_t;

#ifdef PHIST_HAVE_TEUCHOS
  //! serial dense matrix from Teuchos, we need this for e.g. the BLAS interface.
  typedef Teuchos::SerialDenseMatrix<phist_lidx,ST> Teuchos_sdMat_t;
#endif

  //! CRS matrices
  typedef void sparseMat_t;

  //! create a Teuchos' view of a local sdMat
  //static Teuchos::RCP<const Teuchos_sdMat_t> CreateTeuchosView(Teuchos::RCP<const sdMat_t> M, int* iflag)
  //{
    //*iflag=0;

    //const ST *M_val = NULL;
    //phist_lidx lda;
    //int nrows, ncols;
    //PHIST_CHK_IERR(SUBR(sdMat_extract_view)(*M,&M_val,lda,iflag),*iflag);
    //PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(*M,&nrows,iflag),*iflag);
    //PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(*M,&ncols,iflag),*iflag);

    //Teuchos::RCP<const Teuchos_sdMat_t> M_view
      //= Teuchos::rcp(new Teuchos_sdMat_t(Teuchos::View,M_val,lda,nrows,ncols));
    //return M_view;
  //}

  //! create a non-const Teuchos' view of a local mvec/sdMat
  //static Teuchos::RCP<Teuchos_sdMat_t> CreateTeuchosViewNonConst(Teuchos::RCP<sdMat_t> M, int* iflag)
  //{
    //*iflag=0;

    //ST *M_val = NULL;
    //phist_lidx lda;
    //int nrows, ncols;
    //PHIST_CHK_IERR(SUBR(sdMat_extract_view)(*M,&M_val,lda,iflag),*iflag);
    //PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(*M,&nrows,iflag),*iflag);
    //PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(*M,&ncols,iflag),*iflag);

    //Teuchos::RCP<Teuchos_sdMat_t> M_view
      //= Teuchos::rcp(new Teuchos_sdMat_t(Teuchos::View,M_val,lda,nrows,ncols));
    //return M_view;
  //}

};
  
#endif
