#ifndef KERNELS_EPETRA_HELPERS_HPP
#define KERNELS_EPETRA_HELPERS_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

//! this file is just for internal use when implementing the
//! Epetra variant of our kernel functions.

typedef Teuchos::SerialDenseMatrix<int,double> Teuchos_sdMat_t;

class Epetra_MultiVector;

  //! create a Teuchos' view of a local mvec/sdMat
  static Teuchos::RCP<const Teuchos_sdMat_t > CreateTeuchosView(Teuchos::RCP<const Epetra_MultiVector> M, int* ierr);

  //! create a non-const Teuchos' view of a local mvec/sdMat
  static Teuchos::RCP<Teuchos_sdMat_t> CreateTeuchosViewNonConst(Teuchos::RCP<Epetra_MultiVector> M, int* ierr);
  
#endif
