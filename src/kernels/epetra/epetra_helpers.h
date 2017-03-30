/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef KERNELS_EPETRA_HELPERS_HPP
#define KERNELS_EPETRA_HELPERS_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

class Epetra_Operator;
class Epetra_CrsMatrix;
class Ifpack_Preconditioner;
namespace ML_Epetra
{
  class MultiLevelPreconditioner;
}

//! this file is just for internal use when implementing the
//! Epetra variant of our kernel functions.

typedef Teuchos::SerialDenseMatrix<int,double> Teuchos_sdMat_t;

class Epetra_MultiVector;

  //! create a Teuchos' view of a local mvec/sdMat
  static Teuchos::RCP<const Teuchos_sdMat_t > CreateTeuchosView(Teuchos::RCP<const Epetra_MultiVector> M, int* iflag);

  //! create a non-const Teuchos' view of a local mvec/sdMat
  static Teuchos::RCP<Teuchos_sdMat_t> CreateTeuchosViewNonConst(Teuchos::RCP<Epetra_MultiVector> M, int* iflag);

namespace phist 
{
  namespace internal 
  {

  //! class to allow computing A-simga*B or A-sigma*I and store it
  //! alongside a preconditioner which we can thus update during an
  //! eigensolve. All the ifpack and ML specific stuff is implemented
  //! in src/precon/tpl
  class prec_and_mat 
  {
    public:
    
    //!
    prec_and_mat(phist_Dconst_sparseMat_ptr A, double sigma,
               phist_Dconst_sparseMat_ptr B);
               
    //!
    ~prec_and_mat();

    //! construct A-sigma*B and store it in Mat
    void UpdateMatrix(phist_Dconst_sparseMat_ptr A, double sigma,
               phist_Dconst_sparseMat_ptr B, int *iflag);
    
    // pointer to preconditioner object
    Teuchos::RCP<Epetra_Operator> Prec;
    // for convenience points to Prec if this is Ifpack
    Teuchos::RCP<Ifpack_Preconditioner> IfpackPrec;
    // for convenience points to Prec if this is ML
    Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> MLPrec;
    // pointer to A-sigma*B
    Teuchos::RCP<Epetra_CrsMatrix> Mat;
  };
  }
}
  
#endif
