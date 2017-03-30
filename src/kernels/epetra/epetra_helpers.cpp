/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "./epetra_helpers.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_DataAccess.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "EpetraExt_MatrixMatrix.h"

#include "Epetra_MultiVector.h"

  //! create a Teuchos' view of a local mvec/sdMat
  static Teuchos::RCP<const Teuchos_sdMat_t > CreateTeuchosView(Teuchos::RCP<const Epetra_MultiVector> M, int* iflag)
  {
    *iflag=0;
    int stride = M->Stride();
    int nrows = M->MyLength();
    int ncols = M->NumVectors();
    int lda; 
    double *M_val;
    *iflag=M->ExtractView(&M_val,&lda);
    if (*iflag<0) return Teuchos::null;
    Teuchos::RCP<const Teuchos_sdMat_t> M_view
                  = Teuchos::rcp(new Teuchos_sdMat_t(Teuchos::View,M_val,stride,nrows,ncols));
    return M_view;
  }

  //! create a non-const Teuchos' view of a local mvec/sdMat
  static Teuchos::RCP<Teuchos_sdMat_t> CreateTeuchosViewNonConst(Teuchos::RCP<Epetra_MultiVector> M, int* iflag)
  {
    *iflag=0;
    int stride = M->Stride();
    int nrows = M->MyLength();
    int ncols = M->NumVectors();
    int lda; 
    double *M_val;
    *iflag=M->ExtractView(&M_val,&lda);
    if (*iflag<0) return Teuchos::null;
    Teuchos::RCP<Teuchos_sdMat_t> M_view
                  = Teuchos::rcp(new Teuchos_sdMat_t(Teuchos::View,M_val,stride,nrows,ncols));
    return M_view;
  }

  
namespace phist 
{
  namespace internal 
  {

    // class to allow computing A-simga*B or A-sigma*I and store it
    // alongside a preconditioner which we can thus update during an
    // eigensolve. All the ifpack and ML specific stuff is implemented
    // in src/precon/tpl
    prec_and_mat::prec_and_mat(phist_Dconst_sparseMat_ptr A, double sigma,
               phist_Dconst_sparseMat_ptr B)
      : Prec(Teuchos::null),
        IfpackPrec(Teuchos::null),
        MLPrec(Teuchos::null)
    {
      int iflag;
      UpdateMatrix(A,sigma,B,&iflag);
      if (iflag!=0) throw iflag;
    }
    
    prec_and_mat::~prec_and_mat()
    {
    }

    void prec_and_mat::UpdateMatrix(phist_Dconst_sparseMat_ptr vA, double sigma,
               phist_Dconst_sparseMat_ptr vB, int *iflag)
    {
      *iflag=0;
      PHIST_CAST_PTR_FROM_VOID(const Epetra_CrsMatrix,A,vA,*iflag);
      const Epetra_CrsMatrix *B=(const Epetra_CrsMatrix*)vB;
      
    // note: the Ifpack class requires a non-const Epetra_RowMatrix*, so at this point we're
    // forced to copy the matrix even if sigma==0
    if (Mat==Teuchos::null)
    {
      Mat = Teuchos::rcp(new Epetra_CrsMatrix(*A));
    }
    else
    {
      *Mat=*A;
    }
    if (sigma!=0.0)
    {
      if (B!=NULL)
      {
        PHIST_CHK_IERR(*iflag=EpetraExt::MatrixMatrix::Add(*B,false,-sigma,*Mat,1.0),*iflag);
      }
      else
      {
        Teuchos::RCP<Epetra_Vector> diag=Teuchos::rcp(new Epetra_Vector(A->RowMap()));
        PHIST_CHK_IERR(*iflag=A->ExtractDiagonalCopy(*diag),*iflag);
        for (int i=0; i<diag->MyLength(); i++) (*diag)[i]-=sigma;
        PHIST_CHK_IERR(*iflag=Mat->ReplaceDiagonalValues(*diag),*iflag);
      }
    }

    }
  }
}    
