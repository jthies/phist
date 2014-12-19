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

  
