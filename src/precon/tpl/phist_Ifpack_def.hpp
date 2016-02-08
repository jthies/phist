#include "phist_config.h"
#include "phist_PreconTraits.hpp"

#ifdef PHIST_HAVE_IFPACK

#ifndef PHIST_KERNEL_LIB_EPETRA
# error "Ifpack only works with Epetra kernel library!"
#endif

#include "phist_macros.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Ifpack.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

typedef struct {
Teuchos::RCP<Ifpack_Preconditioner> P;
} phist_IfpackWrapper;

class PreconTraits<double,IFPACK>
{

  static void Usage()
  {
    PHIST_SOUT(PHIST_INFO,"Ifpack: accepts name of an XML parameter file as 'options' string.\n"
                          "        The parameter list may contain two entries, \n"
                          "        \"Method\" (passed as PrecType ot the Ifpack factory)\n"
                          "        \"Overlap\" (set to >0 to use overlapping Additive Schwarz)\n"
                          "The remaining parameter list is passed to the factory unchanged.\n");
);
  }

  static void Create(DlinearOp_t* P, 
        const void* vA, double sigma, const void* vB, std::string options, int* iflag)
  {
    PHIST_ENTER_FCN(__FUNCTION__);
    *iflag=0;
    CAST_PTR_FROM_VOID(const Epetra_CrsMatrix, A, vA,iflag);
    const Epetra_CrsMatrix* B = (const Epetra_CrsMatrix*)vB;
    
    Teuchos::RCP<Teuchos::ParameterList> ifpack_list=Teuchos::rcp(new Teuchos::ParameterList());
    
    Teuchos::Comm<int> comm=TODO;
    TRY_CATCH(updateParametersFromXmlFileAndBroadcast(options,*ifpack_list.ptr(),comm),*iflag);
    
    Ifpack Factory;

    std::string PrecType = ifpack_list->get("Method","ILU");
    ifpack_list->remove("Method");
    int OverlapLevel = ifpack_list->get("Overlap",0);
    ifpack_list->remove("Overlap");

    // computing A-sigma*B is possible in Epetra but not implemented here
    PHIST_CHK_IERR(*iflag= (sigma!=0.0)? -99:0,*iflag);

    Ifpack_Preconditioner* Prec = Factory.Create(PrecType, A, OverlapLevel);
    PHIST_CHK_IERR(iflag=Prec!=Teuchos::null?0:PHIST_BAD_CAST,*iflag);

    IFPACK_CHK_ERR(Prec->SetParameters(*ifpack_list));
    IFPACK_CHK_ERR(Prec->Initialize());
    IFPACK_CHK_ERR(Prec->Compute());    
    
    return iflag;
  }

  

};

#endif
