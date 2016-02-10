#include "phist_config.h"

#ifdef PHIST_HAVE_IFPACK

#ifndef PHIST_KERNEL_LIB_EPETRA
# error "Ifpack only works with Epetra kernel library!"
#endif

#include "phist_void_aliases.h"
#include "phist_trilinos_macros.h"
#include "phist_macros.h"

#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Ifpack.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

namespace phist {

class PreconTraits<double,IFPACK>
{

  public:

  static void Usage()
  {
    PHIST_SOUT(PHIST_INFO,"Ifpack: accepts name of an XML parameter file as 'options' string.\n"
                          "        The parameter list may contain two entries, \n"
                          "        \"Method\" (passed as PrecType ot the Ifpack factory)\n"
                          "        \"schwarz: overlap level\" (set to >0 to use overlapping Additive Schwarz)\n"
                          "The remaining parameter list is passed to the factory unchanged.\n");
  }

  static void Create(void** P, 
        const void* vA, double sigma, const void* vB, 
        Dconst_mvec_ptr_t Vkern, Dconst_mvec_ptr_t BVkern,
        std::string options, int* iflag)
  {
    PHIST_ENTER_FCN(__FUNCTION__);
    *iflag=0;
    PHIST_CAST_PTR_FROM_VOID(const Epetra_CrsMatrix, A, vA,*iflag);
    const Epetra_CrsMatrix* B = (const Epetra_CrsMatrix*)vB;
    
    Teuchos::RCP<Teuchos::ParameterList> ifpack_list=Teuchos::rcp(new Teuchos::ParameterList());
    
    PHIST_TRY_CATCH(updateParametersFromXmlFile(options,ifpack_list.ptr()),*iflag);
    
    Ifpack Factory;

    std::string PrecType = ifpack_list->get("Method","ILU");
    ifpack_list->remove("Method");
    int OverlapLevel = ifpack_list->get("schwarz: overlap level",0);
    ifpack_list->remove("schwarz: overlap level");
    
    // this parameter is (correctly) int in Ifpack2 but double in Ifpack, we want to
    // be more compatible so we convert it to double if it is found as an int
    if (ifpack_list->isType<int>("fact: ilut level-of-fill"))
    {
      int lof=ifpack_list->get<int>("fact: ilut level-of-fill");
      ifpack_list->remove("fact: ilut level-of-fill");
      ifpack_list->set("fact: ilut level-of-fill",(double)lof);
    }

    // computing A-sigma*B is possible in Epetra but not implemented here
    PHIST_CHK_IERR(*iflag= (sigma!=0.0)? -99:0,*iflag);
    
    // note: the Ifpack class requires a non-const Epetra_RowMatrix*, so technically at this point we're
    // forced to copy the matrix anyway. Just for testing things out we avoid this right now, we'd have
    // to decide who's responsible for deleting the copied matrix and when (TODO).

    Ifpack_Preconditioner* Prec = Factory.Create(PrecType, (Epetra_CrsMatrix*)A, OverlapLevel);
    PHIST_CHK_IERR(*iflag=Prec!=NULL?0:PHIST_BAD_CAST,*iflag);

    PHIST_CHK_IERR(*iflag=Prec->SetParameters(*ifpack_list),*iflag);
    PHIST_CHK_IERR(*iflag=Prec->Initialize(),*iflag);
    PHIST_CHK_IERR(*iflag=Prec->Compute(),*iflag);
    
    // return created object as void pointer
    *P=(void*)Prec;
    
    return;
  }

  static void Delete(void* vP, int *iflag)
  {
    PHIST_ENTER_FCN(__FUNCTION__);
    *iflag=0;
    PHIST_CAST_PTR_FROM_VOID(Ifpack_Preconditioner, P, vP,*iflag);
    delete P;
  }
  
  static void Apply(double alpha, void const* vP, Dconst_mvec_ptr_t vX, double beta, Dmvec_ptr_t vY, int* iflag)
  {
    PHIST_ENTER_FCN(__FUNCTION__);
    *iflag=0;
    PHIST_CAST_PTR_FROM_VOID(const Ifpack_Preconditioner, P, vP,*iflag);
    PHIST_CAST_PTR_FROM_VOID(const Epetra_MultiVector, X, vX,*iflag);
    PHIST_CAST_PTR_FROM_VOID(      Epetra_MultiVector, Y, vY,*iflag);
    if (alpha!=1.0||beta!=0.0)
    {
      PHIST_SOUT(PHIST_ERROR,"Ifpack preconditioner can only be applied as Y=inv(P)X up to now\n");
      *iflag=PHIST_NOT_IMPLEMENTED;
      return;
    }
    PHIST_CHK_IERR(*iflag=P->ApplyInverse(*X,*Y),*iflag);
  }
  
  static void ApplyT(double alpha, void const* P, Dconst_mvec_ptr_t X, double beta, Dmvec_ptr_t Y, int* iflag)
  {
    PHIST_ENTER_FCN(__FUNCTION__);
    // we currently don't need to apply the transpose of a preconditioner, and in Ifpack
    // it would mean we have to call SetUseTranspose on the operator, which is not possible
    // here as it is const. So we do the default thing, return -99 and leave the problem
    // for whoever stumbles on it first.
    *iflag=PHIST_NOT_IMPLEMENTED;
    return;
  }
  
  static void ApplyShifted(double alpha, const void* vP, double const * sigma,
          Dconst_mvec_ptr_t vX, double beta,  Dmvec_ptr_t vY, int* iflag)
  {
    // As we are talking about preconditioning here, we have the freedom to simply apply the same operator
    // to each column (ignoring the shift altogether or using a single shift). We could decide what to do 
    // depending on the Ifpack method used, for instance polynomial preconditioners can handle a different
    // shift for each column.
    PHIST_CHK_IERR(Apply(alpha,vP,vX,beta,vY,iflag),*iflag);
    
  }

};

} // namespace phist
#endif
