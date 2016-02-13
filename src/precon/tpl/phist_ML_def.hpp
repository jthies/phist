#include "phist_config.h"

#ifdef PHIST_HAVE_ML

#ifndef PHIST_KERNEL_LIB_EPETRA
# error "ML only works with Epetra kernel library!"
#endif

#include "phist_void_aliases.h"
#include "phist_trilinos_macros.h"
#include "phist_macros.h"

#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"

#include "ml_MultiLevelPreconditioner.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

namespace phist {

template<>
class PreconTraits<double,phist_ML>
{

  public:

  static void Usage()
  {
    PHIST_SOUT(PHIST_INFO,"ML: accepts name of an XML parameter file as 'options' string,\n"
                          "    with all the valid ML parameters and optionally a string  \n"
                          "    'default values' which may be 'SA', 'NSSA', 'DD' etc.     \n"
                          "    default values are only used for entries not present in the list.\n");
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
    
    Teuchos::RCP<Teuchos::ParameterList> ml_list=Teuchos::rcp(new Teuchos::ParameterList());
    
    PHIST_TRY_CATCH(updateParametersFromXmlFile(options,ml_list.ptr()),*iflag);
    if (ml_list->isParameter("default values"))
    {
      std::string defaults=ml_list->get("default values","SA");
      ML_Epetra::SetDefaults(defaults,*ml_list,NULL,NULL,false);
      ml_list->remove("default values");
    }
    
    // computing A-sigma*B is possible in Epetra but not implemented here
    PHIST_CHK_IERR(*iflag= (sigma!=0.0)? -99:0,*iflag);
    
    ML_Epetra::MultiLevelPreconditioner* Prec = new ML_Epetra::MultiLevelPreconditioner(*A, *ml_list);
    PHIST_CHK_IERR(*iflag=Prec!=NULL?0:PHIST_BAD_CAST,*iflag);
    
    // return created object as void pointer
    *P=(void*)Prec;
    
    return;
  }

  static void Delete(void* vP, int *iflag)
  {
    PHIST_ENTER_FCN(__FUNCTION__);
    *iflag=0;
    PHIST_CAST_PTR_FROM_VOID(ML_Epetra::MultiLevelPreconditioner, P, vP,*iflag);
    delete P;
  }
  
  static void Apply(double alpha, void const* vP, Dconst_mvec_ptr_t vX, double beta, Dmvec_ptr_t vY, int* iflag)
  {
    PHIST_ENTER_FCN(__FUNCTION__);
    *iflag=0;
    PHIST_CAST_PTR_FROM_VOID(const ML_Epetra::MultiLevelPreconditioner, P, vP,*iflag);
    PHIST_CAST_PTR_FROM_VOID(const Epetra_MultiVector, X, vX,*iflag);
    PHIST_CAST_PTR_FROM_VOID(      Epetra_MultiVector, Y, vY,*iflag);
    if (alpha!=1.0||beta!=0.0)
    {
      PHIST_SOUT(PHIST_ERROR,"ML preconditioner can only be applied as Y=inv(P)X up to now\n");
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
