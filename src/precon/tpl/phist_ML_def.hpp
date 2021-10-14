/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_config.h"

#ifdef PHIST_HAVE_ML

#ifndef PHIST_KERNEL_LIB_EPETRA
# error "ML only works with Epetra kernel library!"
#endif

#ifndef DOXYGEN
#include "phist_void_aliases.h"
#endif //DOXYGEN

#include "phist_trilinos_macros.hpp"
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

  static void Wrap(void** P, 
        const void* vA, double sigma, const void* vB, 
        phist_Dconst_mvec_ptr Vkern, phist_Dconst_mvec_ptr BVkern,
        void* last_arg, int* iflag)
  {
    PHIST_ENTER_FCN(__FUNCTION__);
    *iflag=0;
    PHIST_CAST_PTR_FROM_VOID(const Epetra_CrsMatrix, A, vA,*iflag);
    const Epetra_CrsMatrix* B = (const Epetra_CrsMatrix*)vB;
    ML_Epetra::MultiLevelPreconditioner* P_in = (ML_Epetra::MultiLevelPreconditioner*)last_arg;
    
    phist::internal::prec_and_mat* PAM=new 
        phist::internal::prec_and_mat((phist_Dconst_sparseMat_ptr)A,sigma,(phist_Dconst_sparseMat_ptr)B);

    PAM->MLPrec = Teuchos::rcp(P_in,false);
    PAM->Prec=PAM->MLPrec;
    PHIST_CHK_IERR(*iflag=PAM->Prec.get()!=NULL?0:PHIST_BAD_CAST,*iflag);
    
    // return created object as void pointer
    *P=(void*)PAM;
    
    return;
  }

  static void Create(void** P, 
        const void* vA, double sigma, const void* vB, 
        phist_Dconst_mvec_ptr Vkern, phist_Dconst_mvec_ptr BVkern,
        std::string options, void* last_arg, int* iflag)
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

    phist::internal::prec_and_mat* PAM=new 
        phist::internal::prec_and_mat((phist_Dconst_sparseMat_ptr)A,sigma,(phist_Dconst_sparseMat_ptr)B);

    PAM->MLPrec = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*(PAM->Mat), *ml_list));
    PAM->Prec=PAM->MLPrec;
    PHIST_CHK_IERR(*iflag=PAM->Prec!=Teuchos::null?0:PHIST_BAD_CAST,*iflag);
    
    // return created object as void pointer
    *P=(void*)PAM;
    
    return;
  }

  static void Update(void* P, const void* A, double sigma, const void* B,
                     phist_Dconst_mvec_ptr Vkern, phist_Dconst_mvec_ptr BVkern,
                     int* iflag) 
  {
    PHIST_CAST_PTR_FROM_VOID(phist::internal::prec_and_mat,PAM,P,*iflag);
    PHIST_CHK_IERR(PAM->UpdateMatrix((phist_Dconst_sparseMat_ptr)A,sigma,(phist_Dconst_sparseMat_ptr)B,iflag),*iflag);
    int dimV=0;
    Epetra_MultiVector* V = (Epetra_MultiVector*)Vkern;
//    if (V) dimV=std::min(2,V->NumVectors());
    if (dimV==0)
    {
      // without null space:
      PHIST_CHK_IERR(*iflag=PAM->MLPrec->ComputePreconditioner(true),*iflag);
    }
    else
    {
      // with null space (adaptive SA)
      PHIST_CHK_IERR(*iflag=PAM->MLPrec->ComputeAdaptivePreconditioner(dimV,V->Values()),*iflag);
    }
  }                                                                             

  static void Delete(void* vP, int *iflag)
  {
    PHIST_ENTER_FCN(__FUNCTION__);
    *iflag=0;
    PHIST_CAST_PTR_FROM_VOID(phist::internal::prec_and_mat, PAM, vP,*iflag);
    delete PAM;
  }
  
  static void Apply(double alpha, void const* vP, phist_Dconst_mvec_ptr vX, double beta, phist_Dmvec_ptr vY, int* iflag)
  {
    PHIST_ENTER_FCN(__FUNCTION__);
    *iflag=0;
    PHIST_CAST_PTR_FROM_VOID(const phist::internal::prec_and_mat, PAM, vP,*iflag);
    PHIST_CAST_PTR_FROM_VOID(const Epetra_MultiVector, X, vX,*iflag);
    PHIST_CAST_PTR_FROM_VOID(      Epetra_MultiVector, Y, vY,*iflag);
    Teuchos::RCP<Epetra_MultiVector> Y2 = Teuchos::rcp(Y,false);
    if (beta!=0.0)
    {
      Y2=Teuchos::rcp(new Epetra_MultiVector(*Y));
    }
    PHIST_CHK_IERR(*iflag=PAM->Prec->ApplyInverse(*X,*Y2),*iflag);
    if (beta!=0.0)
    {
      Y->Update(alpha,*Y2,beta);
    }
    else if (alpha!=1.0)
    {
      Y->Scale(alpha);
    }
  }
  
  static void ApplyT(double alpha, void const* P, phist_Dconst_mvec_ptr X, double beta, phist_Dmvec_ptr Y, int* iflag)
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
          phist_Dconst_mvec_ptr vX, double beta,  phist_Dmvec_ptr vY, int* iflag)
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
