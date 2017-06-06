/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_config.h"

#ifdef PHIST_HAVE_IFPACK

#ifndef PHIST_KERNEL_LIB_EPETRA
# error "Ifpack only works with Epetra kernel library!"
#endif

#include "phist_void_aliases.h"
#include "phist_trilinos_macros.h"
#include "phist_macros.h"
#include "../../kernels/epetra/epetra_helpers.h"

#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Ifpack.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

namespace phist {

template<>
class PreconTraits<double,phist_IFPACK>
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

  static void Wrap(void** P, 
        const void* vA, double sigma, const void* vB, 
        phist_Dconst_mvec_ptr Vkern, phist_Dconst_mvec_ptr BVkern,
        void* last_arg, int* iflag)
  {
    PHIST_ENTER_FCN(__FUNCTION__);
    *iflag=0;
    PHIST_CAST_PTR_FROM_VOID(const Epetra_CrsMatrix, A, vA,*iflag);
    const Epetra_CrsMatrix* B = (const Epetra_CrsMatrix*)vB;
    Ifpack_Preconditioner* P_in = (Ifpack_Preconditioner*)last_arg;
    
    phist::internal::prec_and_mat* PAM=new phist::internal::prec_and_mat(A,sigma,B);

    PAM->IfpackPrec = Teuchos::rcp(P_in,false);
    PAM->Prec=PAM->IfpackPrec;
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
    
    phist::internal::prec_and_mat* PAM=new phist::internal::prec_and_mat(A,sigma,B);

    PAM->IfpackPrec = Teuchos::rcp(Factory.Create(PrecType, PAM->Mat.get(), OverlapLevel));
    PAM->Prec=PAM->IfpackPrec;
    PHIST_CHK_IERR(*iflag=PAM->Prec.get()!=NULL?0:PHIST_BAD_CAST,*iflag);

    PHIST_CHK_IERR(*iflag=PAM->IfpackPrec->SetParameters(*ifpack_list),*iflag);
    PHIST_CHK_IERR(*iflag=PAM->IfpackPrec->Initialize(),*iflag);
    PHIST_CHK_IERR(*iflag=PAM->IfpackPrec->Compute(),*iflag);
    
    // return created object as void pointer
    *P=(void*)PAM;
    
    return;
  }

  static void Update(void* vP, const void* vA, double sigma, const void* vB,
        const void* Vkern, const void* BVkern,
        int* iflag)
  {
    PHIST_ENTER_FCN(__FUNCTION__);
    *iflag=0;
    PHIST_CAST_PTR_FROM_VOID(phist::internal::prec_and_mat, PAM, vP,*iflag);
    PHIST_CHK_IERR(PAM->UpdateMatrix(vA,sigma,vB,iflag),*iflag);
    PHIST_CHK_IERR(*iflag=PAM->IfpackPrec->Compute(),*iflag);
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
    PHIST_CHK_IERR(*iflag=PAM->IfpackPrec->ApplyInverse(*X,*Y2),*iflag);
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
