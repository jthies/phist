/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_config.h"

#ifdef PHIST_HAVE_IFPACK2

#ifndef PHIST_KERNEL_LIB_TPETRA
# error "Ifpack2 only works with Tpetra kernel library!"
#endif

#ifndef DOXYGEN
#include "phist_void_aliases.h"
#endif //DOXYGEN

#include "phist_trilinos_macros.hpp"
#include "phist_macros.h"

#include "phist_ScalarTraits.hpp"
#include "phist_tpetra_typedefs.hpp"

#include "Ifpack2_Factory.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "phist_DisallowPreconTraits.hpp"

#include "Teuchos_config.h"
#include "TpetraCore_config.h"

namespace phist {

/* do not allow instantiating Ifpac2 objects with data types not supported by the 
   Trilinos instaallation.
 */
#if !defined(HAVE_TEUCHOS_FLOAT)||!defined(HAVE_TPETRA_INST_FLOAT)
PHIST_DISALLOW_PRECON_TRAITS(float,phist_IFPACK)
#endif
#if !defined(HAVE_TEUCHOS_FLOAT)||!defined(HAVE_TEUCHOS_COMPLEX)||!defined(HAVE_TPETRA_INST_COMPLEX_FLOAT)
PHIST_DISALLOW_PRECON_TRAITS(phist_s_complex,phist_IFPACK)
#endif
#if !defined(HAVE_TEUCHOS_COMPLEX)||!defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE)
PHIST_DISALLOW_PRECON_TRAITS(phist_d_complex,phist_IFPACK)
#endif

template<typename ST>
class PreconTraits<ST,phist_IFPACK>
{

  typedef typename tpetra::Traits<ST>::mvec_t mvec_type;
  typedef typename tpetra::Traits<ST>::sparseMat_t sparseMat_type;
  typedef typename sparseMat_type::node_type node_type;
  typedef typename phist::ScalarTraits<ST>::mvec_t* phist_mvec_ptr;
  typedef typename phist::ScalarTraits<ST>::mvec_t const* phist_const_mvec_ptr;
  typedef Ifpack2::Preconditioner<ST,phist_lidx,phist_gidx,node_type> prec_type;
  typedef phist::ScalarTraits<ST> st;

  public:

  static void Usage()
  {
    PHIST_SOUT(PHIST_INFO,"Ifpack2: accepts name of an XML parameter file as 'options' string.\n"
                          "         The parameter list may contain an entry \n"
                          "         \"Method\" (passed as PrecType ot the Ifpack2 factory).\n"
                          "The remaining parameter list is passed to the factory unchanged.\n");
  }

  static void Wrap(void** P, 
        const void* vA, ST sigma, const void* vB, 
        void const* Vkern, void const* BVkern,
        void* last_arg, int* iflag)
  {
    PHIST_ENTER_FCN(__FUNCTION__);
    *iflag=PHIST_NOT_IMPLEMENTED;
  }

  static void Create(void** P, 
        const void* vA, ST sigma, const void* vB, 
        void const* Vkern, void const* BVkern,
        std::string options, void* last_arg, int* iflag)
  {
    PHIST_ENTER_FCN(__FUNCTION__);
    *iflag=0;
    PHIST_CAST_PTR_FROM_VOID(const sparseMat_type, A, vA,*iflag);
    const sparseMat_type* B = (const sparseMat_type*)vB;
    
    Teuchos::RCP<Teuchos::ParameterList> ifpack_list=Teuchos::rcp(new Teuchos::ParameterList());
    
    PHIST_TRY_CATCH(updateParametersFromXmlFile(options,ifpack_list.ptr()),*iflag);
    
    Ifpack2::Factory Factory;

    std::string PrecType = ifpack_list->get("Method","ILU");
    ifpack_list->remove("Method");

    // computing A-sigma*B is possible in Epetra but not implemented here
    PHIST_CHK_IERR(*iflag= (sigma!=st::zero())? -99:0,*iflag);
    
    Teuchos::RCP<const sparseMat_type> A_ptr = Teuchos::rcp(A,false);

    Teuchos::RCP<prec_type> Prec = Factory.create(PrecType, A_ptr);
    PHIST_CHK_IERR(*iflag=Prec!=Teuchos::null?0:PHIST_BAD_CAST,*iflag);

    PHIST_CHK_IERR(Prec->setParameters(*ifpack_list),*iflag);
    PHIST_CHK_IERR(Prec->initialize(),*iflag);
    PHIST_CHK_IERR(Prec->compute(),*iflag);
    
    // return created object as void pointer
    *P=(void*)Prec.release().get();
    
    return;
  }

  static void Update(void* P, const void* A, ST sigma, const void* B,
        const void* Vkern, const void* BVkern,
        int* iflag)
  {
    PHIST_CHK_IERR(*iflag=PHIST_NOT_IMPLEMENTED,*iflag);
  }                                                                             

  static void Delete(void* vP, int *iflag)
  {
    PHIST_ENTER_FCN(__FUNCTION__);
    *iflag=0;
    PHIST_CAST_PTR_FROM_VOID(prec_type, P, vP,*iflag);
    delete P;
  }
  
  static void Apply(ST alpha, void const* vP, phist_const_mvec_ptr vX, ST beta, phist_mvec_ptr vY, int* iflag)
  {
    PHIST_ENTER_FCN(__FUNCTION__);
    *iflag=0;
    PHIST_CAST_PTR_FROM_VOID(const prec_type, P, vP,*iflag);
    PHIST_CAST_PTR_FROM_VOID(const mvec_type, X, vX,*iflag);
    PHIST_CAST_PTR_FROM_VOID(      mvec_type, Y, vY,*iflag);

    PHIST_CHK_IERR(P->apply(*X,*Y,Teuchos::NO_TRANS,alpha,beta),*iflag);
  }
  
  static void ApplyT(ST alpha, void const* vP, phist_const_mvec_ptr vX, ST beta, phist_mvec_ptr vY, int* iflag)
  {
    PHIST_ENTER_FCN(__FUNCTION__);
    *iflag=0;
    PHIST_CAST_PTR_FROM_VOID(const prec_type, P, vP,*iflag);
    PHIST_CAST_PTR_FROM_VOID(const mvec_type, X, vX,*iflag);
    PHIST_CAST_PTR_FROM_VOID(      mvec_type, Y, vY,*iflag);

    PHIST_CHK_IERR(P->apply(*X,*Y,Teuchos::TRANS,alpha,beta),*iflag);
  }
  
  static void ApplyShifted(ST alpha, const void* vP, ST const * sigma,
          phist_const_mvec_ptr vX, ST beta,  phist_mvec_ptr vY, int* iflag)
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
