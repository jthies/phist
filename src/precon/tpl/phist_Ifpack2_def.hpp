#include "phist_config.h"

#ifdef PHIST_HAVE_IFPACK2

#ifndef PHIST_KERNEL_LIB_TPETRA
# error "Ifpack only works with Tpetra kernel library!"
#endif

#include "phist_void_aliases.h"
#include "phist_trilinos_macros.h"
#include "phist_macros.h"

#include "phist_ScalarTraits.hpp"
#include "phist_tpetra_typedefs.hpp"

#include "Ifpack2_Factory.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

namespace phist {

template<typename ST>
class PreconTraits<ST,IFPACK>
{

  typedef Tpetra::MultiVector<ST,lidx_t,gidx_t,node_t> mvec_t;
  typedef Tpetra::CrsMatrix<ST,lidx_t,gidx_t,node_t> sparseMat_t;
  typedef Ifpack2::Preconditioner<ST,lidx_t,gidx_t,node_t> prec_type;

  public:

  static void Usage()
  {
    PHIST_SOUT(PHIST_INFO,"Ifpack2: accepts name of an XML parameter file as 'options' string.\n"
                          "         The parameter list may contain an entry \n"
                          "         \"Method\" (passed as PrecType ot the Ifpack2 factory).\n"
                          "The remaining parameter list is passed to the factory unchanged.\n");
  }

  static void Create(void** P, 
        const void* vA, ST sigma, const void* vB, 
        void const* Vkern, void const* BVkern,
        std::string options, int* iflag)
  {
    PHIST_ENTER_FCN(__FUNCTION__);
    *iflag=0;
    PHIST_CAST_PTR_FROM_VOID(const sparseMat_t, A, vA,*iflag);
    const sparseMat_t* B = (const sparseMat_t*)vB;
    
    Teuchos::RCP<Teuchos::ParameterList> ifpack_list=Teuchos::rcp(new Teuchos::ParameterList());
    
    PHIST_TRY_CATCH(updateParametersFromXmlFile(options,ifpack_list.ptr()),*iflag);
    
    Ifpack2::Factory Factory;

    std::string PrecType = ifpack_list->get("Method","ILU");
    ifpack_list->remove("Method");

    // computing A-sigma*B is possible in Epetra but not implemented here
    PHIST_CHK_IERR(*iflag= (sigma!=0.0)? -99:0,*iflag);
    
    Teuchos::RCP<const sparseMat_t> A_ptr = Teuchos::rcp(A,false);

    Teuchos::RCP<prec_type> Prec = Factory.create(PrecType, A_ptr);
    PHIST_CHK_IERR(*iflag=Prec!=Teuchos::null?0:PHIST_BAD_CAST,*iflag);

    PHIST_CHK_IERR(Prec->setParameters(*ifpack_list),*iflag);
    PHIST_CHK_IERR(Prec->initialize(),*iflag);
    PHIST_CHK_IERR(Prec->compute(),*iflag);
    
    // return created object as void pointer
    *P=(void*)Prec.release().get();
    
    return;
  }

  static void Delete(void* vP, int *iflag)
  {
    PHIST_ENTER_FCN(__FUNCTION__);
    *iflag=0;
    PHIST_CAST_PTR_FROM_VOID(prec_type, P, vP,*iflag);
    delete P;
  }
  
  static void Apply(ST alpha, void const* vP, Dconst_mvec_ptr_t vX, ST beta, Dmvec_ptr_t vY, int* iflag)
  {
    PHIST_ENTER_FCN(__FUNCTION__);
    *iflag=0;
    PHIST_CAST_PTR_FROM_VOID(const prec_type, P, vP,*iflag);
    PHIST_CAST_PTR_FROM_VOID(const mvec_t, X, vX,*iflag);
    PHIST_CAST_PTR_FROM_VOID(      mvec_t, Y, vY,*iflag);
    if (alpha!=1.0||beta!=0.0)
    {
      PHIST_SOUT(PHIST_ERROR,"Ifpack preconditioner can only be applied as Y=inv(P)X up to now\n");
      *iflag=PHIST_NOT_IMPLEMENTED;
      return;
    }
    PHIST_CHK_IERR(P->apply(*X,*Y),*iflag);
  }
  
  static void ApplyT(ST alpha, void const* P, Dconst_mvec_ptr_t X, ST beta, Dmvec_ptr_t Y, int* iflag)
  {
    PHIST_ENTER_FCN(__FUNCTION__);
    // we currently don't need to apply the transpose of a preconditioner, and in Ifpack
    // it would mean we have to call SetUseTranspose on the operator, which is not possible
    // here as it is const. So we do the default thing, return -99 and leave the problem
    // for whoever stumbles on it first.
    *iflag=PHIST_NOT_IMPLEMENTED;
    return;
  }
  
  static void ApplyShifted(ST alpha, const void* vP, ST const * sigma,
          Dconst_mvec_ptr_t vX, ST beta,  Dmvec_ptr_t vY, int* iflag)
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
