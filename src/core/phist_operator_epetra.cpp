#include "phist_config.h"

#ifdef PHIST_KERNEL_LIB_EPETRA

#include "phist_tools.h"
#include "phist_kernels.h"
#include "phist_operator.h"
#include "Epetra_Operator.h"
#include "phist_gen_d.h"

#include "Epetra_Operator.h"
#include "Epetra_MultiVector.h"

namespace phist {
namespace internal {

  //! pointer to function for computing Y=alpha*A*X+beta*Y
  void Epetra_Operator_apply(_ST_ alpha, const void* vOp,
        TYPE(const_mvec_ptr) vX, _ST_ beta,  TYPE(mvec_ptr) vY, int* iflag)
  {
    PHIST_CAST_PTR_FROM_VOID(const Epetra_Operator,Op,vOp,*iflag);
    PHIST_CAST_PTR_FROM_VOID(const Epetra_MultiVector,X,vX,*iflag);
    PHIST_CAST_PTR_FROM_VOID(      Epetra_MultiVector,Y,vY,*iflag);
    int iflag1=0,iflag2=0;
    if (beta==0.0)
    {
      iflag1=Op->Apply(*X,*Y);
      if (alpha!=1.0) iflag2=Y->Scale(alpha);
    }
    else
    {
      Epetra_MultiVector AX(X->Map(),X->NumVectors());
      iflag1=Op->Apply(*X,AX);
      iflag2=Y->Update(alpha,AX,beta);
    }
    PHIST_CHK_IERR(*iflag=iflag1,*iflag);
    PHIST_CHK_IERR(*iflag=iflag2,*iflag);
  }

  //! apply transpose
  void Epetra_Operator_applyT(_ST_ alpha, const void* vOp,
        TYPE(const_mvec_ptr) X, _ST_ beta,  TYPE(mvec_ptr) Y, int* iflag)
  {
    // note: Epetra_Operator has a construction SetUseTranspose(true)/Apply/SetUseTranspose(false), but
    // we can't use it here because Op is const, so we just return an error for now.
    PHIST_CHK_IERR(*iflag=PHIST_NOT_IMPLEMENTED,*iflag);
  }
 //! pointer to function for computing Y=(A-sigma[j]B)*X[j]+beta*Y[j]
 void Epetra_Operator_apply_shifted(_ST_ alpha, const void* vOp, _ST_ const * sigma,
        TYPE(const_mvec_ptr) vX, _ST_ beta,  TYPE(mvec_ptr) vY, int* iflag)
  {
    PHIST_CAST_PTR_FROM_VOID(const Epetra_Operator,Op,vOp,*iflag);
    PHIST_CAST_PTR_FROM_VOID(const Epetra_MultiVector,X,vX,*iflag);
    PHIST_CAST_PTR_FROM_VOID(      Epetra_MultiVector,Y,vY,*iflag);
    int iflag1=0,iflag2=0,iflag3=0;

    _ST_ minus_shifts[X->NumVectors()];
    for (int i=0; i<X->NumVectors(); i++) minus_shifts[i]=-sigma[i];

    if (beta==0.0)
    {
      iflag1=Op->Apply(*X,*Y);
      SUBR(mvec_vadd_mvec)(minus_shifts, vX,1.0,vY,&iflag2);
      if (alpha!=1.0) iflag3=Y->Scale(alpha);
    }
    else
    {
      Epetra_MultiVector AX(X->Map(),X->NumVectors());
      iflag1=Op->Apply(*X,AX);
      SUBR(mvec_vadd_mvec)(minus_shifts, vX,1.0,&AX,&iflag2);
      iflag3=Y->Update(alpha,AX,beta);
    }

    PHIST_CHK_IERR(*iflag=iflag1,*iflag);
    PHIST_CHK_IERR(*iflag=iflag2,*iflag);
    PHIST_CHK_IERR(*iflag=iflag3,*iflag);
  }

//! apply operator and compute inner products with in- and output vector
  void Epetra_Operator_fused_apply_mvTmv(_ST_ alpha,   const void* vOp, TYPE(const_mvec_ptr)  V,
                            _ST_ beta,                 TYPE(mvec_ptr)        W,
                            TYPE(sdMat_ptr) WtW, TYPE(sdMat_ptr) VtW,
                            int* iflag)
  {
    // not implemented yet
    PHIST_CHK_IERR(*iflag=PHIST_NOT_IMPLEMENTED,*iflag);
  }

}
}
// this function can be used to create an operator which encapsulates an Epetra_Operator.
// It does not allocate memory for the op struct, the caller has to do that beforehand.
extern "C" void SUBR(linearOp_wrap_epetra)(TYPE(linearOp_ptr) op, Epetra_Operator const* A, int* iflag)
{
  *iflag=0;
  op->A = A;
  op->aux=NULL;
  op->range_map=&(A->OperatorRangeMap());
  op->domain_map=&(A->OperatorRangeMap());

  op->apply = &phist::internal::Epetra_Operator_apply;
  op->applyT = &phist::internal::Epetra_Operator_applyT;
  op->apply_shifted = &phist::internal::Epetra_Operator_apply_shifted;
  op->fused_apply_mvTmv = &phist::internal::Epetra_Operator_fused_apply_mvTmv;
  op->update=NULL;
  op->destroy=NULL;
  return;
}



#endif /* Epetra */
