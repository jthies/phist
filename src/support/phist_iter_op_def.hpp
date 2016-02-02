typedef struct {
#ifndef IS_COMPLEX
  TYPE(feastCorrectionSolver) *fCS;
#endif
  int maxIters;
  _MT_ tol;
} TYPE(iter_solver_op);

extern "C" void SUBR(private_iter_op_apply)
        (_ST_ alpha, const void* arg, TYPE(const_mvec_ptr) X,
        _ST_ beta,  TYPE(mvec_ptr) Y, int* iflag)
{
#include "phist_std_typedefs.hpp"
#ifdef IS_COMPLEX
  // there is no complex feastCorrectionSolver right now, I think.
  *iflag=PHIST_NOT_IMPLEMENTED;
  return;
#else

  if (alpha!=st::one()||beta!=st::zero())
  {
    *iflag=PHIST_INVALID_INPUT;
    return;
  }
  PHIST_CAST_PTR_FROM_VOID(TYPE(iter_solver_op),op,arg,*iflag);
  PHIST_CHK_IERR(SUBR(feastCorrectionSolver_run)
        (op->fCS,X,op->tol,op->maxIters, &Y, NULL, iflag),*iflag);
  return;
#endif
}

//! wrap up an iterative solver as an operator, that
//! is op->apply(alpha,op,X,beta,Y) will actually 
//! approximate Y=(A-sigma*I)\X. alpha and beta are
//! not used, resp. must be alpha=1, beta=0. 
void SUBR(op_wrap_solver)(TYPE(linearOp_ptr) Ainv_op,TYPE(const_sparseMat_ptr) A, _ST_ shift,
        linSolv_t method,int block_size, _MT_ tol,int maxIter,int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
#ifdef IS_COMPLEX
  // there is no complex feastCorrectionSolver right now, I think.
  *iflag=PHIST_NOT_IMPLEMENTED;
  return;
#else

  //! create a feastCorrectionSolver object. You have to specify the number of rhs
  //! per shift to be treated simultaneously (blockSize), the number of shifts
  //! (numShifts), and the complex shifts
  TYPE(iter_solver_op) *op = new TYPE(iter_solver_op);
  Ainv_op->A = (void*)op;
  MT shift_r=st::real(shift);
  MT shift_i=st::imag(shift);
  PHIST_CHK_IERR(SUBR(feastCorrectionSolver_create)(&op->fCS, A, method,
        block_size, 1, &shift_r,&shift_i,iflag),*iflag);

// set remaining params and pointers

  op->maxIters=maxIter;
  op->tol=tol;

  //only defined for square matrices, range=domain map
  PHIST_CHK_IERR(SUBR(sparseMat_get_range_map)(A,&Ainv_op->range_map,iflag),*iflag);
Ainv_op->domain_map=Ainv_op->range_map;
Ainv_op->apply=&SUBR(private_iter_op_apply);
  return;    
#endif
  }
