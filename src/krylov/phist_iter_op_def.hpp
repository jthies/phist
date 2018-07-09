/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/


typedef struct {
  TYPE(linearOp_ptr) A_op, // the operator
  TYPE(linearOp_ptr) P_op, // optional right preconditioner (only for certain iterative schemes)
  phist_ElinSolv method; // selects the iterative scheme to be used
  int maxIters; // number of iterations to perform
  _MT_ tol, // required residual
} TYPE(iter_solver_op);

extern "C" void SUBR(private_iter_op_apply)
        (_ST_ alpha, const void* arg, TYPE(const_mvec_ptr) X,
        _ST_ beta,  TYPE(mvec_ptr) Y, int* iflag)
{
#include "phist_std_typedefs.hpp"

  if (alpha!=st::one()||beta!=st::zero())
  {
    *iflag=PHIST_INVALID_INPUT;
    return;
  }
  PHIST_CAST_PTR_FROM_VOID(TYPE(iter_solver_op),op,arg,*iflag);

  *iflag=0;
  int nIter=op->maxIters;
  _MT_ tol=_MT_(0);

  if (op->method==phist_PCG)
  {
    PHIST_CHK_NEG_IERR(SUBR(PCG)(op->A_op, op->P_op, X, Y,
                                nIter, tol, iflag), *iflag);
  }
  else if (op->method==phist_GMRES)
  {
    *iflag=PHIST_NOT_IMPLEMENTED;
  /*
    PHIST_CHK_NEG_IERR(SUBR(restartedBlockedGMRES)(op->A_op, op->P_op, X, Y,
        int numSys,
        int* nIter, _MT_ const tol[], int block_size, int max_blocks, int* iflag);
  */
  }
  else if (op->method==phist_BICGSTAB)
  {
    PHIST_CHK_NEG_IERR(SUBR(BiCGStab)(op->A_op, op->P_op,
                                X, Y, &nIter, tol, iflag), *iflag);
  }
  else if (op->method==phist_CARP_CG)
  {
    PHIST_CHK_NEG_IERR(SUBR(carp_cg)( op->A_op->A, X, Y,
        &nIter, tol, iflag), *iflag);
  }
  else
  {
    PHIST_CHK_IERR(*iflag=PHIST_INVALID_INPUT,*iflag);
  }
  return;
}

//! wrap up an iterative solver as an operator, that
//! is op->apply(alpha,op,X,beta,Y) will actually 
//! approximate Y=(A-sigma*I)\X. alpha and beta are
//! not used, resp. must be alpha=1, beta=0. 
void SUBR(linearOp_wrap_solver)(TYPE(linearOp_ptr) Ainv_op,TYPE(const_sparseMat_ptr) A, _ST_ shift,
        phist_ElinSolv method,int block_size, _MT_ tol, int numIter,int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);

  //! create a feastCorrectionSolver object. You have to specify the number of rhs
  //! per shift to be treated simultaneously (blockSize), the number of shifts
  //! (numShifts), and the complex shifts
  TYPE(iter_solver_op) *op = new TYPE(iter_solver_op);
  Ainv_op->A = (void*)op;

  op->method=method;
  op->maxIters=maxIter;
  op->tol=tol;

  //only defined for square matrices, range=domain map
  PHIST_CHK_IERR(SUBR(sparseMat_get_range_map)(A,&Ainv_op->range_map,iflag),*iflag);
Ainv_op->domain_map=Ainv_op->range_map;
Ainv_op->apply=&SUBR(private_iter_op_apply);
  return;    
  }
