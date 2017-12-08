/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/

//! \defgroup core Core functionality used by algorithms
//@{

//!\name basic operator concept, "something that only provides Y=alpha*Op*X + beta*Y"
//@{

//!
typedef struct TYPE(linearOp) {
 void const* A; //! data structure needed for representing A
 void * aux;    //! non-const data field that can e.g. be used
                //! for storing preconditioner data to be updated
                //! when calling update(). To implement this one
                //! can e.g. have aux point to A (providing a non-
                //! const path to the preconditioner to the update
                //! function).
 phist_const_map_ptr range_map; //! map for vectors Y in Y=A*X
 phist_const_map_ptr domain_map; //! map for vectors X in Y=A*X
 //! pointer to function for computing Y=alpha*A*X+beta*Y
 void (*apply)(_ST_ alpha, const void* A, 
        TYPE(const_mvec_ptr) X, _ST_ beta,  TYPE(mvec_ptr) Y, int* iflag);
//! apply transpose
 void (*applyT)(_ST_ alpha, const void* A, 
        TYPE(const_mvec_ptr) X, _ST_ beta,  TYPE(mvec_ptr) Y, int* iflag);
 //! pointer to function for computing Y=(A-sigma[j]B)*X[j]+beta*Y[j]
 void (*apply_shifted)(_ST_ alpha, const void* A, _ST_ const * sigma,
        TYPE(const_mvec_ptr) X, _ST_ beta,  TYPE(mvec_ptr) Y, int* iflag);
  //! apply operator and compute inner products with in- and output vector
  void (*fused_apply_mvTmv)(_ST_ alpha, const void* A, TYPE(const_mvec_ptr)  V,
                            _ST_ beta,                 TYPE(mvec_ptr)        W,
                            TYPE(sdMat_ptr) WtW, TYPE(sdMat_ptr) VtW,
                            int* iflag);
  // given an existing operator, update it for a new shift sigma and (near) kernel Vkern.
  // this function is mainly intended for implementing custom preconditioners, before actually
  // calling it on any linearOp you should check if it is not NULL.
  void (*update)(const void* A, void* aux, _ST_ sigma,
                        TYPE(const_mvec_ptr) Vkern,
                        TYPE(const_mvec_ptr) BVkern,
                        int* iflag);
  //! this function can be used to clean up any data the operator may *own*,
  //! if the operator is just a wrapper for some other object that is created
  //! and deleted separately, this function should not do anything.
  //! The me object itself should *not* be free'd.
  void (*destroy)(struct TYPE(linearOp)* me, int* iflag);
  
} TYPE(linearOp);

typedef TYPE(linearOp)* TYPE(linearOp_ptr);
typedef const TYPE(linearOp)* TYPE(const_linearOp_ptr);

//! this function can be used to create an operator which encapsulates a CRS matrix.
//! It does not allocate memory for the op struct, the caller has to do that beforehand.
void SUBR(linearOp_wrap_sparseMat)(TYPE(linearOp_ptr) op, TYPE(const_sparseMat_ptr) A, int* iflag);

//! given two sparse matrices A and B, this operator acts as Y=alpha*A*X+beta*Y (apply)
//! or Y=alpha*(AX+sigma[j]B)X_j + beta*Y (apply_shifted)
void SUBR(linearOp_wrap_sparseMat_pair)(TYPE(linearOp_ptr) op, 
                                        TYPE(const_sparseMat_ptr) A, TYPE(const_sparseMat_ptr) B, 
                                        int* iflag);

//! create a 'product operator' from opA and opB which acts as Y <- alpha*A*B*X + beta*Y (apply).
//! The other functions will at the moment return an error (-99, not implemented) because they are
//! typically not needed and their behavior has not been defined yet for product operators.
//! An exception is 'destroy', which will delete temporary storage used by the operator but not the
//! operators A and B, which are only wrapped and not owned.
//!
//! CAVEAT: as the resulting wrapped operator cannot access details of the implementation of A and B,
//!         we cannot use a fused kernel here. For optimal performance e.g. if A and B are sparse row
//!         matrices, a fused kernel should be used.
//!
void SUBR(linearOp_wrap_linearOp_product)(TYPE(linearOp_ptr) op,
        TYPE(const_linearOp_ptr) A, TYPE(const_linearOp_ptr) B, int* iflag);
		
//! create a 'product triple operator' from opA, opB and opC which acts as Y <- alpha*A*B*C*X + beta*Y (apply).
//! The other functions will at the moment return an error (-99, not implemented) because they are
//! typically not needed and their behavior has not been defined yet for product operators.
//! An exception is 'destroy', which will delete temporary storage used by the operator but not the
//! operators A, B and C, which are only wrapped and not owned.
//!
//! CAVEAT: as the resulting wrapped operator cannot access details of the implementation of A, B and C,
//!         we cannot use a fused kernel here. For optimal performance e.g. if A, B and C are sparse row
//!         matrices, a fused kernel should be used.
//!
void SUBR(linearOp_wrap_linearOp_product_triple)(TYPE(linearOp_ptr) op,
        TYPE(const_linearOp_ptr) A, TYPE(const_linearOp_ptr) B, TYPE(const_linearOp_ptr) C, int* iflag);
		
//! create a 'product k operator' which acts as Y = alpha*A1_apply1*...*Ak_applyk*X + beta*Y (apply).
//! k_ops gives us the k operators A1,...,Ak, (k_ops[i]=Ai the i-th operator) and 
//! which_apply tells us, how we want to apply the operators (which_apply[i]=0 if apply i-th operator,
//! which_apply[i]=1 if applyT and which_apply[i]=2 if apply_shifted)
//! The other functions will at the moment return an error (-99, not implemented) because they are
//! typically not needed and their behavior has not been defined yet for product operators.
//! An exception is 'destroy', which will delete temporary storage used by the operator but not the
//! operators A, B and C, which are only wrapped and not owned.
//!
//! CAVEAT: as the resulting wrapped operator cannot access details of the implementation of A, B and C,
//!         we cannot use a fused kernel here. For optimal performance e.g. if A, B and C are sparse row
//!         matrices, a fused kernel should be used.
//!		
void SUBR(linearOp_wrap_linearOp_product_k)(TYPE(linearOp_ptr) op,
int k, TYPE(const_linearOp_ptr)* k_ops, int* which_apply, const _ST_* sigma, int num_shifts, int* iflag);		
		

#if defined(__cplusplus)&&defined(PHIST_KERNEL_LIB_EPETRA)&&defined(IS_DOUBLE)&&!defined(IS_COMPLEX)
// forward declaration
class Epetra_Operator;

//! this function can be used to create an operator which encapsulates an Epetra_Operator.
//! It does not allocate memory for the op struct, the caller has to do that beforehand.
void SUBR(linearOp_wrap_epetra)(TYPE(linearOp_ptr) op, Epetra_Operator const* A, int* iflag);
#endif

//! create the identity operator that returns Y=alpha*X+beta*Y
void SUBR(linearOp_identity)(TYPE(linearOp_ptr) op, 
                             phist_const_map_ptr  range_map,
                             phist_const_map_ptr domain_map, int* iflag);

//@}

//! \name wrappers that simply call the corresponding function in the linearOp struct
//! These are particularly useful for Fortran users, for which it is awkward to use  
//! the c_funptr members of the linearOp types (cf. https://bitbucket.org/essex/phist_fort)
//@{
 //! pointer to function for computing Y=alpha*A*X+beta*Y
 void SUBR(linearOp_apply)(_ST_ alpha, TYPE(const_linearOp_ptr) A_op, 
        TYPE(const_mvec_ptr) X, _ST_ beta,  TYPE(mvec_ptr) Y, int* iflag);
//! apply transpose
 void SUBR(linearOp_applyT)(_ST_ alpha, TYPE(const_linearOp_ptr) A_op, 
        TYPE(const_mvec_ptr) X, _ST_ beta,  TYPE(mvec_ptr) Y, int* iflag);
 //! pointer to function for computing Y=(A-sigma[j]B)*X[j]+beta*Y[j]
 void SUBR(linearOp_apply_shifted)(_ST_ alpha, TYPE(const_linearOp_ptr) A_op, _ST_ const * sigma,
        TYPE(const_mvec_ptr) X, _ST_ beta,  TYPE(mvec_ptr) Y, int* iflag);
  //! apply operator and compute inner products with in- and output vector
  void SUBR(linearOp_fused_apply_mvTmv)(_ST_ alpha, TYPE(const_linearOp_ptr) A_op, TYPE(const_mvec_ptr)  V,
                            _ST_ beta,                 TYPE(mvec_ptr)        W,
                            TYPE(sdMat_ptr) WtW, TYPE(sdMat_ptr) VtW,
                            int* iflag);
  
  //! this function can be used to clean up any data the operator may *own*,
  //! if the operator is just a wrapper for some other object that is created
  //! and deleted separately, this function should not do anything.
  //! The me object itself should *not* be free'd.
  void SUBR(linearOp_destroy)(TYPE(linearOp_ptr) A_op, int* iflag);

  //!
  void SUBR(linearOp_update)(TYPE(linearOp_ptr) A_op, _ST_ sigma,
                        TYPE(const_mvec_ptr) Vkern,
                        TYPE(const_mvec_ptr) BVkern,
                        int* iflag);

//@}

//@}


