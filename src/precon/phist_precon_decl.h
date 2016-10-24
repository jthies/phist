
//! \defgroup precon methods to construct, update and destroy a preconditioner
//@{

//! create a preconditioner for an iterative linear solver

//! this function can be used to create an operator that can be used to precondition linear systems     
//! with A-sigma*B. If sigma is 0, B is not touched. If sigma!=0 but B==NULL, B=I (identity matrix)     
//! is assumed.                                                                                         
//!
//! Applying the preconditioner is done via the usual apply member functions of the linearOp struct.    
//!
//! For some preconditioners it may be useful to provide (an approximation of) the kernel of            
//! A-sigma*B, this can be done via Vkern and BVkern. If they are NULL, they are not used anyway.       
//!
//! The preconditioning method is selected based on the string <method>. Options are                    
//! passed via the options string. The methods available depend on the kernel lib used and the          
//! optional libraries found while building phist. Passing in method="usage". The option string         
//! in turn depends on the input <method>. Passing in a supported method and options="usage" prints     
//! usage information for that particular preconditioner.                                               
//!                                                                                                     
//! The special string method="update" can be used if a preconditioner has already                      
//! been created and only the values (but not the pattern) of A, B and sigma changed.                   
//! method="update shift" indicates that only the value of sigma changed. In both cases, the options    
//! string is ignored.                                                                                  
//!
//! Example:                                                                                            
//!                                                                                                     
//! If you use the Epetra kernel library and the Trilinos library Ifpack (incomplete factorization      
//! preconditioners) is available, you could use                                                        
//!                                                                                                     
//! DlinearOp_t P;                                                                                      
//! phist_Dprecon_create(&P, &A, 0, NULL, NULL, NULL, "ifpack", "ifpack_params.xml", &iflag);           
//!                                                                                                     
//! which will read the preconditioner settings from an XML file compatible with the Teuchos            
//! ParameterList.                                                                                      
//!                                                                                                     
void SUBR(precon_create)(TYPE(linearOp_ptr) op, TYPE(const_sparseMat_ptr) A, 
                         _ST_ sigma, TYPE(const_sparseMat_ptr) B,
                         TYPE(const_mvec_ptr) Vkern, TYPE(const_mvec_ptr) BVkern,
                         const char* method, const char* options, int* iflag);

//! given an existing preconditioner, recompute it for a new shift sigma and (near) kernel Vkern.

//! The matrices A and B and preconditioning options remain unchanged from the cann SUBR(precon_create).
//! This means that the preconditioner needs to store pointers to A and B at create() time and can assume
//! those matrices are still there and unchanged. If they aren't (there or unchanged), a new preconditioner
//! must be created instead.
void SUBR(precon_update)(TYPE(linearOp_ptr) op, _ST_ sigma,
                         TYPE(const_mvec_ptr) Vkern,
                         TYPE(const_mvec_ptr) BVkern,
                         int* iflag);

//! destroy preconditioner
void SUBR(precon_delete)(TYPE(linearOp_ptr) op, int* iflag);

//! apply preconditioner to an mvec

//! This should be done via the member of the linearOp, not by calling this subroutine directly.
void SUBR(precon_apply)(_ST_ alpha, void const* P, TYPE(const_mvec_ptr) X, _ST_ beta, TYPE(mvec_ptr) Y, int* iflag);

//! apply transposed preconditioner to an mvec

//! This should be done via the member of the linearOp, not by calling this subroutine directly.
void SUBR(precon_applyT)(_ST_ alpha, void const* P, TYPE(const_mvec_ptr) X, _ST_ beta, TYPE(mvec_ptr) Y, int* iflag);

//! apply preconditioner to an mvec.

//! The shifts indicate that the operator preconditioned is A-sigma[j]B for column j of the mvec. It
//! is up to the preconditioner wether he exploits this info somehow.
//! This should be done via the member of the linearOp, not by calling this subroutine directly.
void SUBR(precon_apply_shifted)(_ST_ alpha, void const* P, _ST_ const* sigma,
        TYPE(const_mvec_ptr) X, 
        _ST_ beta, TYPE(mvec_ptr) Y, int* iflag);

//@}



