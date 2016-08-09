//! a simple Lanczos process to compute the largest and smallest eigenvalue 
//! of a linear operator. This code is mainly intended as a testbed 
//! for fault tolerance (FT) capabilities in PHIST right now. Only block size 1 
//! is implemented.
//!
//! input:
//!
//! A: linearOp whose largest eigenpair is sought
//! *numIter, maximum number of iterations allowed
//!
//! output:
//!
//! *lambda_min, lambda_max: approximations to the largest and smallest eigenvalue
//! *numIter, number of iterations performed
//!
void SUBR(simple_lanczos_ft)(TYPE(const_linearOp_ptr) A_op,
        _MT_* lambda_min, _MT_* lambda_max, int *numIter, Cp_Options * cpOpt, int* iflag);



