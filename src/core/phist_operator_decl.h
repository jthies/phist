
//!\name basic operator concept, "something that only provides Y=alpha*Op*X + beta*Y"
//@{

//!
typedef struct TYPE(op) {
 const void* A; //! data structure needed for representing A
 const_map_ptr_t range_map; //! map for vectors Y in Y=A*X
 const_map_ptr_t domain_map; //! map for vectors X in Y=A*X
 //! pointer to function for computing Y=alpha*A*X+beta*Y
 void (*apply)(_ST_ alpha, const void* A, 
        TYPE(const_mvec_ptr) X, _ST_ beta,  TYPE(mvec_ptr) Y, int* ierr);
 //! pointer to function for computing Y=(A-sigma[j]B)*X[j]+beta*Y[j]
 void (*apply_shifted)(_ST_ alpha, const void* A, _ST_ const * sigma,
        TYPE(const_mvec_ptr) X, _ST_ beta,  TYPE(mvec_ptr) Y, int* ierr);
} TYPE(op);

typedef TYPE(op)* TYPE(op_ptr);
typedef const TYPE(op)* TYPE(const_op_ptr);

//! this function can be used to create an operator which encapsulates a CRS matrix.
//! It does not allocate memory for the op struct, the caller has to do that beforehand.
void SUBR(op_wrap_crsMat)(TYPE(op_ptr) op, TYPE(const_crsMat_ptr) A, int* ierr);

//! create the identity operator that returns Y=alpha*X+beta*Y
void SUBR(op_identity)(TYPE(op_ptr) op, int* ierr);

//! given a crsMat, create the operator I-DCSWP(A,x), double CARP sweep.
//! This operator can be passed to CG for solving Ax=b via the normal   
//! equations AA'y=b, x=A'y. We allow for an array of shifts here, so   
//! that the operator acting on each column i is in fact                
//! I-DCSWP(A-s[i]I,x[i]). If shift is NULL, s=0 is assumed.            
void SUBR(op_carp)(TYPE(op_ptr) op, TYPE(const_crsMat_ptr) A,
        _MT_ omega, int* ierr);

//@}


