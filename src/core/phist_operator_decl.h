
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
} TYPE(op);

typedef TYPE(op)* TYPE(op_ptr);
typedef const TYPE(op)* TYPE(const_op_ptr);

//! this function can be used to create an operator which encapsulates a CRS matrix.
//! It does not allocate memory for the op struct, the caller has to do that beforehand.
void SUBR(op_wrap_crsMat)(TYPE(op_ptr) op, TYPE(const_crsMat_ptr) A, int* ierr);

//! create the identity operator that returns Y=alpha*X+beta*Y
void SUBR(op_identity)(TYPE(op_ptr) op, int* ierr);

//@}


//! This function allows us to put the operation Y=Op*X in a task buffer and thus perform it
//! on an array of vectors without loading the matrix k times into cache. Note that we do   
//! not support the parameters alpha and beta that op->apply accepts, we assume alpha=1 and 
//! beta=0 for now. The in_args of the argList should be the vectors X, the out_args the    
//! vectors Y. The shared_arg is the operator itself. X and Y are allowed to contain 
//! pointers to non-contiguous memory, in which case the vectors are copied before and after
//! the operation.
void SUBR(op_apply_multi)(argList_t* args);
