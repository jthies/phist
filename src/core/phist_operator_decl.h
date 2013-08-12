typedef struct TYPE(op) {
 const void* A_; // data structure needed for representing A
 const_map_ptr_t range_map_; //! map for vectors Y in Y=A*X
 const_map_ptr_t domain_map_; //! map for vectors X in Y=A*X
 // pointer to function for computing Y=alpha*A*X+beta*Y
 void (*apply)(_ST_ alpha, const void* A, 
        TYPE(const_mvec_ptr) X, _ST_ beta,  TYPE(mvec_ptr) Y, int* ierr);
} TYPE(op);

typedef TYPE(op)* TYPE(op_ptr);
typedef const TYPE(op)* TYPE(const_op_ptr);

//! this function can be used to create an operator which encapsulates a CRS matrix.
//! It does not allocate memory for the op struct, the caller has to do that beforehand.
void SUBR(op_wrap_crsMat)(TYPE(op_ptr) op, TYPE(const_crsMat_ptr) A, int* ierr);
