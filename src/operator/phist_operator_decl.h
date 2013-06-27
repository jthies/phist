typedef struct _TYPE_(op) {
 const void* A_; // data structure needed for representing A
 const_map_ptr_t range_map_; //! map for vectors Y in Y=A*X
 const_map_ptr_t domain_map_; //! map for vectors X in Y=A*X
 // pointer to function for computing Y=alpha*A*X+beta*Y
 void (*apply)(_ST_ alpha, const void* A, 
        _TYPE_(const_mvec_ptr) X, _ST_ beta,  _TYPE_(mvec_ptr) Y, int* ierr);
} _TYPE_(op);

typedef _TYPE_(op)* _TYPE_(op_ptr);
typedef const _TYPE_(op)* _TYPE_(const_op_ptr);

void _SUBR_(Op_wrap_crsMat)(_TYPE_(op_ptr) op, _TYPE_(const_crsMat_ptr) A, int* ierr);
