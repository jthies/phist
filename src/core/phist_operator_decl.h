
//! \defgroup core Core functionality used by algorithms
//@{

//!\name basic operator concept, "something that only provides Y=alpha*Op*X + beta*Y"
//@{

//!
typedef struct TYPE(linearOp) {
 const void* A; //! data structure needed for representing A
 phist_const_map_ptr range_map; //! map for vectors Y in Y=A*X
 phist_const_map_ptr domain_map; //! map for vectors X in Y=A*X
 void const* aux; //! This field can be used to carry along
                  //! additional info like a space which
                  //! is projected out of the operator etc.,
                  //! it is currently ignored in phist but
                  //! used in HYMLS to implement a custom
                  //! residual evaluation. In the future we
                  //! want to provide some mechanism to add
                  //! user-defined projections
 //! pointer to function for computing Y=alpha*A*X+beta*Y
 void (*apply)(_ST_ alpha, const void* A, 
        TYPE(const_mvec_ptr) X, _ST_ beta,  TYPE(mvec_ptr) Y, int* iflag);
//! apply transpose
 void (*applyT)(_ST_ alpha, const void* A, 
        TYPE(const_mvec_ptr) X, _ST_ beta,  TYPE(mvec_ptr) Y, int* iflag);
 //! pointer to function for computing Y=(A-sigma[j]B)*X[j]+beta*Y[j]
 void (*apply_shifted)(_ST_ alpha, const void* A, _ST_ const * sigma,
        TYPE(const_mvec_ptr) X, _ST_ beta,  TYPE(mvec_ptr) Y, int* iflag);
} TYPE(linearOp);

typedef TYPE(linearOp)* TYPE(linearOp_ptr);
typedef const TYPE(linearOp)* TYPE(const_linearOp_ptr);

//! this function can be used to create an operator which encapsulates a CRS matrix.
//! It does not allocate memory for the op struct, the caller has to do that beforehand.
void SUBR(linearOp_wrap_sparseMat)(TYPE(linearOp_ptr) op, TYPE(const_sparseMat_ptr) A, int* iflag);

//! create the identity operator that returns Y=alpha*X+beta*Y
void SUBR(linearOp_identity)(TYPE(linearOp_ptr) op, int* iflag);

//@}

//@}


